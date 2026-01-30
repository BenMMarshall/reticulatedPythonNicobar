
library(here)
library(dplyr)
library(sf)
library(CoordinateCleaner)
library(geodata)
library(ggplot2)
library(terra)

nicoGrids <- read_sf(here("data", "Nicobar 3x3 grids.gpkg"))
plot(nicoGrids)
nicoIslands <- read_sf(here("data", "Nicobar_all_islands_UTM.gpkg"))
plot(nicoIslands)

pythonLocs <- read.csv(here("data", "sixty_Pythons_locations_UTM.csv")) %>%
  rename("x" = X, "y" = Y)

pythonSF <- st_as_sf(pythonLocs,
         coords = c("x", "y"),
         crs = 32648,
         remove = FALSE) %>%
  st_transform(4326) %>%
  st_coordinates() %>%
  as.data.frame() %>%
  mutate(species = "Malayopython reticulatus",
         source = "orig60locs")


occGBIFdata <- read.csv(here("data", "gbif_malayopythonreticulatus.csv"),
                        sep = "\t")


seed <- 2026


institutionCode

occData <- rbind(occGBIFdata %>%
                   ### FILTER OUT INACCURATE LOCATIONS
                   filter(coordinateUncertaintyInMeters <= 250) %>%
                   mutate(source = paste0("GBIF_", institutionCode)) %>%
                   dplyr::select(species, decimalLatitude, decimalLongitude, source),
                 pythonSF %>%
                   dplyr::select(species,
                                 decimalLatitude = Y,
                                 decimalLongitude = X,
                                 source))


ggplot(occData) +
  geom_point(aes(x = decimalLongitude, y = decimalLatitude, colour = source)) +
  coord_equal()

# Clean coords ------------------------------------------------------------

flags <- clean_coordinates(x = occData,
                           lon = "decimalLongitude",
                           lat = "decimalLatitude",
                           species = "species",
                           tests = c("capitals", "centroids",
                                     "equal", "zeros", "gbif",
                                     "outliers"),
                           capitals_rad = 10000,
                           centroids_rad = 1000,
                           outliers_mtp = 5,
                           outliers_method = "quantile",
                           outliers_size = 7)

flags[!flags$.summary,]
summary(flags)
plot(flags, lon = "decimalLongitude", lat = "decimalLatitude")

occData <- occData[flags$.summary,]

occDataSF <- st_as_sf(occData,
                      coords = c("decimalLongitude", "decimalLatitude"),
                      crs = 4326,
                      remove = FALSE)


occExt <- ext(c(range(occDataSF$decimalLongitude)+c(-5,5), range(occDataSF$decimalLatitude)+c(-5,5)))

# Read in raster data -----------------------------------------------------
dir.create(here::here("data", "worldclim"))
worldclim_global(var = "bio", res = 10, path = here::here("data", "worldclim"))

climFiles <- list.files(here::here("data", "worldclim", "climate", "wc2.1_10m"))

climRast <- rast(here::here("data", "worldclim", "climate", "wc2.1_10m", climFiles))

climCrop <- crop(climRast, occExt)


# Build bias layer --------------------------------------------------------

hfRast <- rast(here::here("data", "humanFootprint", "hfp2022.tif"))
hfCrop <- project(hfRast, climCrop)

rescale <- function(x){(x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))}
hfDataBNG <- hfDataBNG %>%
  mutate(hfp2022 = rescale(hfp2022))

weightedRandom <- spatSample(# x = hfData25m,
  x = hfMasked,
  # size = 10,
  size = nrow(occData)*nPointMultiplier*nReps,
  method = "weights",
  na.rm = TRUE, as.df = FALSE, values = FALSE,
  xy = TRUE)

# Generate pseudo-absences ------------------------------------------------

#### NEEDS TEARING APART
# pseudoAbs <- create_psuedo_abs(occDataSF,
#                                hfBiasLayer = here("data", "humanFootprint", "hfp2022.tif"),
#                                nPointMultiplier = 3, nReps = 10,
#                                envLayers = read_stack_layers(layerLoc = here("data", "rasterLayers"),
#                                                              tar_sdm_layers))

# Format data for model ---------------------------------------------------

biomodData <- BIOMOD_FormatingData(resp.var = pseudoAbs$sp,
                                   expl.var = pseudoAbs$env,
                                   resp.xy = pseudoAbs$xy,
                                   resp.name = "Malayopython.reticulatus",
                                   #     # advice from biomod2’s team:
                                   # # - random selection of PA when high specificity is valued over high sensitivity
                                   # # - number of PA = 3 times the number of presences
                                   # # - 10 repetitions
                                   # PA.strategy = "random",
                                   # PA.nb.rep = 2,
                                   # PA.nb.absences = 1000,
                                   filter.raster = TRUE,
                                   # PA.nb.rep = 2,
                                   PA.strategy = "user.defined",
                                   PA.user.table = pseudoAbs$pa.tab
)


# Run suite of models -----------------------------------------------------

biomodModels <- BIOMOD_Modeling(bm.format = biomodData,
                                modeling.id = "AllModels",
                                models = c("ANN",
                                           "GBM",
                                           "GLM",
                                           # "MAXNET",
                                           # "XGBOOST",
                                           "RF"
                                ),
                                CV.strategy = "user.defined",
                                CV.user.table = bm_CrossValidation(bm.format = biomodData,
                                                                   strategy = "strat",
                                                                   k = 4,
                                                                   balance = "presences",
                                                                   strat = "both"),
                                # CV.strategy = "random",
                                # CV.nb.rep = 2,
                                # CV.perc = 0.8,
                                OPT.strategy = "bigboss",
                                var.import = 3,
                                metric.eval = c("TSS","ROC"),
                                seed.val = seed,
                                nb.cpu = 6)


# Run with best models combined -------------------------------------------


biomodEns <- BIOMOD_EnsembleModeling(bm.mod = biomodModels,
                                     models.chosen = "all",
                                     em.by = "all",
                                     em.algo = c("EMmean", #"EMcv", "EMci",
                                                 "EMmedian", #"EMca",
                                                 "EMwmean"),
                                     metric.select = c("TSS"),
                                     metric.select.thresh = c(0.25),
                                     metric.eval = c("TSS", "ROC"),
                                     var.import = 3,
                                     EMci.alpha = 0.05,
                                     EMwmean.decay = "proportional",
                                     seed.val = seed,
                                     nb.cpu = 6)

# Project models to area of interest --------------------------------------

biomodForecast <- BIOMOD_EnsembleForecasting(bm.em = biomodModels,
                                             proj.name = "CurrentEM_",
                                             ########### CROP STUFF TO NICOBAR
                                             new.env = read_stack_layers(layerLoc = here("data", "GIS data", "SDM Layers"),
                                                                         tar_sdm_layers) %>%
                                               crop(tar_patchList$Wessex),
                                             ########### CROP STUFF TO NICOBAR
                                             models.chosen = get_built_models(biomodEns)[
                                               stringr::str_detect(get_built_models(biomodEns), "EMwmeanByTSS")],
                                             metric.binary = "all",
                                             metric.filter = "all",
                                             output.format = ".tif",
                                             nb.cpu = 6)

# Save final suitability raster -------------------------------------------

save_proj_layer(tar_biomodForecast_fallow)
