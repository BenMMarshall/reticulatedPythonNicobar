
library(here)
library(dplyr)
library(sf)
library(CoordinateCleaner)
library(geodata)
library(ggplot2)
library(ggrepel)
library(terra)
library(tidyterra)
library(biomod2)
library(stringr)

seed <- 2026

nicoIslands <- read_sf(here("data", "Nicobar_all_islands_UTM.gpkg")) %>%
  st_transform(4326)

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
worldclim_global(var = "bio", res = 0.5, path = here::here("data", "worldclim"))

climFiles <- list.files(here::here("data", "worldclim", "climate", "wc2.1_10m"))
# climFiles <- list.files(here::here("data", "worldclim", "climate", "wc2.1_30s"))
climRast <- rast(here::here("data", "worldclim", "climate", "wc2.1_10m", climFiles))
# climRast <- rast(here::here("data", "worldclim", "climate", "wc2.1_30s", climFiles))

climCrop <- crop(climRast, occExt)

# Build bias layer --------------------------------------------------------

hfRast <- rast(here::here("data", "humanFootprint", "hfp2022.tif"))
hfCrop <- project(hfRast, climCrop)

rescale <- function(x){(x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))}
hfDataBNG <- hfCrop %>%
  mutate(hfp2022 = rescale(hfp2022))

nPointMultiplier = 3; nReps = 3

occDataSimple <- occData %>%
  mutate(x = decimalLongitude,
         y = decimalLatitude) %>%
  mutate(resp = 1) %>%
  select(x, y, resp)

weightedRandom <- spatSample(# x = hfData25m,
  x = hfDataBNG,
  # size = 10,
  size = nrow(occData)*nPointMultiplier*nReps,
  method = "weights",
  na.rm = TRUE, as.df = FALSE, values = FALSE,
  xy = TRUE)

fullRespData <- rbind(occDataSimple, weightedRandom %>%
                        mutate(resp = 0))

PAtable <- as.data.frame(matrix(FALSE, nrow = nrow(fullRespData), ncol = nReps))
names(PAtable) <- paste0("PA", 1:nReps)
# known locations always included
PAtable[which(fullRespData$resp == 1),] <- TRUE

# Generate pseudo-absences ------------------------------------------------

for (i in 1:ncol(PAtable)) PAtable[sample(which(PAtable[, i] == FALSE),
                                          nrow(occData)*nPointMultiplier), i] = TRUE

fullRespData <- as_spatvector(fullRespData, geom = c("x", "y")) %>%
  select(resp)

print("- bm_PseudoAbsences start")

pseudoAbs <- bm_PseudoAbsences(resp.var = fullRespData,
                               expl.var = climCrop,
                               strategy = "user.defined",
                               user.table = PAtable)

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
                                           "MAXNET",
                                           "XGBOOST",
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
                                     metric.select.thresh = c(0.5),
                                     metric.eval = c("TSS", "ROC"),
                                     var.import = 3,
                                     EMci.alpha = 0.05,
                                     EMwmean.decay = "proportional",
                                     seed.val = seed,
                                     nb.cpu = 6)

# Project models to area of interest --------------------------------------

biomodForecast <- BIOMOD_EnsembleForecasting(bm.em = biomodEns,
                                             proj.name = "CurrentEM_",
                                             ########### CROP STUFF TO NICOBAR
                                             new.env = climCrop %>%
                                               crop(st_bbox(nicoIslands)+c(-5,-5,5,5)),
                                             ########### CROP STUFF TO NICOBAR
                                             models.chosen = get_built_models(biomodEns)[
                                               stringr::str_detect(get_built_models(biomodEns), "EMwmeanByTSS")],
                                             metric.binary = "all",
                                             metric.filter = "all",
                                             output.format = ".tif",
                                             nb.cpu = 6)

# Save final suitability raster -------------------------------------------

forecastTerra <- rast(here::here(biomodForecast@sp.name,
                                 paste0("proj_", biomodForecast@proj.name),
                                 paste0("proj_", biomodForecast@proj.name, "_", biomodForecast@sp.name, "_ensemble.tif")))

plot(forecastTerra)


# Output plots ------------------------------------------------------------
## can also be done for single models as well as the ensemble

species <- biomodEns@sp.name

ggplotThemeCombo <-
  theme_bw() +
  theme(
    text = element_text(colour = "grey25"),
    line = element_line(colour = "grey25"),
    plot.title = element_text(face = 2),
    axis.title = element_text(face = 2),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_text(angle = 0, face = 4, hjust = 0, vjust = 1)
    # legend.position = "none"
  )

# Single models -----------------------------------------------------------

evalModels <- bm_PlotEvalMean(bm.out = biomodModels)

evalPlot <- evalModels$tab %>%
  ggplot() +
  geom_errorbar(aes(x = mean1, ymin = mean2-sd2, ymax = mean2+sd2, colour = name),
                width = 0.001)+
  geom_errorbarh(aes(y = mean2, xmin = mean1-sd1, xmax = mean1+sd1, colour = name),
                 height = 0.001)+
  geom_point(aes(x = mean1, y = mean2, colour = name), size = 3) +
  labs(y = "TSS", x = "ROC", colour = "Model") +
  ggplotThemeCombo

evalPlot

varImportByModel <- bm_PlotVarImpBoxplot(bm.out = biomodModels, group.by = c('expl.var', 'algo', 'algo'))

varImportPlot <- varImportByModel$tab %>%
  mutate(expl.var = reorder(expl.var, var.imp, FUN = max)) %>%
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.85) +
  geom_boxplot(aes(y = expl.var, x = var.imp, colour = algo, fill = algo)) +
  facet_grid(cols = vars(algo)) +
  ggplotThemeCombo +
  theme(legend.position = "none")

varImportPlot

respCurve <- bm_PlotResponseCurves(bm.out = biomodModels,
                                   models.chosen = get_built_models(biomodModels),
                                   fixed.var = 'mean')

respDistPlot <- respCurve$tab %>%
  mutate(model.type = str_extract(pred.name, "[^_]+$")) %>%
  ggplot() +
  # geom_rect(data = wessexRanges, aes(xmin = min, xmax = max,
  #                                    ymin = -Inf, ymax = Inf), fill = "#85AB7A", alpha = 0.2) +
  geom_path(aes(x = expl.val, y = pred.val, group = pred.name, colour = model.type)) +
  facet_wrap(facet = vars(expl.name), scales = "free") +
  ggplotThemeCombo

varImportPlot

# Ensemble models ---------------------------------------------------------

evalEns <- bm_PlotEvalMean(bm.out = biomodEns, group.by = 'full.name')

evalEnsPlot <- evalEns$tab %>%
  mutate(name = str_extract_all(name, "EM.*?ByTSS", simplify = TRUE)) %>%
  ggplot() +
  geom_point(aes(x = mean1, y = mean2, colour = name, fill = name),
             size = 4, pch = 21) +
  geom_text_repel(aes(x = mean1, y = mean2, colour = name, label = name),
                  hjust = 0, box.padding = unit(1.2, "lines")) +
  labs(y = "TSS", x = "ROC") +
  ggplotThemeCombo +
  theme(legend.position = "none")

evalEnsPlot

varImport <- bm_PlotVarImpBoxplot(bm.out = biomodEns,
                                  group.by = c('expl.var', 'algo', 'merged.by.run'))

varImportPlot <- varImport$tab %>%
  filter(algo %in% c("EMwmean", "EMmean", "EMmedian", "EMca")) %>%
  mutate(expl.var = reorder(expl.var, var.imp, FUN = max)) %>%
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.85) +
  # geom_hline(yintercept = seq(1.5,11.5,1), linetype = "dashed", alpha = 0.25) +
  geom_boxplot(aes(y = expl.var, x = var.imp, colour = algo, fill = algo)) +
  scale_colour_manual(values = c("#DCBD0A",
                                 "#CD602A",
                                 "#9F7E93",
                                 "#85AB7A")) +
  scale_fill_manual(values = c("#DCBD0A",
                               "#CD602A",
                               "#9F7E93",
                               "#85AB7A")) +
  ggplotThemeCombo +
  theme_bw() +
  theme(
    # panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom")

varImportPlot

respCurve <- bm_PlotResponseCurves(bm.out = biomodEns,
                                   models.chosen = get_built_models(biomodEns)[
                                     stringr::str_detect(get_built_models(biomodEns),
                                                         "EMwmeanByTSS")],
                                   fixed.var = 'mean')

respDistPlot <- respCurve$tab %>%
  ggplot() +
  geom_path(aes(x = expl.val, y = pred.val, group = expl.name)) +
  facet_wrap(facet = vars(expl.name), scales = "free") +
  ggplotThemeCombo

respDistPlot
