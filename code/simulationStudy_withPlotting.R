######################
## Simulation study ##
######################
# This simulation study is in parts based on the simulation study of
# Meyer and Pebesma (2021): 
# - https://doi.org/10.1111/2041-210X.13650
# - https://github.com/HannaMeyer/MEE_AOA

rm(list=ls())

# data handling
library(CAST)
library(CASTvis)
library(virtualspecies)
library(caret)
library(terra)
library(sf)

# visualization
library(viridis)
library(see)
library(gridExtra)
library(stringr)
library(tidyterra)
library(ggplot2)
library(plotly)
library(reshape2)
library(ggpubr)
library(cowplot)
library(knitr)
library(xtable)

# rmse helper function
rmse <- function(pred,obs){sqrt( mean((pred - obs)^2, na.rm = TRUE) )}

# ploting styles for ggplot
element_textbox <- function(...) {
  el <- element_text(...)
  class(el) <- c("element_textbox", class(el))
  el
}
element_grob.element_textbox <- function(element, ...) {
  text_grob <- NextMethod()
  rect_grob <- element_grob(calc_element("strip.background", theme_gray()))

  ggplot2:::absoluteGrob(
    grid::gList(
      element_grob(calc_element("strip.background",theme_gray())),
      text_grob
    ),
    height = grid::grobHeight(text_grob),
    width = grid::unit(1, "npc")
  )
}
# plotting theme for ggplot
theme = theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  panel.grid = element_blank(),
  legend.title = element_blank(),
  plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)),
  strip.text = element_text(size = 14),
  legend.text = element_text(size = 14),
)

# Presettings for the simulation study
# design: random, clustered or biasedWithOutlier (all calculated)
npoints <- 100                                                                  # npoints: number of training samples
nclusters <- 10                                                                 # nclusters: number of clusters if design==clustered
maxdist <- 0.8                                                                  # maxdist: maxdist for clustered samples if design==clustered
countries <- c("Germany","Ireland","France", "Sweden")                          # countries: main sample countries if design==biased
countriesOutlier <- "Norway"                                              # countriesOutlier: outlier country if design==biasedWithOutlier (a single point is set here)
meansPCA <- c(3, -1)                                                            # meansPCA: means of the gaussian response functions to the 2 axes
sdPCA <- c(2, 2)                                                                # sdPCA: sds of the gaussian response functions to the 2 axes
simulateResponse <- c("bio_2","bio_5","bio_10", "bio_13", "bio_14","bio_19")    # simulateResponse: variables used to simulate the response
studyarea <- c(-15, 65, 30, 75)                                                 # studyarea: extent of study area. Default: Europe c(-15, 65, 30, 75), South America: c(-82,-34,-56,13)
seed <- 20                                                                      # seed: seed for case study reproduction



####################
## Get predictors ##
####################

# Download Worldclim data
geodata::worldclim_global(var = 'bio', path = 'data/raw/', res = 10)

# Load predictors
predictors <- rast(
  list (
    rast("./data/raw/wc2.1_10m/wc2.1_10m_bio_1.tif"),
    rast("./data/raw/wc2.1_10m/wc2.1_10m_bio_2.tif"),
    rast("./data/raw/wc2.1_10m/wc2.1_10m_bio_3.tif"),
    rast("./data/raw/wc2.1_10m/wc2.1_10m_bio_4.tif"),
    rast("./data/raw/wc2.1_10m/wc2.1_10m_bio_5.tif"),
    rast("./data/raw/wc2.1_10m/wc2.1_10m_bio_6.tif"),
    rast("./data/raw/wc2.1_10m/wc2.1_10m_bio_7.tif"),
    rast("./data/raw/wc2.1_10m/wc2.1_10m_bio_8.tif"),
    rast("./data/raw/wc2.1_10m/wc2.1_10m_bio_9.tif"),
    rast("./data/raw/wc2.1_10m/wc2.1_10m_bio_10.tif"),
    rast("./data/raw/wc2.1_10m/wc2.1_10m_bio_11.tif"),
    rast("./data/raw/wc2.1_10m/wc2.1_10m_bio_12.tif"),
    rast("./data/raw/wc2.1_10m/wc2.1_10m_bio_13.tif"),
    rast("./data/raw/wc2.1_10m/wc2.1_10m_bio_14.tif"),
    rast("./data/raw/wc2.1_10m/wc2.1_10m_bio_15.tif"),
    rast("./data/raw/wc2.1_10m/wc2.1_10m_bio_16.tif"),
    rast("./data/raw/wc2.1_10m/wc2.1_10m_bio_17.tif"),
    rast("./data/raw/wc2.1_10m/wc2.1_10m_bio_18.tif"),
    rast("./data/raw/wc2.1_10m/wc2.1_10m_bio_19.tif")
  )
)

# crop predictors to target area
aoi <- ext(studyarea, xy = FALSE)
predictors <- crop(predictors,aoi)
names(predictors) = names(predictors) |> str_remove(pattern = "wc2.1_10m_")
predictor_names <- names(predictors)

# create a mask for the land area:
mask <- predictors[[1]]
mask[!is.na(mask)] <- 1
names(mask) <- "mask"


#######################
## Generate Response ##
#######################
# The virtual response variable is created based on the PCA of a subset of Worldclim predictors.
# See the virtualspecies package for further information.

response_vs <- generateSpFromPCA(predictors[[simulateResponse]],
                                 means = meansPCA,sds = sdPCA, plot=F)

response <- response_vs$suitab.raster
names(response) <- "response"


#################################
## Simulate sampling locations ##
#################################

# transform mask
maskVect <- as.polygons(mask, aggregate = TRUE)
mask <- as.polygons(mask, aggregate = FALSE)
mask <- st_as_sf(mask)
mask <- st_make_valid(mask)
modeldomain <- st_transform(mask, "+init=epsg:4326")
modeldomain_vect <- maskVect

# random
set.seed(seed)
samplepoints_random <- st_as_sf(st_sample(mask, size = npoints, "random"))
samplepoints_random <- st_transform(samplepoints_random, "+init=epsg:4326")
st_write(samplepoints_random, "./data/samples/samplepoints_sim_study_random.geojson", driver = "GeoJSON", append = FALSE, delete_dsn = TRUE)

# clustered
set.seed(9352)
mask <- st_transform(mask, "+init=epsg:3857")
samplepoints_clustered <- clustered_sample(mask, npoints, nclusters,radius=200000)
samplepoints_clustered <- st_transform(samplepoints_clustered, "+init=epsg:4326")
mask <- st_transform(mask, "+init=epsg:4326")
st_write(samplepoints_clustered, "./data/samples/samplepoints_sim_study_clustered.geojson", driver = "GeoJSON", append = FALSE, delete_dsn = TRUE)

# biased with outlier
set.seed(7)
countryboundaries <- geodata::world(res = 4, path = './data/countryboundaries/')
countryboundariesOut <- st_as_sf(countryboundaries[countryboundaries$NAME_0%in%c(countriesOutlier),])
countryboundaries <- st_as_sf(countryboundaries[countryboundaries$NAME_0%in%c(countries),])
samplepoints_biased <- st_as_sf(st_sample(countryboundaries,npoints-1,"random"))
samplepoints_out <- st_as_sf(st_sample(countryboundariesOut,1,"random"))
samplepoints_biasedWithOutlier <- rbind(samplepoints_biased, samplepoints_out)
samplepoints_biasedWithOutlier <- st_transform(samplepoints_biasedWithOutlier, st_crs(mask))
st_write(samplepoints_biasedWithOutlier, "./data/samples/samplepoints_sim_study_biasedWithOutlier.geojson", driver = "GeoJSON", append = FALSE, delete_dsn = TRUE)


##############################################
##  Plot all predictors and sampling designs #
##############################################
# add response to predictors for plotting
predictors$response <- response

plot0 = ggplot() +
  geom_spatraster(data = stretch(predictors,0,1)) +
  facet_wrap(~lyr, ncol = 5) +
  scale_fill_viridis(na.value = "transparent") +
  theme +
  theme(legend.key.height = unit(0.17, "npc"))

plot1 = ggplot() +
  geom_spatvector(data = maskVect, show.legend = F, fill = "grey") +
  labs(title = "Random") +
  geom_sf(data = samplepoints_random, color = "red", size = 0.5) +
  theme

plot2 = ggplot() +
  geom_spatvector(data = maskVect, show.legend = F, fill = "grey") +
  labs(title = "Clustered") +
  geom_sf(data = samplepoints_clustered, color = "red", size = 0.5) +
  theme

plot3 = ggplot() +
  geom_spatvector(data = maskVect, show.legend = F, fill = "grey") +
  labs(title = "Biased with outlier") +
  geom_sf(data = samplepoints_biasedWithOutlier, color = "red", size = 0.5) +
  theme

plot4 = grid.arrange(plot1,plot2,plot3, nrow=3)

# plot only sampling designs
pdf(file = "code/figures/simulation_study/sample_designs.pdf", width = 15, height = 5)
plot_grid(plot1, plot2, plot3, labels = c("(a)", "(b)", "(c)"), ncol = 3)
invisible(dev.off())

# plot predictors and sampling designs
pdf(file = "code/figures/simulation_study/predictors_and_sample_designs.pdf", width = 12, height = 8)
plot_grid(plot0, plot4, ncol = 2, nrow = 1,rel_heights = c(0.8), rel_widths = c(735,265), labels = c("(a)", "(b)"))
invisible(dev.off())


###########################################
## Plot only PCA predictors and response ##
###########################################
plot1 = ggplot() +
  geom_spatraster(data = stretch(predictors[[simulateResponse]],0,1), show.legend = F) +
  facet_wrap(~lyr, ncol = 3) +
  scale_fill_viridis(na.value = "transparent") +
  labs(title = "Predictors") +
  theme

plot2 = ggplot() +
  geom_spatraster(data = response) +
  scale_fill_viridis(na.value = "transparent") +
  labs(title = "Response") +
  theme +
  theme(legend.key.height = unit(0.16, "npc")
  )

pdf(file = "code/figures/simulation_study/generate_response.pdf", width = 14, height = 6)
plot_grid(plot1, plot2, nrow = 1,rel_widths = c(530,470), rel_heights = c(1), labels = c("(a)", "(b)"))
invisible(dev.off())


# remove response from predictos
predictors <- subset(predictors, predictor_names)



################################
## Prepare training data sets ##
################################
# extract predictors for the sampling locations

# random
trainDat_random <- extract(predictors, samplepoints_random, df=TRUE)
trainDat_random$response <- extract(response, samplepoints_random, ID = FALSE)$response
trainDat_random <- trainDat_random[complete.cases(trainDat_random),]

# clustered
trainDat_clustered <- extract(predictors, samplepoints_clustered, df=TRUE)
trainDat_clustered$response <- extract(response, samplepoints_clustered, ID = FALSE)$response
trainDat_clustered$clstrID <- samplepoints_clustered$parent
trainDat_clustered <- trainDat_clustered[complete.cases(trainDat_clustered),]

# biased with outlier
trainDat_biasedWithOutlier <- extract(predictors, samplepoints_biasedWithOutlier, df=TRUE)
trainDat_biasedWithOutlier$response <- extract(response, samplepoints_biasedWithOutlier, ID = FALSE)$response
trainDat_biasedWithOutlier <- trainDat_biasedWithOutlier[complete.cases(trainDat_biasedWithOutlier),]



########################
## Create kNNDM folds ##
########################

set.seed(seed)

# random
knndm_folds_random <- knndm(samplepoints_random, modeldomain = modeldomain, k = 10)

samplepoints_random$fold = knndm_folds_random$clusters

plot_sample_cv_design_random = ggplot() +
  geom_spatvector(data = modeldomain_vect, fill = "white") +
  labs(title = "CV design (random)") +
  geom_sf(mapping = aes(col = factor(fold)), data = samplepoints_random, size = 1) +
  scale_color_manual(name = "knndm folds", values = c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695')) +
  theme +
  theme(legend.position.inside=c(.8,.16),
        legend.text=element_text(size = 8),
        legend.title=element_text(size = 10),
        legend.background=element_blank(),
        legend.key.height = unit(0.03, "npc"),
        legend.key.width = unit(0.015, "npc"))

# clustered
knndm_folds_clustered <- knndm(samplepoints_clustered, modeldomain = modeldomain, k = 10)

samplepoints_clustered$fold = knndm_folds_clustered$clusters

plot_sample_cv_design_clustered = ggplot() +
  geom_spatvector(data = modeldomain_vect, fill = "white") +
  labs(title = "CV design (clustered)") +
  geom_sf(mapping = aes(col = factor(fold)), data = samplepoints_clustered, size = 1) +
  scale_color_manual(name = "knndm folds", values = c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695')) +
  theme +
  theme(legend.position.inside=c(.8,.16),
        legend.text=element_text(size = 8),
        legend.title=element_text(size = 10),
        legend.background=element_blank(),
        legend.key.height = unit(0.03, "npc"),
        legend.key.width = unit(0.015, "npc"))

# biased with outlier
knndm_folds_biasedWithOutlier <- knndm(samplepoints_biasedWithOutlier, modeldomain = modeldomain, k = 10)

samplepoints_biasedWithOutlier$fold = knndm_folds_biasedWithOutlier$clusters

plot_sample_cv_design_biasedWithOutlier = ggplot() +
  geom_spatvector(data = modeldomain_vect, fill = "white") +
  labs(title = "CV design (biased with outliers)") +
  geom_sf(mapping = aes(col = factor(fold)), data = samplepoints_biasedWithOutlier, size = 1) +
  scale_color_manual(name = "knndm folds", values = c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695')) +
  theme +
  theme(legend.position.inside=c(.8,.16),
        legend.text=element_text(size = 8),
        legend.title=element_text(size = 10),
        legend.background=element_blank(),
        legend.key.height = unit(0.03, "npc"),
        legend.key.width = unit(0.015, "npc"))


pdf(file = "code/figures/simulation_study/sample_designs_knndm_folds.pdf", width = 15, height = 5)
plot_grid(plot_sample_cv_design_random, plot_sample_cv_design_clustered, plot_sample_cv_design_biasedWithOutlier, labels = c("(a)", "(b)", "(c)"), ncol = 3)
invisible(dev.off())


####################
## Model training ##
####################
# Model training is then performed using the caret package.
# The model output gives information on the general estimated
# model performance based on the cross-validation.

set.seed(seed)

# random with 10-fold random CV
model_random <- train(
  trainDat_random[, names(predictors)],
  trainDat_random$response,
  method = "rf",
  importance = TRUE,
  tuneGrid = expand.grid(mtry = c(2:length(names(
    predictors
  )))),
  trControl = trainControl(
    method = "cv",
    index = knndm_folds_random$indx_train,
    savePredictions = TRUE
  )
)

# clustered with spatial leave one cluster out CV
model_clustered <- train(
  trainDat_clustered[, names(predictors)],
  trainDat_clustered$response,
  method = "rf",
  importance = TRUE,
  tuneGrid = expand.grid(mtry = c(2:length(names(
    predictors
  )))),
  trControl = trainControl(
    method = "cv",
    index = knndm_folds_clustered$indx_train,
    savePredictions = TRUE
  )
)

# biased with outlier with 10-fold NNDM CV
model_biasedWithOutlier <-
  train(
    trainDat_biasedWithOutlier[, names(predictors)],
    trainDat_biasedWithOutlier$response,
    method = "rf",
    importance = TRUE,
    tuneGrid = expand.grid(mtry = c(2:length(names(
      predictors
    )))),
    trControl = trainControl(
      method = "cv",
      index = knndm_folds_biasedWithOutlier$indx_train,
      savePredictions = TRUE
    )
  )

# arrange CV results
cv_results = rbind(global_validation(model_random),
                   global_validation(model_clustered),
                   global_validation(model_biasedWithOutlier)) |>
  as.data.frame() |> mutate("CV" = c("10-fold random CV", "10-fold LOCO CV", "10-fold NNDM CV"),
                            "predictors" = c(ncol(model_random$trainingData)-1,
                                             ncol(model_clustered$trainingData)-1,
                                             ncol(model_biasedWithOutlier$trainingData)-1))
# print CV results
knitr::kable(cv_results)

# print variable importances
plot1 = plot(varImp(model_random,scale = F), col="black")
plot2 = plot(varImp(model_clustered,scale = F), col="black")
plot3 = plot(varImp(model_biasedWithOutlier,scale = F), col="black")
pdf(file = "code/figures/simulation_study/varImps.pdf", width = 14, height = 6)
plot_grid(plot1, plot2, plot3, ncol = 3, nrow = 1, labels = c("(a)", "(b)", "(c)"))
invisible(dev.off())


######################################
## Prediction and error calculation ##
######################################
# The trained models are used to make predictions for the entire target area.
# The true absolute error between prediction and reference is calculated for later
# comparison with the local data point density (LPD).

# random
prediction_random <- predict(predictors, model_random, na.rm = T)
truediff_random <- abs(prediction_random-response)
names(truediff_random) <- "absError"

# clustered
prediction_clustered <- predict(predictors, model_clustered, na.rm = T)
truediff_clustered <- abs(prediction_clustered-response)
names(truediff_clustered) <- "absError"

# biased with outlier
prediction_biasedWithOutlier <- predict(predictors, model_biasedWithOutlier, na.rm = T)
truediff_biasedWithOutlier <- abs(prediction_biasedWithOutlier-response)
names(truediff_biasedWithOutlier) <- "absError"


###############################
## Calculate DI, LPD and AOA ##
###############################
# The area of applicability, the dissimilarity index and the local data point
# density are then calculated.

# random
AOA_random <- aoa(newdata = predictors, model = model_random, method = "L2", LPD = TRUE, maxLPD = 1)

# clustered
AOA_clustered <- aoa(newdata = predictors, model = model_clustered, method = "L2", LPD = TRUE, maxLPD = 1)

# biased with outlier
AOA_biasedWithOutlier <- aoa(newdata = predictors, model = model_biasedWithOutlier, method = "L2", LPD = TRUE, maxLPD = 1)

# exploreAOA
# exploreAOA(AOA_random)
# exploreAOA(AOA_clustered)
# exploreAOA(AOA_biasedWithOutlier)


#############################
## LPD ~ SD of RF ensemble ##
#############################
# For camparison to what is often used as uncertainty information,
# the standard deviations of the individual predictions from the
# 500 developed trees within the Random Forest model are calculated.

# RF sd helper function
RFsd <- function(predictors,model){
  prep <- as.data.frame(predictors)
  prep[is.na(prep)] <- -9999
  pred_all <- predict(model$finalModel,prep,predict.all=TRUE)
  sds <-  apply(pred_all$individual,1,sd)
  predsd <- predictors[[1]]
  predsd[!is.na(predsd)] <- sds
  names(predsd) <- "sd"
  return(predsd)
}

# random
predsd_random <- RFsd(predictors,model_random)

# clustered
predsd_clustered <- RFsd(predictors,model_clustered)

# biased with outlier
predsd_biasedWithOutlier <- RFsd(predictors,model_biasedWithOutlier)


#######################
## LPD ~ RMSE tables ##
#######################
# random
LPD_RMSE_random <- data.frame()
for (i in seq(0, AOA_random$parameters$maxLPD)) {
  RMSE = rmse(values(prediction_random)[values(AOA_random$LPD)==i],
              values(response)[values(AOA_random$LPD)==i])

  new_row <- data.frame(LPD = i, RMSE = RMSE)
  LPD_RMSE_random <- rbind(LPD_RMSE_random, new_row)
}

# clustered
LPD_RMSE_clustered <- data.frame()
for (i in seq(0, AOA_clustered$parameters$maxLPD)) {
  RMSE = rmse(values(prediction_clustered)[values(AOA_clustered$LPD)==i],
              values(response)[values(AOA_clustered$LPD)==i])

  new_row <- data.frame(LPD = i, RMSE = RMSE)
  LPD_RMSE_clustered <- rbind(LPD_RMSE_clustered, new_row)
}

# biased with outlier
LPD_RMSE_biasedWithOutlier <- data.frame()
for (i in seq(0, AOA_biasedWithOutlier$parameters$maxLPD)) {
  RMSE = rmse(values(prediction_biasedWithOutlier)[values(AOA_biasedWithOutlier$LPD)==i],
              values(response)[values(AOA_biasedWithOutlier$LPD)==i])

  new_row <- data.frame(LPD = i, RMSE = RMSE)
  LPD_RMSE_biasedWithOutlier <- rbind(LPD_RMSE_biasedWithOutlier, new_row)
}

## RMSE inside and outside AOA
# random
rmse(values(prediction_random)[values(AOA_random$LPD)==0],
     values(response)[values(AOA_random$LPD)==0])
rmse(values(prediction_random)[values(AOA_random$LPD)>0],
     values(response)[values(AOA_random$LPD)>0])

# clustered
rmse(values(prediction_clustered)[values(AOA_clustered$LPD)==0],
     values(response)[values(AOA_clustered$LPD)==0])
rmse(values(prediction_clustered)[values(AOA_clustered$LPD)>0],
     values(response)[values(AOA_clustered$LPD)>0])

# biased with outlier
rmse(values(prediction_biasedWithOutlier)[values(AOA_biasedWithOutlier$LPD)==0],
     values(response)[values(AOA_biasedWithOutlier$LPD)==0])
rmse(values(prediction_biasedWithOutlier)[values(AOA_biasedWithOutlier$LPD)>0],
     values(response)[values(AOA_biasedWithOutlier$LPD)>0])


max_ln <- max(c(length(LPD_RMSE_random$LPD), length(LPD_RMSE_clustered$LPD), length(LPD_RMSE_biasedWithOutlier$LPD)))

# generate table for latex
LPD_RMSE_complete = data.frame(
  LPD_random = c(LPD_RMSE_random$LPD, rep(NA, max_ln - length(
    LPD_RMSE_random$LPD
  ))),
  RMSE_random = c(LPD_RMSE_random$RMSE, rep(NA, max_ln - length(
    LPD_RMSE_random$RMSE
  ))),
  LPD_clustered = c(LPD_RMSE_clustered$LPD, rep(
    NA, max_ln - length(LPD_RMSE_clustered$LPD)
  )),
  RMSE_clustered = c(LPD_RMSE_clustered$RMSE, rep(
    NA, max_ln - length(LPD_RMSE_clustered$RMSE)
  )),
  LPD_biasedWithOutlier = c(LPD_RMSE_biasedWithOutlier$LPD, rep(
    NA, max_ln - length(LPD_RMSE_biasedWithOutlier$LPD)
  )),
  RMSE_biasedWithOutlier = c(LPD_RMSE_biasedWithOutlier$RMSE, rep(
    NA, max_ln - length(LPD_RMSE_biasedWithOutlier$RMSE)
  ))
)
print(xtable(LPD_RMSE_complete), include.rownames=FALSE)



###############################
## LPD ~ Mean Pred SD tables ##
###############################
# random
LPD_predsd_random <- data.frame()
for (i in seq(0, AOA_random$parameters$maxLPD)) {
  meanPredSd = mean(values(predsd_random)[values(AOA_random$LPD)==i], na.rm = T)

  new_row <- data.frame(LPD = i, meanPredSd = meanPredSd)
  LPD_predsd_random <- rbind(LPD_predsd_random, new_row)
}

# clustered
LPD_predsd_clustered <- data.frame()
for (i in seq(0, AOA_clustered$parameters$maxLPD)) {
  meanPredSd = mean(values(predsd_clustered)[values(AOA_clustered$LPD)==i], na.rm = T)

  new_row <- data.frame(LPD = i, meanPredSd = meanPredSd)
  LPD_predsd_clustered <- rbind(LPD_predsd_clustered, new_row)
}

# biased with outlier
LPD_predsd_biasedWithOutlier <- data.frame()
for (i in seq(0, AOA_biasedWithOutlier$parameters$maxLPD)) {
  meanPredSd = mean(values(predsd_biasedWithOutlier)[values(AOA_biasedWithOutlier$LPD)==i], na.rm = T)

  new_row <- data.frame(LPD = i, meanPredSd = meanPredSd)
  LPD_predsd_biasedWithOutlier <- rbind(LPD_predsd_biasedWithOutlier, new_row)
}

max_ln <- max(c(length(LPD_predsd_random$LPD), length(LPD_predsd_clustered$LPD), length(LPD_predsd_biasedWithOutlier$LPD)))

# generate table for latex
LPD_predsd_complete = data.frame(
  LPD_random = c(LPD_predsd_random$LPD, rep(NA, max_ln - length(
    LPD_predsd_random$LPD
  ))),
  meanPredSd_random = c(LPD_predsd_random$meanPredSd, rep(NA, max_ln - length(
    LPD_predsd_random$meanPredSd
  ))),
  LPD_clustered = c(LPD_predsd_clustered$LPD, rep(
    NA, max_ln - length(LPD_predsd_clustered$LPD)
  )),
  meanPredSd_clustered = c(LPD_predsd_clustered$meanPredSd, rep(
    NA, max_ln - length(LPD_predsd_clustered$meanPredSd)
  )),
  LPD_biasedWithOutlier = c(LPD_predsd_biasedWithOutlier$LPD, rep(
    NA, max_ln - length(LPD_predsd_biasedWithOutlier$LPD)
  )),
  meanPredSd_biasedWithOutlier = c(LPD_predsd_biasedWithOutlier$meanPredSd, rep(
    NA, max_ln - length(LPD_predsd_biasedWithOutlier$meanPredSd)
  ))
)
print(xtable(LPD_predsd_complete), include.rownames=FALSE)



##############################################################
## Plot LPD ~ True absolute error & LPD ~ SD of RF ensemble ##
##############################################################
# visualized via data bins

# random
dat_all_random <- data.frame(absError = values(truediff_random, na.rm=T), type = rep("LPD",length(values(AOA_random$LPD, na.rm = T))),
                      LPD = values(AOA_random$LPD, na.rm = T), DI = values(AOA_random$DI, na.rm = T), predsd = values(predsd_random, na.rm = T))
th_random <- cv_results$RMSE[1]


plot1_random = ggplot(dat_all_random, aes(x = LPD, y = absError)) +
  stat_bin_2d(breaks=list(x = seq(-0.25,max(dat_all_random$LPD)+1,0.5), y = seq(0, ceiling(max(dat_all_random$absError)), 0.01))) +
  ylab("RMSE") +
  xlab("Local data point density (LPD)") +
  scale_fill_gradientn(name = "Data points",
                       trans = "log",
                       breaks = 10^(0:3),
                       colors=viridis(10, alpha=.7)) +
  geom_point(data = LPD_RMSE_random, aes(LPD, RMSE, color = "RMSE"), size = .5) +
  scale_color_manual(name = "",values = c("RMSE" = "black")) +
  geom_hline(aes(yintercept=th_random,linetype="CV RMSE"), color = "red") +
  scale_linetype_manual(name = "", values = c("CV RMSE" = "solid")) +
  labs(title = "LPD ~ RMSE (random)") +
  theme_bw() +
  theme(legend.title = element_text( size = 10))


plot2_random = ggplot(dat_all_random, aes(LPD,sd)) +
  stat_bin_2d(breaks=list(x = seq(-0.25,max(dat_all_random$LPD)+1,0.5), y = seq(0, ceiling(max(dat_all_random$sd)), 0.01)))+
  ylab("Prediction SD")+
  xlab("Local data point density (LPD)")+
  scale_fill_gradientn(name = "Data points",
                       trans = "log",
                       breaks = 10^(0:3),
                       colors=viridis(10, alpha=.7))+
  geom_point(data = LPD_predsd_random, aes(LPD, meanPredSd, color = "Mean Pred. SD"), size = .5) +
  scale_color_manual(name = "",values = c("Mean Pred. SD" = "black")) +
  labs(title = "LPD ~ Prediction SD (random)") +
  theme_bw()+
  theme(legend.title = element_text( size = 10))


# clustered
dat_all_clustered <- data.frame(absError = values(truediff_clustered, na.rm=T), type = rep("LPD",length(values(AOA_clustered$LPD, na.rm = T))),
                      LPD = values(AOA_clustered$LPD, na.rm = T), DI = values(AOA_clustered$DI, na.rm = T), predsd = values(predsd_clustered, na.rm = T))
th_clustered <- cv_results$RMSE[2]


plot1_clustered = ggplot(dat_all_clustered, aes(x = LPD, y = absError)) +
  stat_bin_2d(breaks=list(x = seq(-0.25,max(dat_all_clustered$LPD)+1,0.5), y = seq(0, ceiling(max(dat_all_clustered$absError)), 0.01))) +
  ylab("RMSE") +
  xlab("Local data point density (LPD)") +
  scale_fill_gradientn(name = "Data points",
                       trans = "log",
                       breaks = 10^(0:3),
                       colors=viridis(10, alpha=.7)) +
  geom_point(data = LPD_RMSE_clustered, aes(LPD, RMSE, color = "RMSE"), size = .5) +
  scale_color_manual(name = "",values = c("RMSE" = "black")) +
  geom_hline(aes(yintercept=th_clustered,linetype="CV RMSE"), color = "red") +
  scale_linetype_manual(name = "", values = c("CV RMSE" = "solid")) +
  labs(title = "LPD ~ RMSE (clustered)") +
  theme_bw() +
  theme(legend.title = element_text( size = 10))


plot2_clustered = ggplot(dat_all_clustered, aes(LPD,sd)) +
  stat_bin_2d(breaks=list(x = seq(-0.25,max(dat_all_clustered$LPD)+1,0.5), y = seq(0, ceiling(max(dat_all_clustered$sd)), 0.01)))+
  ylab("Prediction SD")+
  xlab("Local data point density (LPD)")+
  scale_fill_gradientn(name = "Data points",
                       trans = "log",
                       breaks = 10^(0:3),
                       colors=viridis(10, alpha=.7))+
  geom_point(data = LPD_predsd_clustered, aes(LPD, meanPredSd, color="Mean Pred. SD"), size = .5) +
  scale_color_manual(name = "",values = c("Mean Pred. SD" = "black")) +
  labs(title = "LPD ~ Prediction SD (clustered)") +
  theme_bw()+
  theme(legend.title = element_text( size = 10))


# biased with outlier
dat_all_biasedWithOutlier <- data.frame(absError = values(truediff_biasedWithOutlier, na.rm=T), type = rep("LPD",length(values(AOA_biasedWithOutlier$LPD, na.rm = T))),
                      LPD = values(AOA_biasedWithOutlier$LPD, na.rm = T), DI = values(AOA_biasedWithOutlier$DI, na.rm = T), predsd = values(predsd_biasedWithOutlier, na.rm = T))
th_biasedWithOutlier <- cv_results$RMSE[3]


plot1_biasedWithOutlier = ggplot(dat_all_biasedWithOutlier, aes(x = LPD, y = absError)) +
  stat_bin_2d(breaks=list(x = seq(-0.25,max(dat_all_biasedWithOutlier$LPD)+1,0.5), y = seq(0, ceiling(max(dat_all_biasedWithOutlier$absError)), 0.01))) +
  ylab("RMSE") +
  xlab("Local data point density (LPD)") +
  scale_fill_gradientn(name = "Data points",
                       trans = "log",
                       breaks = 10^(0:3),
                       colors=viridis(10, alpha=.7)) +
  geom_point(data = LPD_RMSE_biasedWithOutlier, aes(LPD, RMSE, color = "RMSE"), size = .5) +
  scale_color_manual(name = "",values = c("RMSE" = "black")) +
  geom_hline(aes(yintercept=th_biasedWithOutlier,linetype="CV RMSE"), color = "red") +
  scale_linetype_manual(name = "", values = c("CV RMSE" = "solid")) +
  labs(title = "LPD ~ RMSE (biased with outlier)") +
  theme_bw() +
  theme(legend.title = element_text( size = 10))


plot2_biasedWithOutlier =  ggplot(dat_all_biasedWithOutlier, aes(LPD,sd)) +
  stat_bin_2d(breaks=list(x = seq(-0.25,max(dat_all_biasedWithOutlier$LPD)+1,0.5), y = seq(0, ceiling(max(dat_all_biasedWithOutlier$sd)), 0.01)))+
  ylab("Prediction SD")+
  xlab("Local data point density (LPD)")+
  scale_fill_gradientn(name = "Data points",
                       trans = "log",
                       breaks = 10^(0:3),
                       colors=viridis(10, alpha=.7))+
  geom_point(data = LPD_predsd_biasedWithOutlier, aes(LPD, meanPredSd, color = "Mean Pred. SD"), size = .5) +
  scale_color_manual(name = "",values = c("Mean Pred. SD" = "black")) +
  labs(title = "LPD ~ Prediction SD (biased with outlier)") +
  theme_bw()+
  theme(legend.title = element_text( size = 10))


# generate pdf plot for LPD ~ True abs. error
pdf(file = "code/figures/simulation_study/LPD_true_error.pdf", width = 10, height = 4)
plot_grid(plot1_random,
          plot1_clustered,
          plot1_biasedWithOutlier,
          ncol=3,
          labels = c("(a)", "(b)", "(c)"))
invisible(dev.off())

# generate pdf plot for LPD ~
pdf(file = "code/figures/simulation_study/LPD_predsd.pdf", width = 14, height = 6)
plot_grid(plot2_random,
          plot2_clustered,
          plot2_biasedWithOutlier,
          ncol=3,
          labels = c("(a)", "(b)", "(c)"))
invisible(dev.off())



################################################
## LPD as a quantitative uncertainty measure? ##
################################################
# only for the random sampling design

# compute error models
DI_errormodel_random <- errorProfiles(model_random, AOA_random, variable = "DI", calib = "scam", k = 4, window.size = 5)
LPD_errormodel_random <- errorProfiles(model_random, AOA_random, variable = "LPD", calib = "scam", k = 6, window.size = 5)


# DI error model:
recl <- attr(DI_errormodel_random, "performance")

reference_perf <- as.data.frame(list(AOA_random$DI,response,prediction_random))
reference_perf <- reference_perf[order(reference_perf$DI),]
names(reference_perf) <- c("DI","obs","pred")

# use same moving window over the prediction data to get true RMSE
slidingw <- attr(DI_errormodel_random, "performance")
reference_metric <- apply(slidingw,1,function(x){
  x_df <- data.frame(t(x))
  subs_ref <- reference_perf[reference_perf$DI>x_df$ll&
                               reference_perf$DI<x_df$ul,]
  rmse(subs_ref[,"pred"],subs_ref[,"obs"])
})
slidingw$true <- reference_metric
slidingw$what <- "truth"

plot_DI_errormodel_random = plot(DI_errormodel_random) +
  geom_point(data = slidingw, mapping = aes_string(x = "DI", y = "true", color = "what")) +
  scale_linetype_manual(name = "", values = c("model" = "dashed")) +
  scale_color_manual(name = "", values = c("truth" = "red", "cross-validation" = "black")) +
  labs(title = "DI ~ metric (RMSE,random)") +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = "bottom")

slidingw$model = predict(DI_errormodel_random, slidingw)
cor(slidingw$model, slidingw$true)
rmse(slidingw$model, slidingw$true)


# LPD error model:

recl <- attr(LPD_errormodel_random, "performance")

reference_perf <- as.data.frame(list(AOA_random$LPD,response,prediction_random))
reference_perf <- reference_perf[order(reference_perf$LPD),]
names(reference_perf) <- c("LPD","obs","pred")

# use same moving window over the prediction data to get true RMSE
slidingw <- attr(LPD_errormodel_random, "performance")
reference_metric <- apply(slidingw,1,function(x){
  x_df <- data.frame(t(x))
  subs_ref <- reference_perf[reference_perf$LPD>=x_df$ll&
                               reference_perf$LPD<=x_df$ul,]
  rmse(subs_ref[,"pred"],subs_ref[,"obs"])
})
slidingw$true <- reference_metric
slidingw$what <- "truth"

plot_LPD_errormodel_random = plot(LPD_errormodel_random) +
  geom_point(data = slidingw, mapping = aes_string(x = "LPD", y = "true", color = "what")) +
  scale_linetype_manual(name = "", values = c("model" = "dashed")) +
  scale_color_manual(name = "", values = c("truth" = "red", "cross-validation" = "black")) +
  labs(title = "LPD ~ metric (RMSE,random)") +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = "bottom")

slidingw$model = predict(LPD_errormodel_random, slidingw)
cor(slidingw$model, slidingw$true)
rmse(slidingw$model, slidingw$true)

# generate pdf plot for DI and LPD error model
pdf(file = "code/figures/simulation_study/DI_LPD_errormodel_random.pdf", width = 10, height = 5)
plot_grid(plot_DI_errormodel_random,
          plot_LPD_errormodel_random,
          ncol=2,
          labels = c("(a)", "(b)"))
invisible(dev.off())



##################################################
## Estimate model performance with error models ##
##################################################
# prediction only makes sense inside aoa, since model was only fitted for data inside AOA
# inside AOA: DI smaller than AOA threshold and LPD greater than 1

DI_error_prediction_random <- predict(AOA_random$DI, DI_errormodel_random)
LPD_error_prediction_random <- predict(AOA_random$LPD, LPD_errormodel_random)

DI_error_prediction_random[AOA_random$AOA == 0] <- NA
LPD_error_prediction_random[AOA_random$AOA == 0] <- NA
truediff_random[AOA_random$AOA == 0] <- NA

colors <- as.character(values(AOA_random$AOA))
colors[colors==0] <- "violetred"
colors[colors==1] <- "transparent"

plot1 = ggplot() +
  geom_spatraster(data = DI_error_prediction_random) +
  scale_fill_viridis_c(na.value = "transparent") +
  geom_spatraster(data = AOA_random$AOA , fill = colors, na.rm = TRUE, show.legend = T) +
  labs(title = "DI predicted RMSE") +
  theme +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.1, "npc")
  )

plot2 = ggplot() +
  geom_spatraster(data = LPD_error_prediction_random) +
  scale_fill_viridis_c(na.value = "transparent") +
  geom_spatraster(data = AOA_random$AOA , fill = colors, na.rm = TRUE, show.legend = T) +
  labs(title = "LPD predicted RMSE") +
  theme +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.1, "npc")
  )

plot3 = ggplot() +
  geom_spatraster(data = truediff_random) +
  scale_fill_viridis_c(na.value = "transparent") +
  geom_spatraster(data = AOA_random$AOA , fill = colors, na.rm = TRUE, show.legend = T) +
  labs(title = "True RMSE") +
  theme +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.1, "npc")
  )

pdf(file = "code/figures/simulation_study/DI_LPD_error_predictions_random.pdf", width = 14, height = 6)
plot_grid(plot1,
          plot2,
          plot3,
          ncol=3,
          nrow=1,
          labels = c("(a)", "(b)", "(c)"))
invisible(dev.off())



##############################
## correlation coefficients ##
##############################

paste0("correlation coefficient (true error~LPD)= ",
       cor(values(truediff_random),values(AOA_random$LPD),use="complete.obs", method = 'pearson'))
paste0("correlation coefficient (true error~LPD)= ",
       cor(values(truediff_random),values(AOA_random$LPD),use="complete.obs", method = 'spearman'))
paste0("correlation coefficient (pred sd~LPD)= ",
       cor(values(predsd_random),values(AOA_random$LPD),use="complete.obs", method = 'pearson'))
paste0("correlation coefficient (pred sd~LPD)= ",
       cor(values(predsd_random),values(AOA_random$LPD),use="complete.obs", method = 'spearman'))
paste0("correlation coefficient (DI~LPD)= ",
       cor(values(AOA_random$DI),values(AOA_random$LPD),use="complete.obs", method = 'pearson'))
paste0("correlation coefficient (DI~LPD)= ",
       cor(values(AOA_random$DI),values(AOA_random$LPD),use="complete.obs", method = 'spearman'))


paste0("correlation coefficient (true error~LPD)= ",
       cor(values(truediff_clustered),values(AOA_clustered$LPD),use="complete.obs", method = 'pearson'))
paste0("correlation coefficient (true error~LPD)= ",
       cor(values(truediff_clustered),values(AOA_clustered$LPD),use="complete.obs", method = 'spearman'))
paste0("correlation coefficient (pred sd~LPD)= ",
       cor(values(predsd_clustered),values(AOA_clustered$LPD),use="complete.obs", method = 'pearson'))
paste0("correlation coefficient (pred sd~LPD)= ",
       cor(values(predsd_clustered),values(AOA_clustered$LPD),use="complete.obs", method = 'spearman'))
paste0("correlation coefficient (DI~LPD)= ",
       cor(values(AOA_clustered$DI),values(AOA_clustered$LPD),use="complete.obs", method = 'pearson'))
paste0("correlation coefficient (DI~LPD)= ",
       cor(values(AOA_clustered$DI),values(AOA_clustered$LPD),use="complete.obs", method = 'spearman'))


paste0("correlation coefficient (true error~LPD)= ",
       cor(values(truediff_biasedWithOutlier),values(AOA_biasedWithOutlier$LPD),use="complete.obs", method = 'pearson'))
paste0("correlation coefficient (true error~LPD)= ",
       cor(values(truediff_biasedWithOutlier),values(AOA_biasedWithOutlier$LPD),use="complete.obs", method = 'spearman'))
paste0("correlation coefficient (pred sd~LPD)= ",
       cor(values(predsd_biasedWithOutlier),values(AOA_biasedWithOutlier$LPD),use="complete.obs", method = 'pearson'))
paste0("correlation coefficient (pred sd~LPD)= ",
       cor(values(predsd_biasedWithOutlier),values(AOA_biasedWithOutlier$LPD),use="complete.obs", method = 'spearman'))
paste0("correlation coefficient (DI~LPD)= ",
       cor(values(AOA_biasedWithOutlier$DI),values(AOA_biasedWithOutlier$LPD),use="complete.obs", method = 'pearson'))
paste0("correlation coefficient (DI~LPD)= ",
       cor(values(AOA_biasedWithOutlier$DI),values(AOA_biasedWithOutlier$LPD),use="complete.obs", method = 'spearman'))




######################################
## LPD to detect outlier-caused AOA ##
######################################

plot_LPD <- ggplot() +
  geom_spatraster(data = AOA_biasedWithOutlier$LPD) +
  scale_fill_viridis_c(na.value = "transparent") +
  geom_sf(data = samplepoints_biasedWithOutlier, color = "red", size = 2) +
  labs(title = "LPD (biased with outlier)") +
  theme +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.1, "npc")
  )


dummyrast <- AOA_biasedWithOutlier$AOA
colors <- as.character(values(AOA_biasedWithOutlier$AOA))
colors[colors==0] <- "violetred"
colors[is.na(colors)] <- "transparent"
colors[colors==1] <- "transparent"

plot_AOA <- ggplot() +
  geom_spatraster(data = prediction_biasedWithOutlier) +
  geom_spatraster(data = dummyrast , fill = colors, na.rm = TRUE) +
  scale_fill_viridis_c(na.value = "transparent") +
  labs(title = "Prediction for AOA (min. LPD = 1)") +
  theme +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.1, "npc")
  )


dummyrast <- AOA_biasedWithOutlier$LPD
colors <- as.character(values(AOA_biasedWithOutlier$LPD))
colors[colors>1] <- "transparent"
colors[is.na(colors)] <- "transparent"
colors[colors<2] <- "violetred"

plot_AOA_LPD <- ggplot() +
  geom_spatraster(data = prediction_biasedWithOutlier) +
  geom_spatraster(data = dummyrast , fill = colors, na.rm = TRUE) +
  scale_fill_viridis_c(na.value = "transparent") +
  labs(title = "Prediction for AOA (min. LPD = 2)") +
  theme +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.1, "npc")
  )

pdf(file = "code/figures/simulation_study/filter_outliers.pdf", width = 14, height = 6)
plot_grid(plot_LPD,
          plot_AOA,
          plot_AOA_LPD,
          ncol=3,
          labels = c("(a)", "(b)", "(c)"))
invisible(dev.off())












#############################################################
## Plot all layers for all sampling scenarios for appendix ##
#############################################################
plot0 = ggplot() +
  geom_spatraster(data = response) +
  scale_fill_viridis_c(na.value = "transparent") +
  labs(title = "Response") +
  theme

# random
dummyrast <- AOA_random$AOA
colors <- as.character(values(AOA_random$AOA))
colors[colors==0] <- "violetred"
colors[is.na(colors)] <- "transparent"
colors[colors==1] <- "transparent"

plot1_random = ggplot() +
  geom_spatraster(data = prediction_random) +
  scale_fill_viridis_c(na.value = "transparent") +
  labs(title = "Prediction (random)") +
  theme +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.1, "npc")
  )

plot2_random <- ggplot() +
  geom_spatraster(data = predsd_random) +
  scale_fill_viridis_c(na.value = "transparent") +
  labs(title = "Prediction SD (random)") +
  theme +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.1, "npc")
  )

plot3_random <- ggplot() +
  geom_spatraster(data = truediff_random) +
  scale_fill_viridis_c(na.value = "transparent") +
  labs(title = "True abs. error (random)") +
  theme +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.1, "npc")
  )

plot4_random <- ggplot() +
  geom_spatraster(data = AOA_random$DI) +
  scale_fill_viridis_c(na.value = "transparent") +
  labs(title = "DI (random)") +
  theme +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.1, "npc")
  )

plot5_random <- ggplot() +
  geom_spatraster(data = AOA_random$LPD) +
  scale_fill_viridis_c(na.value = "transparent") +
  labs(title = "LPD (random)") +
  theme +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.1, "npc")
  )

plot6_random <- ggplot() +
  geom_spatraster(data = prediction_random) +
  geom_spatraster(data = dummyrast , fill = colors, na.rm = TRUE) +
  scale_fill_viridis_c(na.value = "transparent") +
  labs(title = "Prediction for AOA (random)") +
  theme +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.1, "npc")
  )

# clustered
dummyrast <- AOA_clustered$AOA
colors <- as.character(values(AOA_clustered$AOA))
colors[colors==0] <- "violetred"
colors[is.na(colors)] <- "transparent"
colors[colors==1] <- "transparent"

plot1_clustered = ggplot() +
  geom_spatraster(data = prediction_clustered) +
  scale_fill_viridis_c(na.value = "transparent") +
  labs(title = "Prediction (clustered)") +
  theme +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.1, "npc")
  )


plot2_clustered <- ggplot() +
  geom_spatraster(data = predsd_clustered) +
  scale_fill_viridis_c(na.value = "transparent") +
  labs(title = "Prediction SD (clustered)") +
  theme +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.1, "npc")
  )

plot3_clustered <- ggplot() +
  geom_spatraster(data = truediff_clustered) +
  scale_fill_viridis_c(na.value = "transparent") +
  labs(title = "True abs. error (clustered)") +
  theme +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.1, "npc")
  )

plot4_clustered <- ggplot() +
  geom_spatraster(data = AOA_clustered$DI) +
  scale_fill_viridis_c(na.value = "transparent") +
  labs(title = "DI (clustered)") +
  theme +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.1, "npc")
  )

plot5_clustered <- ggplot() +
  geom_spatraster(data = AOA_clustered$LPD) +
  scale_fill_viridis_c(na.value = "transparent") +
  labs(title = "LPD (clustered)") +
  theme +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.1, "npc")
  )

plot6_clustered <- ggplot() +
  geom_spatraster(data = prediction_clustered) +
  geom_spatraster(data = dummyrast , fill = colors, na.rm = TRUE) +
  scale_fill_viridis_c(na.value = "transparent") +
  labs(title = "Prediction for AOA (clustered)") +
  theme +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.1, "npc")
  )

# biased with outlier
dummyrast <- AOA_biasedWithOutlier$AOA
colors <- as.character(values(AOA_biasedWithOutlier$AOA))
colors[colors==0] <- "violetred"
colors[is.na(colors)] <- "transparent"
colors[colors==1] <- "transparent"

plot1_biasedWithOutlier = ggplot() +
  geom_spatraster(data = prediction_biasedWithOutlier) +
  scale_fill_viridis_c(na.value = "transparent") +
  labs(title = "Prediction (biased with outlier)") +
  theme +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.1, "npc")
  )

plot2_biasedWithOutlier <- ggplot() +
  geom_spatraster(data = predsd_biasedWithOutlier) +
  scale_fill_viridis_c(na.value = "transparent") +
  labs(title = "Prediction SD (biased with outlier)") +
  theme +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.1, "npc")
  )

plot3_biasedWithOutlier <- ggplot() +
  geom_spatraster(data = truediff_biasedWithOutlier) +
  scale_fill_viridis_c(na.value = "transparent") +
  labs(title = "True abs. error (biased with outlier)") +
  theme +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.1, "npc")
  )

plot4_biasedWithOutlier <- ggplot() +
  geom_spatraster(data = AOA_biasedWithOutlier$DI) +
  scale_fill_viridis_c(na.value = "transparent") +
  labs(title = "DI (biased with outlier)") +
  theme +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.1, "npc")
  )

plot5_biasedWithOutlier <- ggplot() +
  geom_spatraster(data = AOA_biasedWithOutlier$LPD) +
  scale_fill_viridis_c(na.value = "transparent") +
  labs(title = "LPD (biased with outlier)") +
  theme +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.1, "npc")
  )

plot6_biasedWithOutlier <- ggplot() +
  geom_spatraster(data = prediction_biasedWithOutlier) +
  geom_spatraster(data = dummyrast , fill = colors, na.rm = TRUE) +
  scale_fill_viridis_c(na.value = "transparent") +
  labs(title = "Prediction for AOA (biased with outlier)") +
  theme +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.1, "npc")
  )

pdf(file = "code/figures/simulation_study/response.pdf", width = 6, height = 6)
plot0
invisible(dev.off())

pdf(file = "code/figures/simulation_study/predictions.pdf", width = 14, height = 6)
plot_grid(plot1_random,
          plot1_clustered,
          plot1_biasedWithOutlier,
          ncol=3,
          labels = c("(a)", "(b)", "(c)"))
invisible(dev.off())

pdf(file = "code/figures/simulation_study/predsd.pdf", width = 14, height = 6)
plot_grid(plot2_random,
          plot2_clustered,
          plot2_biasedWithOutlier,
          ncol=3,
          labels = c("(a)", "(b)", "(c)"))
invisible(dev.off())

pdf(file = "code/figures/simulation_study/true_abs_error.pdf", width = 14, height = 6)
plot_grid(plot3_random,
          plot3_clustered,
          plot3_biasedWithOutlier,
          ncol=3,
          labels = c("(a)", "(b)", "(c)"))
invisible(dev.off())

pdf(file = "code/figures/simulation_study/DI.pdf", width = 14, height = 6)
plot_grid(plot4_random,
          plot4_clustered,
          plot4_biasedWithOutlier,
          ncol=3,
          labels = c("(a)", "(b)", "(c)"))
invisible(dev.off())

pdf(file = "code/figures/simulation_study/LPD.pdf", width = 14, height = 6)
plot_grid(plot5_random,
          plot5_clustered,
          plot5_biasedWithOutlier,
          ncol=3,
          labels = c("(a)", "(b)", "(c)"))
invisible(dev.off())

pdf(file = "code/figures/simulation_study/predsAOA.pdf", width = 14, height = 6)
plot_grid(plot6_random,
          plot6_clustered,
          plot6_biasedWithOutlier,
          ncol=3,
          labels = c("(a)", "(b)", "(c)"))
invisible(dev.off())


