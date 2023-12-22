######################
## Simulation study ##
######################
# This simulation study is in parts based on the simulation study of
# Meyer and Pebesma (2021): 
# - https://doi.org/10.1111/2041-210X.13650
# - https://github.com/HannaMeyer/MEE_AOA

rm(list=ls())

# data handling
remotes::install_github("fab-scm/CAST")
library(CAST)
library(virtualspecies)
library(caret)
library(terra)
library(sf)

# visualization
library(viridis)
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

# Presettings for the case study
design <- "random"                                                              # design: either random, clustered, biasedWithOutlier
npoints <- 100                                                                  # npoints: number of training samples
nclusters <- 10                                                                 # nclusters: number of clusters if design==clustered
maxdist <- 0.8                                                                  # maxdist: maxdist for clustered samples if design==clustered
countries <- c("Germany","Ireland","France", "Sweden")                          # countries: main sample countries if design==biasedWithOutlier
countriesOutlier <- "Turkmenistan"                                              # countriesOutlier: outlier country if design==biasedWithOutlier (a single point is set here)
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
modeldomain <- mask

if (design=="random"){
  set.seed(seed)
  samplepoints <- st_as_sf(st_sample(mask, size = npoints, "random"))
  st_write(samplepoints, "./data/samples/samplepoints_sim_study_random.geojson", driver = "GeoJSON", append = FALSE, delete_dsn = TRUE)
}
if (design=="clustered"){
  set.seed(9352)
  mask <- st_transform(mask, "+init=epsg:3857")
  samplepoints <- clustered_sample(mask, npoints, nclusters,radius=200000)
  samplepoints <- st_transform(samplepoints, "+init=epsg:4326")
  mask <- st_transform(mask, "+init=epsg:4326")
  st_write(samplepoints, "./data/samples/samplepoints_sim_study_clustered.geojson", driver = "GeoJSON", append = FALSE, delete_dsn = TRUE)
}
if (design=="biasedWithOutlier"){
  set.seed(7)
  countryboundaries <- geodata::world(res = 4, path = './data/countryboundaries/')
  countryboundariesOut <- st_as_sf(countryboundaries[countryboundaries$NAME_0%in%c(countriesOutlier),])
  countryboundaries <- st_as_sf(countryboundaries[countryboundaries$NAME_0%in%c(countries),])
  samplepoints <- st_as_sf(st_sample(countryboundaries,npoints-1,"random"))
  samplePointsOut <- st_as_sf(st_sample(countryboundariesOut,1,"random"))
  samplepoints <- rbind(samplepoints, samplePointsOut)
  samplepoints <- st_transform(samplepoints, st_crs(mask))
  st_write(samplepoints, "./data/samples/samplepoints_sim_study_biasedWithOutlier.geojson", driver = "GeoJSON", append = FALSE, delete_dsn = TRUE)
}

##############################################
##  Plot all predictors and sampling designs #
##############################################
# add response to predictors for plotting
predictors$response <- response

ggplot() +
  geom_spatraster(data = stretch(predictors,0,1)) +
  facet_wrap(~lyr, ncol = 5) +
  scale_fill_viridis(na.value = "transparent") +
  theme

ggplot() +
  geom_spatvector(data = maskVect, show.legend = F, fill = "grey") +
  labs(title = "Sampling design") +
  geom_sf(data = samplepoints, color = "red", size = 1) +
  theme

# remove response from predictos
predictors <- subset(predictors, predictor_names)

################################
## Prepare training data sets ##
################################
# extract predictors for the sampling locations
trainDat <- extract(predictors, samplepoints, df=TRUE)
trainDat$response <- extract(response, samplepoints, ID = FALSE)$response
if (design == "clustered") {
  trainDat$clstrID <- samplepoints$parent
}
trainDat <- trainDat[complete.cases(trainDat),]


####################
## Model training ##
####################
# Model training is then performed using the caret package.
# The model output gives information on the general estimated
# model performance based on the cross-validation.

set.seed(seed)
if(design =="random"){
  # random with random cv
  model <- train(trainDat[,names(predictors)],
                 trainDat$response,
                 method="rf",
                 importance=TRUE,
                 tuneGrid = expand.grid(mtry = c(2:length(names(predictors)))),
                 trControl = trainControl(method="cv",savePredictions = TRUE))
}
if(design=="clustered"){
  # clustered with spatial leave one cluster out cv
  folds <- CreateSpacetimeFolds(trainDat, spacevar="clstrID", k=nclusters)
  model <- train(trainDat[,names(predictors)],
                 trainDat$response,
                 method="rf",
                 importance=TRUE,
                 tuneGrid = expand.grid(mtry = c(2:length(names(predictors)))),
                 trControl = trainControl(method="cv",index=folds$index,savePredictions = TRUE))
}
if(design == "biasedWithOutlier"){
  # biased with outlier with kNNDM CV
  knndm_folds <- knndm(samplepoints, modeldomain = mask, k = 10)
  model <- train(trainDat[,names(predictors)],
                 trainDat$response,
                 method="rf",
                 importance=TRUE,
                 tuneGrid = expand.grid(mtry = c(2:length(names(predictors)))),
                 trControl = trainControl(method="cv",index = knndm_folds$indx_train,savePredictions = TRUE))
}

if (design == "random") {
  cv <- "random"
} else if (design == "clustered") {
  cv <- "spatial"
} else if (design == "biasedWithOutliers") {
  cv <- "knndm"
}
# print cv results
cv_results = rbind(global_validation(model)) |>
  as.data.frame() |> mutate("CV" = c(cv),
                            "predictors" = c(ncol(model$trainingData)-1))
# print CV results
knitr::kable(cv_results)

# plot var importance
plot(varImp(model,scale = F), col="black")


######################################
## Prediction and error calculation ##
######################################
# The trained models are used to make predictions for the entire target area.
# The true absolute error between prediction and reference is calculated for later
# comparison with the local data point density (LPD).

prediction <- predict(predictors, model, na.rm = T)
truediff <- abs(prediction-response)
names(truediff) <- "absError"

ggplot() +
  geom_spatraster(data = response) +
  scale_fill_viridis(na.value = "transparent") +
  labs(title = "Response") +
  theme

ggplot() +
  geom_spatraster(data = prediction) +
  scale_fill_viridis(na.value = "transparent") +
  labs(title = "Prediction") +
  theme

ggplot() +
  geom_spatraster(data = truediff) +
  scale_fill_viridis(na.value = "transparent", option = "rocket", begin = 0.2) +
  labs(title = "True abs. error") +
  theme


###############################
## Calculate DI, LPD and AOA ##
###############################
# The area of applicability, the dissimilarity index and the local data point
# density are then calculated.

AOA <- aoa(newdata = predictors, model = model, method = "L2", LPD = TRUE, maxLPD = 1) # max LPD = 1 => use 100% of the training data to calculate LPD
# plot(AOA)[[1]]
# plot(AOA)[[2]]


ggplot() +
  geom_spatraster(data = AOA$LPD) +
  scale_fill_viridis_c(na.value = "transparent") +
  labs(title = "LPD") +
  theme



# # exploreAOA
# library(leaflet)
# library(shiny)
# library(shinycssloaders)
# library(rlist)
# library(bslib)
# library(shinyWidgets)
# exploreAOA(AOA)


#################################################################
## calculate prediction standard deviations of the RF ensemble ##
#################################################################
# For camparison to what is often used as uncertainty information,
# the standard deviations of the individual predictions from the
# 500 developed trees within the Random Forest model are calculated.

# pred sd helper function
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

predsd <- RFsd(predictors,model)

ggplot() +
  geom_spatraster(data = predsd) +
  scale_fill_viridis(na.value = "transparent") +
  labs(title = "Prediction sd") +
  theme


## Comparison between LPD and true error/prediction standard deviations/DI
paste0("correlation coefficient (true error~LPD)= ",
       cor(values(truediff),values(AOA$LPD),use="complete.obs", method = 'pearson'))
paste0("correlation coefficient (true error~LPD)= ",
       cor(values(truediff),values(AOA$LPD),use="complete.obs", method = 'spearman'))
paste0("correlation coefficient (DI~LPD)= ",
       cor(values(predsd),values(AOA$LPD),use="complete.obs", method = 'pearson'))
paste0("correlation coefficient (DI~LPD)= ",
       cor(values(predsd),values(AOA$LPD),use="complete.obs", method = 'spearman'))
paste0("correlation coefficient (DI~LPD)= ",
       cor(values(AOA$DI),values(AOA$LPD),use="complete.obs", method = 'pearson'))
paste0("correlation coefficient (DI~LPD)= ",
       cor(values(AOA$DI),values(AOA$LPD),use="complete.obs", method = 'spearman'))

compare <- rast(list(
  response,
  prediction,
  predsd,
  truediff,
  AOA$DI,
  AOA$LPD
))
names(compare) <- c("response", "prediction", "sd", "absError", "DI", "LPD")
summary(values(compare))




###########################
## LPD ~ RMSE data frame ##
###########################
# Create data frame of RMSE of prediction locations summarized by their LPD value
LPD_RMSE <- data.frame()
for (i in seq(0, AOA$parameters$maxLPD)) {
  RMSE = rmse(values(prediction)[values(AOA$LPD)==i],
              values(response)[values(AOA$LPD)==i])

  new_row <- data.frame(LPD = i, RMSE = RMSE)
  LPD_RMSE <- rbind(LPD_RMSE, new_row)
}

## RMSE inside and outside AOA
rmse(values(prediction)[values(AOA$LPD)==0],
     values(response)[values(AOA$LPD)==0])
rmse(values(prediction)[values(AOA$LPD)>0],
     values(response)[values(AOA$LPD)>0])


####################################
## LPD ~ Mean pred. SD data frame ##
####################################
# Create data frame of mean prediction standard deviations of prediction locations
# summarized by their LPD values
LPD_predsd <- data.frame()
for (i in seq(0, AOA$parameters$maxLPD)) {
  meanPredSd = mean(values(predsd)[values(AOA$LPD)==i], na.rm = T)

  new_row <- data.frame(LPD = i, meanPredSd = meanPredSd)
  LPD_predsd <- rbind(LPD_predsd, new_row)
}


# Create data frame from: True abs. error | LPD | Pred. SD | DI
dat_all <- data.frame(absError = values(truediff, na.rm=T), LPD = values(AOA$LPD, na.rm = T), predsd = values(predsd, na.rm = T), DI = round(values(AOA$DI, na.rm = TRUE),digits = 3))
th <- cv_results$RMSE[1]


####################################
## Plot LPD ~ True absolute error ##
####################################
# Relationship between the LPD and the absolute error on a data-point-level
ggplot(dat_all, aes(x = LPD, y = absError)) +
  stat_bin_2d(breaks=list(x = seq(-0.25,max(dat_all$LPD)+1,0.5), y = seq(0, ceiling(max(dat_all$absError)), 0.01))) +
  ylab("True absolute error / RMSE") +
  xlab("Local data point density (LPD)") +
  scale_fill_gradientn(name = "Data points",
                       trans = "log",
                       breaks = 10^(0:3),
                       colors=viridis(10, alpha=.7)) +
  geom_point(data = LPD_RMSE, aes(LPD, RMSE, color = "RMSE"), size = .5) +
  scale_color_manual(name = "",values = c("RMSE" = "black")) +
  geom_hline(aes(yintercept=th,linetype="CV RMSE"), color = "red") +
  scale_linetype_manual(name = "", values = c("CV RMSE" = "solid")) +
  labs(title = "LPD ~ true abs. error") +
  theme_bw() +
  theme(legend.title = element_text( size = 10))
# stat_smooth(data = dat_all, aes(LPD, absError), method = "nls", formula = y~a*exp(b*x), method.args = list(start=c(a=1, b=0)), se = FALSE, col = "red")


#############################
## LPD ~ SD of RF ensemble ##
#############################
# Relationship between the LPD and the SD of RF ensemble on a data-point-level
ggplot(dat_all, aes(LPD,sd)) +
  stat_bin_2d(breaks=list(x = seq(-0.25,max(dat_all$LPD)+1,0.5), y = seq(0, ceiling(max(dat_all$sd)), 0.01)))+
  ylab("Prediction SD")+
  xlab("Local data point density (LPD)")+
  scale_fill_gradientn(name = "Data points",
                       trans = "log",
                       breaks = 10^(0:3),
                       colors=viridis(10, alpha=.7))+
  geom_point(data = LPD_predsd, aes(LPD, meanPredSd, color = "Mean Pred. SD"), size = .5) +
  scale_color_manual(name = "",values = c("Mean Pred. SD" = "black")) +
  labs(title = "LPD ~ Prediction SD (random)") +
  theme_bw()+
  theme(legend.title = element_text( size = 10))
# stat_smooth(data = dat_all, aes(LPD, sd), method = "nls", formula = y~a*exp(b*x), method.args = list(start=c(a=1, b=0)), se = FALSE, col = "red")


###################
## Plot LPD ~ DI ##
###################
# Relationship between the LPD and DI on a data-point-level
dat_all[dat_all$LPD == 0 & dat_all$DI <= AOA$parameters$threshold, "DI"] = AOA$parameters$threshold + 0.01 # shift values for clear bins
ggplot(dat_all, aes(x = LPD, y = DI)) +
  stat_bin_2d(breaks = list(x = seq(-0.5, max(dat_all$LPD)+1, 1), y = seq(0, ceiling(max(dat_all$DI)), 0.05))) +
  scale_fill_viridis(begin = 0.1) +
  geom_hline(aes(yintercept = AOA$parameters$threshold, linetype = "AOA_threshold")) +
  scale_linetype_manual(name = "", values = c(AOA_threshold = "dashed")) +
  labs(title = "LPD ~ DI") +
  theme_bw() +
  stat_smooth(data = dat_all, aes(LPD, DI, color = "LPD~DI"), method = "gam", se = FALSE) +
  scale_color_manual(name = "",values = c("LPD~DI" = "red"))



################################################
## LPD as a quantitative uncertainty measure? ##
################################################

DI_errormodel <- DItoErrormetric(model, AOA, calib = "scam", k = 4)
LPD_errormodel <- LPDtoErrormetric(model, AOA, calib = "scam")
DI_LPD_errormodel <- DI_LPDtoErrormetric(model, AOA, calib = "scam")

# DI error model:
recl <- attr(DI_errormodel, "performance")

reference_perf <- as.data.frame(list(AOA$DI,response,prediction))
reference_perf <- reference_perf[order(reference_perf$DI),]
names(reference_perf) <- c("DI","obs","pred")

# use same moving window over the prediction data to get true RMSE
slidingw <- attr(DI_errormodel, "performance")
reference_metric <- apply(slidingw,1,function(x){
  x_df <- data.frame(t(x))
  subs_ref <- reference_perf[reference_perf$DI>x_df$ll&
                               reference_perf$DI<x_df$ul,]
  rmse(subs_ref[,"pred"],subs_ref[,"obs"])
})
slidingw$true <- reference_metric
slidingw$what <- "truth"

plot(DI_errormodel) +
  geom_point(data = slidingw, mapping = aes_string(x = "DI", y = "true", color = "what")) +
  scale_linetype_manual(name = "", values = c("model" = "dashed")) +
  scale_color_manual(name = "", values = c("truth" = "red", "cross-validation" = "black")) +
  labs(title = "DI ~ metric (RMSE)") +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = "bottom")

slidingw$model = predict(DI_errormodel, slidingw)
cor(slidingw$model, slidingw$true)
rmse(slidingw$model, slidingw$true)


# LPD error model:
recl <- attr(LPD_errormodel, "performance")

reference_perf <- as.data.frame(list(AOA$LPD,response,prediction))
reference_perf <- reference_perf[order(reference_perf$LPD),]
names(reference_perf) <- c("LPD","obs","pred")

# use same moving window over the prediction data to get true RMSE
slidingw <- attr(LPD_errormodel, "performance")
reference_metric <- apply(slidingw,1,function(x){
  x_df <- data.frame(t(x))
  subs_ref <- reference_perf[reference_perf$LPD>=x_df$ll&
                               reference_perf$LPD<=x_df$ul,]
  rmse(subs_ref[,"pred"],subs_ref[,"obs"])
})
slidingw$true <- reference_metric
slidingw$what <- "truth"

plot(LPD_errormodel) +
  geom_point(data = slidingw, mapping = aes_string(x = "LPD", y = "true", color = "what")) +
  scale_linetype_manual(name = "", values = c("model" = "dashed")) +
  scale_color_manual(name = "", values = c("truth" = "red", "cross-validation" = "black")) +
  labs(title = "LPD ~ metric (RMSE)") +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = "bottom")

slidingw$model = predict(LPD_errormodel, slidingw)
cor(slidingw$model, slidingw$true)
rmse(slidingw$model, slidingw$true)


# DI + LPD error model:
recl <- attr(DI_LPD_errormodel, "performance")

reference_perf_DI <- as.data.frame(list(AOA$DI,response,prediction))
reference_perf_DI <- reference_perf_DI[order(reference_perf_DI$DI),]
names(reference_perf_DI) <- c("DI","obs","pred")

# use same moving window over the prediction data to get true RMSE
slidingw_DI <- attr(DI_LPD_errormodel, "performance")
reference_metric_DI <- apply(slidingw_DI,1,function(x){
  x_df <- data.frame(t(x))
  subs_ref <- reference_perf_DI[reference_perf_DI$DI>x_df$ll.DI&
                                  reference_perf_DI$DI<x_df$ul.DI,]
  rmse(subs_ref[,"pred"],subs_ref[,"obs"])
})
slidingw_DI$true.DI <- reference_metric_DI

reference_perf_LPD <- as.data.frame(list(AOA$LPD,response,prediction))
reference_perf_LPD <- reference_perf_LPD[order(reference_perf_LPD$LPD),]
names(reference_perf_LPD) <- c("LPD","obs","pred")

# use same moving window over the prediction data to get true RMSE
slidingw_LPD <- attr(DI_LPD_errormodel, "performance")
reference_metric_LPD <- apply(slidingw_LPD,1,function(x){
  x_df <- data.frame(t(x))
  subs_ref <- reference_perf_LPD[reference_perf_LPD$LPD>=x_df$ll.LPD&
                                   reference_perf_LPD$LPD<=x_df$ul.LPD,]
  rmse(subs_ref[,"pred"],subs_ref[,"obs"])
})
slidingw_LPD$true.LPD <- reference_metric_LPD

# calculate true RMSE from DI and LPD moving window
slidingw <- merge(slidingw_DI, slidingw_LPD)
slidingw$true <- (slidingw$true.DI + slidingw$true.LPD) / 2

plot(DI_LPD_errormodel) %>%
  layout(title = list(text = "DI + LPD ~ metric (RMSE)", x = 0.2, y = 0.8)) %>%
  add_markers(x = ~DI,
              y = ~LPD,
              z = ~true,
              data = slidingw,
              color = I("red"),
              size = I(40),
              name = "truth",
              opacity = 1,
              hovertemplate = paste0("DI: %{x}<br>LPD: %{y}<br>metric: %{z}<br>"))

slidingw$model = predict(DI_LPD_errormodel_random, slidingw)
cor(slidingw$model, slidingw$true)
rmse(slidingw$model, slidingw$true)



##################################################
## Estimate model performance with error models ##
##################################################
# prediction only makes sense inside aoa, as model was only fitted for data inside AOA
# DI: smaller than AOA threshold
# LPD: bigger than 1

DI_error_prediction <- predict(AOA$DI, DI_errormodel)
LPD_error_prediction <- predict(AOA$LPD, LPD_errormodel)
DI_LPD_error_prediction <- predict(rast(list(AOA$DI, AOA$LPD)), DI_LPD_errormodel)

# mask outside AOA values
DI_error_prediction[AOA$AOA == 0] <- NA
LPD_error_prediction[AOA$AOA == 0] <- NA
DI_LPD_error_prediction[AOA$AOA == 0] <- NA
truediff[AOA$AOA == 0] <- NA

colors <- as.character(values(AOA$AOA))
colors[colors==0] <- "violetred"
colors[colors==1] <- "transparent"

ggplot() +
  geom_spatraster(data = DI_error_prediction) +
  scale_fill_viridis_c(na.value = "transparent") +
  geom_spatraster(data = AOA$AOA , fill = colors, na.rm = TRUE, show.legend = T) +
  labs(title = "DI predicted RMSE") +
  theme

ggplot() +
  geom_spatraster(data = LPD_error_prediction) +
  scale_fill_viridis_c(na.value = "transparent") +
  geom_spatraster(data = AOA$AOA , fill = colors, na.rm = TRUE, show.legend = T) +
  labs(title = "LPD predicted RMSE") +
  theme

ggplot() +
  geom_spatraster(data = DI_LPD_error_prediction) +
  scale_fill_viridis_c(na.value = "transparent") +
  geom_spatraster(data = AOA$AOA , fill = colors, na.rm = TRUE, show.legend = T) +
  labs(title = "DI and LPD predicted RMSE") +
  theme

# compare performance predictions with true perfromance (true abs. error)
ggplot() +
  geom_spatraster(data = truediff) +
  scale_fill_viridis_c(na.value = "transparent") +
  geom_spatraster(data = AOA$AOA , fill = colors, na.rm = TRUE, show.legend = T) +
  labs(title = "True RMSE (True abs. error)") +
  theme



######################################
## LPD to detect outlier-caused AOA ##
######################################
if (design == "biasedWithOutlier") {
  plot_LPD <- ggplot() +
    geom_spatraster(data = AOA$LPD) +
    scale_fill_viridis_c(na.value = "transparent") +
    geom_sf(data = samplepoints, color = "red", size = 1) +
    labs(title = "LPD with samples") +
    theme


  dummyrast <- AOA$AOA
  colors <- as.character(values(AOA$AOA))
  colors[colors==0] <- "violetred"
  colors[is.na(colors)] <- "transparent"
  colors[colors==1] <- "transparent"

  ggplot() +
    geom_spatraster(data = prediction) +
    geom_spatraster(data = dummyrast , fill = colors, na.rm = TRUE) +
    scale_fill_viridis_c(na.value = "transparent") +
    labs(title = "Prediction for AOA (min. LPD = 1)") +
    theme


  dummyrast <- AOA$LPD
  colors <- as.character(values(AOA$LPD))
  colors[colors>1] <- "transparent"
  colors[is.na(colors)] <- "transparent"
  colors[colors<2] <- "violetred"

  ggplot() +
    geom_spatraster(data = prediction) +
    geom_spatraster(data = dummyrast , fill = colors, na.rm = TRUE) +
    scale_fill_viridis_c(na.value = "transparent") +
    labs(title = "Prediction for AOA (min. LPD = 2)") +
    theme
}

