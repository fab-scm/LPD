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
library(reshape2)
library(ggpubr)
library(cowplot)
library(knitr)
library(xtable)

# r_squared, rmse and mae helper functions
rmse <- function(pred,obs){sqrt( mean((pred - obs)^2, na.rm = T) )}
r_squared <- function(pred,obs){
  mean_obs <- mean(obs, na.rm = T)                                              # Calculate mean of the true values
  ss_total <- sum((obs - mean_obs)^2, na.rm = T)                                # Calculate the total sum of squares
  ss_residual <- sum((obs - pred)^2, na.rm = T)                                 # Calculate the residual sum of squares
  res <- 1 - (ss_residual / ss_total)                                           # Calculate R^2
  return(res)
}
mae <- function(pred,obs){mean(abs(obs - pred), na.rm = T)}

get_errors <- function(pred,obs){
  data.frame(RMSE = rmse(pred, obs),
             Rsquared = r_squared(pred, obs),
             MAE = mae(pred, obs))
}

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
maxdist <- 200000                                                               # maxdist: maxdist for clustered samples if design==clustered (in metres)
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
samplepoints_random <- st_as_sf(st_sample(mask, size = npoints, type = "random"))
samplepoints_random <- st_transform(samplepoints_random, "+init=epsg:4326")
st_write(samplepoints_random, "./data/samples/samplepoints_sim_study_random.geojson", driver = "GeoJSON", append = FALSE, delete_dsn = TRUE)

# clustered
set.seed(9352)
mask <- st_transform(mask, "+init=epsg:3857")
samplepoints_clustered <- clustered_sample(mask, npoints, nclusters,radius=maxdist)
samplepoints_clustered <- st_transform(samplepoints_clustered, "+init=epsg:4326")
mask <- st_transform(mask, "+init=epsg:4326")
st_write(samplepoints_clustered, "./data/samples/samplepoints_sim_study_clustered.geojson", driver = "GeoJSON", append = FALSE, delete_dsn = TRUE)


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
  labs(title = "random") +
  geom_sf(data = samplepoints_random, color = "red", size = 0.5) +
  theme

plot2 = ggplot() +
  geom_spatvector(data = maskVect, show.legend = F, fill = "grey") +
  labs(title = "clustered") +
  geom_sf(data = samplepoints_clustered, color = "red", size = 0.5) +
  theme

# plot predictors and sampling designs
pdf(file = "code/figures/simulation_study/predictors_and_sample_designs.pdf", width = 12, height = 8)
ggarrange(plot0, grid.arrange(plot1,plot2, nrow=2), ncol = 2, heights = c(0.8), widths = c(735,265), labels = c("(a)", "(b)"))
invisible(dev.off())


################################
## Prepare training data sets ##
################################
# remove response from predictos
predictors <- subset(predictors, predictor_names)

# extract predictors for the sampling locations

# random
trainDat_random <- terra::extract(predictors, samplepoints_random, df=TRUE)
trainDat_random$response <- terra::extract(response, samplepoints_random, ID = FALSE)$response
trainDat_random <- trainDat_random[complete.cases(trainDat_random),]

# clustered
trainDat_clustered <- terra::extract(predictors, samplepoints_clustered, df=TRUE)
trainDat_clustered$response <- terra::extract(response, samplepoints_clustered, ID = FALSE)$response
trainDat_clustered <- trainDat_clustered[complete.cases(trainDat_clustered),]



########################
## Create kNNDM folds ##
########################

fold_colors <- c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695')

set.seed(seed)

# random
knndm_folds_random <- knndm(samplepoints_random, modeldomain = modeldomain, k = 10)

samplepoints_random$fold = knndm_folds_random$clusters

plot_sample_cv_design_random = ggplot() +
  geom_spatvector(data = modeldomain_vect, fill = "white") +
  labs(title = "CV design (random)") +
  geom_sf(mapping = aes(col = factor(fold)), data = samplepoints_random, size = 1) +
  scale_color_manual(name = "knndm folds", values = fold_colors) +
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
  scale_color_manual(name = "knndm folds", values = fold_colors) +
  theme +
  theme(legend.position.inside=c(.8,.16),
        legend.text=element_text(size = 8),
        legend.title=element_text(size = 10),
        legend.background=element_blank(),
        legend.key.height = unit(0.03, "npc"),
        legend.key.width = unit(0.015, "npc"))

pdf(file = "code/figures/simulation_study/sample_designs_knndm_folds.pdf", width = 15, height = 5)
ggarrange(plot_sample_cv_design_random,
          plot_sample_cv_design_clustered,
          labels = c("(a)", "(b)"),
          ncol = 2,
          common.legend = T,
          legend = "bottom")
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

# print variable importances
plot1 = plot(varImp(model_random,scale = F), col="black")
plot2 = plot(varImp(model_clustered,scale = F), col="black")
pdf(file = "code/figures/simulation_study/varImps.pdf", width = 14, height = 6)
ggarrange(plot1, plot2, ncol = 2, nrow = 1, labels = c("(a)", "(b)"))
invisible(dev.off())


# arrange CV errors
cv_errors = rbind(global_validation(model_random),
                  global_validation(model_clustered)) |> 
            as.data.frame() |> 
            mutate("Sampling design" = c("random", "clustered"), .before = 1) |>
            mutate("CV" = c("", ""), .before = 2)


###############################
## Calculate DI, LPD and AOA ##
###############################
# The area of applicability, the dissimilarity index and the local data point
# density are then calculated.

# random
AOA_random <- aoa(newdata = predictors, model = model_random, method = "L2", LPD = TRUE)

# clustered
AOA_clustered <- aoa(newdata = predictors, model = model_clustered, method = "L2", LPD = TRUE)


###########################################
## Prediction and true error calculation ##
###########################################
# The trained models are used to make predictions for the entire target area.
# The true absolute error between prediction and reference is calculated for later
# comparison with the local data point density (LPD).

# response (observation) values
obs <- values(response, na.rm = T)

# random
prediction_random <- predict(predictors, model_random, na.rm = T)
truediff_random <- abs(prediction_random-response)
names(truediff_random) <- "absError"
pred_random <- values(prediction_random, na.rm = T)

# clustered
prediction_clustered <- predict(predictors, model_clustered, na.rm = T)
truediff_clustered <- abs(prediction_clustered-response)
names(truediff_clustered) <- "absError"
pred_clustered <- values(prediction_clustered, na.rm = T)

# arrange true prediction errors
true_errors = rbind(get_errors(pred_random, obs),
                    get_errors(pred_clustered, obs)) |> 
  mutate("Truth" = c("", ""), .before = 0)

# arrange true errors inside AOA
true_errors_inside_AOA = rbind(get_errors(prediction_random[AOA_random$AOA == 1], response[AOA_random$AOA == 1]),
                               get_errors(prediction_clustered[AOA_clustered$AOA == 1], response[AOA_clustered$AOA == 1])) |> 
  mutate("Truth inside AOA" = c("", ""), .before = 0)

errors = cbind(cv_errors, true_errors, true_errors_inside_AOA)

knitr::kable(errors)


#####################################################
## Print parameters of training and prediction LPD ##
#####################################################

# helper functions to arrange parameters
get_parameters_train <- function(aoa){
  data.frame(
     Similarity_threshold = aoa$parameters$threshold[[1]],
     avrgLPD = round(mean(aoa$parameters$trainLPD[aoa$parameters$trainLPD > 0])),
     sd = round(sd(aoa$parameters$trainLPD[aoa$parameters$trainLPD > 0])),
     maxLPD = max(aoa$parameters$trainLPD)
  )
}

get_parameters_pred <- function(aoa){
  data.frame(
    avrgLPD = round(mean(values(aoa$LPD, na.rm = T)[values(aoa$LPD, na.rm = T) > 0])),
    sd = round(sd(values(aoa$LPD, na.rm = T)[values(aoa$LPD, na.rm = T) > 0])),
    maxLPD = max(values(aoa$LPD, na.rm = T))
  )
}

# arrange parameters
train_parameters_random = get_parameters_train(AOA_random)
train_parameters_clustered = get_parameters_train(AOA_clustered)

pred_parameters_random = get_parameters_pred(AOA_random)
pred_parameters_clustered = get_parameters_pred(AOA_clustered)

percentage_aoa_random = round(length(AOA_random$AOA[AOA_random$AOA == 1]) * 100 / length(AOA_random$AOA[!is.na(AOA_random$AOA)]), digits = 1)
percentage_aoa_clustered = round(length(AOA_clustered$AOA[AOA_clustered$AOA == 1]) *100 / length(AOA_clustered$AOA[!is.na(AOA_clustered$AOA)]), digits = 1)

train_parameters = rbind(train_parameters_random,
                        train_parameters_clustered) |> 
                   mutate("Sampling design" = c("random", "clustered"), .before = 1) |>
                   mutate("Training (CV)" = c("", ""), .before = 2)

pred_parameters = rbind(pred_parameters_random,
                      pred_parameters_clustered) |>
                  mutate("Prediction" = c("", ""), .before = 1)

parameters = cbind(train_parameters, pred_parameters)

parameters["% inside AOA"] <- c(percentage_aoa_random, percentage_aoa_clustered)

knitr::kable(parameters)


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


###########################################
## Plot LPD ~ True absolute error (RMSE) ##
###########################################
# visualized via data bins

# random
dat_all_random <- data.frame(absError = values(truediff_random, na.rm=T), 
                             type = rep("LPD",length(values(AOA_random$LPD, na.rm = T))),
                             LPD = values(AOA_random$LPD, na.rm = T),
                             DI = values(AOA_random$DI, na.rm = T))
th_random <- cv_errors$RMSE[1]


plot_random = ggplot(dat_all_random, aes(x = LPD, y = absError)) +
  stat_bin_2d(breaks=list(x = seq(-0.25,max(dat_all_random$LPD)+1,0.5), y = seq(0, ceiling(max(dat_all_random$absError)), 0.01))) +
  ylab("RMSE") +
  xlab("Local data point density (LPD)") +
  scale_fill_gradientn(name = "Data points",
                       trans = "log",
                       breaks = 10^(0:3),
                       colors=viridis(10, alpha=.7)) +
  geom_point(data = LPD_RMSE_random, aes(LPD, RMSE, color = "RMSE"), size = .5, show.legend = F) +
  scale_color_manual(name = "",values = c("RMSE" = "black")) +
  geom_hline(aes(yintercept=th_random,linetype="CV RMSE"), color = "red") +
  scale_linetype_manual(name = "", values = c("CV RMSE" = "solid")) +
  labs(title = "random") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.title = element_text( size = 10),
    legend.position = "none")


# clustered
dat_all_clustered <- data.frame(absError = values(truediff_clustered, na.rm=T), 
                                type = rep("LPD",length(values(AOA_clustered$LPD, na.rm = T))),
                                LPD = values(AOA_clustered$LPD, na.rm = T), 
                                DI = values(AOA_clustered$DI, na.rm = T))
th_clustered <- cv_errors$RMSE[2]


plot_clustered = ggplot(dat_all_clustered, aes(x = LPD, y = absError)) +
  stat_bin_2d(breaks=list(x = seq(-0.25,max(dat_all_clustered$LPD)+1,0.5), y = seq(0, ceiling(max(dat_all_clustered$absError)), 0.01))) +
  ylab("RMSE") +
  xlab("Local data point density (LPD)") +
  scale_fill_gradientn(name = "Data points",
                       trans = "log",
                       breaks = 10^(0:3),
                       colors=viridis(10, alpha=.7)) +
  geom_point(data = LPD_RMSE_clustered, aes(LPD, RMSE, color = "RMSE"), size = .5, show.legend = F) +
  scale_color_manual(name = "",values = c("RMSE" = "black")) +
  geom_hline(aes(yintercept=th_clustered,linetype="CV RMSE"), color = "red") +
  scale_linetype_manual(name = "", values = c("CV RMSE" = "solid")) +
  labs(title = "clustered") +
  ylab("") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.title = element_text( size = 10),
    legend.position = c(0.8, 0.7),
    legend.key.size =unit(0.3, "cm"),
    legend.background = element_rect(fill = "transparent", 
                                     linewidth = 0.5,
                                     linetype="solid")
    )


# generate pdf plot for LPD ~ True abs. error
pdf(file = "code/figures/simulation_study/LPD_true_error.pdf", width = 10, height = 4)
ggarrange(plot_random,
          plot_clustered,
          ncol=2,
          labels = c("(a)", "(b)"))
invisible(dev.off())



################################################
## LPD as a quantitative uncertainty measure? ##
################################################

# compute error models:

# random #
DI_errormodel_random <- errorProfiles(model_random, AOA_random, variable = "DI")
LPD_errormodel_random <- errorProfiles(model_random, AOA_random, variable = "LPD")

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
  scale_linetype_manual(name = "", values = c("model" = "dashed"), guide = "none") +
  scale_color_manual(name = "", values = c("truth" = "red",
                                           "model" = "blue",
                                           "cross-validation" = "black")) +
  labs(title = "random") +
  ylim(c(0, 0.16)) +
  theme_bw() +
  theme(legend.title = element_blank(), 
        legend.position = "none",
        panel.grid = element_blank())
plot_DI_errormodel_random$layers <- plot_DI_errormodel_random$layers[c(1,3,2)]

slidingw$model = predict(DI_errormodel_random, slidingw)
r_squared(slidingw$model, slidingw$true)
rmse(slidingw$model, slidingw$true)

# get residuals
residuals = abs(slidingw["model"]-slidingw["true"])
colnames(residuals) <- "residuals"
residuals$rowIndex <- row.names(residuals)

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
  scale_linetype_manual(name = "", values = c("model" = "dashed"), guide = "none") +
  scale_color_manual(name = "", values = c("truth" = "red",
                                           "model" = "blue",
                                           "cross-validation" = "black")) +
  labs(title = "") +
  ylim(c(0, 0.16)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())
plot_LPD_errormodel_random$layers <- plot_LPD_errormodel_random$layers[c(1,3,2)]

slidingw$model = predict(LPD_errormodel_random, slidingw)
r_squared(slidingw$model, slidingw$true)
rmse(slidingw$model, slidingw$true)

# residual analysis
slidingw$rowIndex = row.names(slidingw)
residuals <- merge(slidingw, residuals, by = "rowIndex")
residuals <- residuals[order(residuals$LPD),]
residuals$what <- "residual"
residuals <- residuals[c("LPD", "residuals", "what")]
cor(residuals$LPD, residuals$residuals, method = "spearman")


plot_residuals_random = ggplot(data = residuals, mapping = aes(x = LPD, y = residuals, shape = what)) + 
                            geom_point() + 
                            # scale_shape_manual(name = "residual", values = c("residual" = 16)) +
                            geom_smooth(method = "gam") +
                            labs(title = "") +
                            ylab("|residual|") +
                            theme_bw() +
                            theme(legend.title = element_blank(),
                                  legend.position = "none",
                                  panel.grid = element_blank())

# clustered
DI_errormodel_clustered <- errorProfiles(model_clustered, AOA_clustered, variable = "DI")
LPD_errormodel_clustered <- errorProfiles(model_clustered, AOA_clustered, variable = "LPD")

# DI error model:
recl <- attr(DI_errormodel_clustered, "performance")

reference_perf <- as.data.frame(list(AOA_clustered$DI,response,prediction_clustered))
reference_perf <- reference_perf[order(reference_perf$DI),]
names(reference_perf) <- c("DI","obs","pred")

# use same moving window over the prediction data to get true RMSE
slidingw <- attr(DI_errormodel_clustered, "performance")
reference_metric <- apply(slidingw,1,function(x){
  x_df <- data.frame(t(x))
  subs_ref <- reference_perf[reference_perf$DI>x_df$ll&
                               reference_perf$DI<x_df$ul,]
  rmse(subs_ref[,"pred"],subs_ref[,"obs"])
})
slidingw$true <- reference_metric
slidingw$what <- "truth"

plot_DI_errormodel_clustered = plot(DI_errormodel_clustered) +
  geom_point(data = slidingw, mapping = aes_string(x = "DI", y = "true", color = "what")) +
  scale_linetype_manual(name = "", values = c("model" = "dashed"), guide = "none") +
  scale_color_manual(name = "", values = c("truth" = "red",
                                           "model" = "blue",
                                           "cross-validation" = "black")) +
  labs(title = "clustered") +
  ylim(c(0, 0.35)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())
plot_DI_errormodel_clustered$layers <- plot_DI_errormodel_clustered$layers[c(1,3,2)]

slidingw$model = predict(DI_errormodel_clustered, slidingw)
r_squared(slidingw$model, slidingw$true)
rmse(slidingw$model, slidingw$true)

# get residuals
residuals = abs(slidingw["model"]-slidingw["true"])
colnames(residuals) <- "residuals"
residuals$rowIndex <- row.names(residuals)


# LPD error model:
recl <- attr(LPD_errormodel_clustered, "performance")

reference_perf <- as.data.frame(list(AOA_clustered$LPD,response,prediction_clustered))
reference_perf <- reference_perf[order(reference_perf$LPD),]
names(reference_perf) <- c("LPD","obs","pred")

# use same moving window over the prediction data to get true RMSE
slidingw <- attr(LPD_errormodel_clustered, "performance")
reference_metric <- apply(slidingw,1,function(x){
  x_df <- data.frame(t(x))
  subs_ref <- reference_perf[reference_perf$LPD>=x_df$ll&
                               reference_perf$LPD<=x_df$ul,]
  rmse(subs_ref[,"pred"],subs_ref[,"obs"])
})
slidingw$true <- reference_metric
slidingw$what <- "truth"

plot_LPD_errormodel_clustered = plot(LPD_errormodel_clustered) +
  geom_point(data = slidingw, mapping = aes_string(x = "LPD", y = "true", color = "what")) +
  scale_linetype_manual(name = "", values = c("model" = "dashed"), guide = "none") +
  scale_color_manual(name = "", values = c("truth" = "red",
                                           "model" = "blue",
                                           "cross-validation" = "black")) +
  labs(title = "") +
  ylim(c(0, 0.35)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())
plot_LPD_errormodel_clustered$layers <- plot_LPD_errormodel_clustered$layers[c(1,3,2)]

slidingw$model = predict(LPD_errormodel_clustered, slidingw)
r_squared(slidingw$model, slidingw$true)
rmse(slidingw$model, slidingw$true)

# residual analysis
slidingw$rowIndex = row.names(slidingw)
residuals <- merge(slidingw, residuals, by = "rowIndex")
residuals <- residuals[order(residuals$LPD),]
residuals$what <- "residual"
residuals <- residuals[c("LPD", "residuals", "what")]
cor(residuals$LPD, residuals$residuals, method = "spearman")

plot_residuals_clustered = ggplot(data = residuals, mapping = aes(x = LPD, y = residuals, shape = what)) + 
  geom_point() + 
  # scale_shape_manual(name = "residual", values = c("residual" = 16)) +
  geom_smooth(method = "gam") +
  labs(title = "") +
  ylab("|residual|") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())


# generate pdf plot for DI and LPD error model
pdf(file = "code/figures/simulation_study/DI_LPD_errormodels.pdf", width = 10, height = 7)
ggarrange(plot_DI_errormodel_random,
          plot_LPD_errormodel_random,
          plot_residuals_random,
          plot_DI_errormodel_clustered,
          plot_LPD_errormodel_clustered,
          plot_residuals_clustered,
          ncol=3,
          nrow = 2,
          #byrow = F,
          legend = "bottom",
          common.legend = T,
          labels = c("(a)", "", "", "(b)", "", ""))
invisible(dev.off())


##################################################
## Estimate model performance with error models ##
##################################################
# prediction only makes sense inside aoa, since model was only fitted for data inside AOA
# inside AOA: DI smaller than AOA threshold and LPD greater than 1

#random
DI_error_prediction_random <- predict(AOA_random$DI, DI_errormodel_random)
LPD_error_prediction_random <- predict(AOA_random$LPD, LPD_errormodel_random)

DI_error_prediction_random[AOA_random$AOA == 0] <- NA
LPD_error_prediction_random[AOA_random$AOA == 0] <- NA

colors <- as.character(values(AOA_random$AOA))
colors[colors==0] <- "violetred" #violetered or darkgoldenrod1
colors[colors==1] <- "transparent"

plot1 = ggplot() +
  geom_spatraster(data = DI_error_prediction_random) +
  scale_fill_viridis_c(na.value = "transparent", limits = c(0, 0.45), option = "viridis", begin = 0.2) +
  geom_spatraster(data = AOA_random$AOA , fill = colors, na.rm = TRUE, show.legend = T) +
  labs(title = "RMSE Prediction (DI)") +
  theme +
  theme(
    legend.position = "none",
    legend.key.width = unit(0.1, "npc")
  )

plot2 = ggplot() +
  geom_spatraster(data = LPD_error_prediction_random) +
  scale_fill_viridis_c(na.value = "transparent", limits = c(0, 0.45), option = "viridis", begin = 0.2) +
  geom_spatraster(data = AOA_random$AOA , fill = colors, na.rm = TRUE, show.legend = T) +
  labs(title = "RMSE Prediction (LPD)") +
  theme +
  theme(
    legend.position = "none",
    legend.key.width = unit(0.1, "npc")
  )

plot3 = ggplot() +
  geom_spatraster(data = truediff_random) +
  scale_fill_viridis_c(na.value = "transparent", limits = c(0, 0.45), option = "viridis", begin = 0.2) +
  geom_spatraster(data = AOA_random$AOA , fill = colors, na.rm = TRUE, show.legend = T) +
  labs(title = "True RMSE") +
  theme +
  theme(
    legend.position = "none",
    legend.key.width = unit(0.1, "npc")
  )

# clustered
DI_error_prediction_clustered <- predict(AOA_clustered$DI, DI_errormodel_clustered)
LPD_error_prediction_clustered <- predict(AOA_clustered$LPD, LPD_errormodel_clustered)

DI_error_prediction_clustered[DI_error_prediction_clustered < 0] <- 0
LPD_error_prediction_clustered[LPD_error_prediction_clustered < 0] <- 0

DI_error_prediction_clustered[AOA_clustered$AOA == 0] <- NA
LPD_error_prediction_clustered[AOA_clustered$AOA == 0] <- NA

colors <- as.character(values(AOA_clustered$AOA))
colors[colors==0] <- "violetred" #violetered or darkgoldenrod1
colors[colors==1] <- "transparent"

plot4 = ggplot() +
  geom_spatraster(data = DI_error_prediction_clustered) +
  scale_fill_viridis_c(na.value = "transparent", limits = c(0, 0.45), option = "viridis", begin = 0.2) +
  geom_spatraster(data = AOA_clustered$AOA , fill = colors, na.rm = TRUE, show.legend = T) +
  labs(title = "RMSE Prediction (DI)") +
  theme +
  theme(
    legend.position = "none",
    legend.key.width = unit(0.1, "npc")
  )

plot5 = ggplot() +
  geom_spatraster(data = LPD_error_prediction_clustered) +
  scale_fill_viridis_c(na.value = "transparent", limits = c(0, 0.45), option = "viridis", begin = 0.2) +
  geom_spatraster(data = AOA_clustered$AOA , fill = colors, na.rm = TRUE, show.legend = T) +
  labs(title = "RMSE Prediction (LPD)") +
  theme +
  theme(
    legend.position = "none",
    legend.key.width = unit(0.1, "npc")
  )

plot6 = ggplot() +
  geom_spatraster(data = truediff_clustered) +
  scale_fill_viridis_c(na.value = "transparent", limits = c(0, 0.45), option = "viridis", begin = 0.2) +
  geom_spatraster(data = AOA_clustered$AOA , fill = colors, na.rm = TRUE, show.legend = T) +
  labs(title = "True RMSE") +
  theme +
  theme(
    legend.position = "none",
    legend.key.width = unit(0.1, "npc")
  )

pdf(file = "code/figures/simulation_study/DI_LPD_error_predictions.pdf", width = 12, height = 10)
ggarrange(plot1,
          plot2,
          plot3,
          plot4,
          plot5,
          plot6,
          ncol=3,
          nrow=2,
          labels = c("(a)  random", "", "", "(b)  clustered", "", ""),
          common.legend = T,
          legend = "bottom",
          # font.label = list(size = 18, color = "black", face = "plain", family = NULL),
          hjust = -0.2)
invisible(dev.off())


#############################################################
## Plot all layers for all sampling scenarios for appendix ##
#############################################################

# random

plot_prediction_random = ggplot() +
  geom_spatraster(data = prediction_random) +
  scale_fill_viridis_c(na.value = "transparent", limits = c(0,1)) +
  labs(title = "Prediction (random)") +
  theme +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.1, "npc")
  )

plot_DI_random <- ggplot() +
  geom_spatraster(data = AOA_random$DI) +
  scale_fill_viridis_c(na.value = "transparent", limits = c(0,5)) +
  labs(title = "DI (random)") +
  theme +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.1, "npc")
  )

plot_LPD_random <- ggplot() +
  geom_spatraster(data = AOA_random$LPD) +
  scale_fill_viridis_c(na.value = "transparent", limits = c(0,60)) +
  labs(title = "LPD (random)") +
  theme +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.1, "npc")
  )

# clustered

plot_prediction_clustered = ggplot() +
  geom_spatraster(data = prediction_clustered) +
  scale_fill_viridis_c(na.value = "transparent", limits = c(0,1)) +
  labs(title = "Prediction (clustered)") +
  theme +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.1, "npc")
  )

plot_DI_clustered <- ggplot() +
  geom_spatraster(data = AOA_clustered$DI) +
  scale_fill_viridis_c(na.value = "transparent", limits = c(0,5)) +
  labs(title = "DI (clustered)") +
  theme +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.1, "npc")
  )

plot_LPD_clustered <- ggplot() +
  geom_spatraster(data = AOA_clustered$LPD) +
  scale_fill_viridis_c(na.value = "transparent", limits = c(0,60)) +
  labs(title = "LPD (clustered)") +
  theme +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.1, "npc")
  )


# create pdf plots
pdf(file = "code/figures/simulation_study/predictions.pdf", width = 14, height = 6)
ggarrange(plot_prediction_random,
          plot_prediction_clustered,
          ncol=2,
          labels = c("(a)", "(b)"),
          common.legend = T,
          legend = "bottom",
          hjust = -2)
invisible(dev.off())

pdf(file = "code/figures/simulation_study/DI.pdf", width = 14, height = 6)
ggarrange(plot_DI_random,
          plot_DI_clustered,
          ncol=2,
          labels = c("(a)", "(b)"),
          common.legend = T,
          legend = "bottom",
          hjust = -2)
invisible(dev.off())

pdf(file = "code/figures/simulation_study/LPD.pdf", width = 14, height = 6)
ggarrange(plot_LPD_random,
          plot_LPD_clustered,
          ncol=2,
          labels = c("(a)", "(b)"),
          common.legend = T,
          legend = "bottom",
          hjust = -2)
invisible(dev.off())
