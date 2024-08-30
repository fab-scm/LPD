################
## Case study ##
################
# This case study is in parts based on the code of the scripts found on
# https://github.com/LOEK-RS/CAST4ecology?tab=readme-ov-file

rm(list=ls())

# for (spatial) data handling
library(terra)
library(sf)
library(tidyverse)

# model training and AOA calculation
library(caret)
library(CAST)
library(ranger)

# visualization
library(viridis)
library(scales)
library(ggpubr)
library(ggplot2)
library(tidyterra)
library(gridExtra)
library(cowplot)
library(reshape2)

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

theme = theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  panel.grid = element_blank(),
  # legend.title = element_blank(),
  plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)),
  strip.text = element_text(size = 14),
  legend.text = element_text(size = 14),
)

#########################
## Data for case study ##
#########################

## Wording:
# training_data: reference samples without coordinates
# predictors: spatially continuous predictor stack of South America
# modeldomain: where we want to predict (the outline of South America)
# predictor_names: names of predictors in the training_data and the predictor stack
# response_name: name of the response variable in plots

#predictors <- rast("data/predictors/predictors_south_america_30s.tif")
predictors <- rast("data/predictors/predictors_south_america_5m.tif")
splotdata <- readRDS("data/samples/training_samples.RDS")
training_data <- splotdata |> st_drop_geometry() # reference samples without coordinates
modeldomain <- st_read("data/modeldomain.gpkg", quiet = TRUE)
rfmodel_ffs <- readRDS("data/model/rfmodel_ffs.rds") # if uncommented the model training can be skipped

predictor_names <- names(predictors)
response_name <- "Species_richness"

# define color palette
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7")


# plot samples with outcome
modeldomain_vect = vect(modeldomain)
plot_samples_with_outcome = ggplot() +
  geom_spatvector(data = modeldomain_vect, fill = "white") +
  labs(title = "Sample locations") +
  geom_sf(mapping = aes(col = Species_richness), data = splotdata, size = 1) +
  scale_color_viridis(name = "Species richness", option = "C", direction = -1, begin = 0.2) +
  theme +
  theme(legend.position.inside=c(.75,.16),
        legend.text=element_text(size = 8),
        legend.title=element_text(size = 10),
        legend.background=element_blank(),
        legend.key.height = unit(0.03, "npc"),
        legend.key.width = unit(0.015, "npc"))

###################################
## Set up kNNDM cross-validation ##
###################################

### Setting up knndm-cv ###
# setting up cv folds such that between-folds distance matches sample-prediction distance
# cv more representative of actual prediction task

knndm_folds = CAST::knndm(tpoints = splotdata,
                          modeldomain = modeldomain, k = 5)

splotdata$fold = knndm_folds$clusters

# plot kNNDM CV folds
plot_sample_cv_design = ggplot() +
  geom_spatvector(data = modeldomain_vect, fill = "white") +
  labs(title = "Sample location CV design") +
  geom_sf(mapping = aes(col = factor(fold)), data = splotdata, size = 1) +
  scale_color_manual(name = "knndm folds", values = Okabe_Ito) +
  theme +
  theme(legend.position=c(.8,.16),
        legend.text=element_text(size = 8),
        legend.title=element_text(size = 10),
        legend.background=element_blank(),
        legend.key.height = unit(0.03, "npc"),
        legend.key.width = unit(0.015, "npc"))

# arrange plots
plot_samples = arrangeGrob(plot_samples_with_outcome, plot_sample_cv_design, nrow=2)

tr_control <- trainControl(method = "cv",
                           number = 5,
                           index = knndm_folds$indx_train,
                           savePredictions = TRUE)

##################################
## RF model training with kNNDM ##
##################################
set.seed(62)
rfmodel_ffs <- CAST::ffs(training_data[,predictor_names],
                         training_data[,response_name],
                         method = "ranger",
                         importance = "permutation",
                         num.trees = 100,
                         trControl = tr_control,
                         verbose = FALSE)

selected_predictors <- rfmodel_ffs$finalModel$xNames
# selected_predictors <- c("bio_1", "bio_8", "bio_12", "bio_14", "elev")

# plot predictors
plot_predictors = ggplot() +
  geom_spatraster(data = stretch(predictors[[selected_predictors]],0,1)) +
  facet_wrap(~lyr, ncol = 3) +
  scale_fill_viridis(na.value = "transparent") +
  theme +
  theme(legend.key.height = unit(0.17, "npc"),
        legend.title = element_blank())

# genrate pdf plot of predictors, sampling locations wit outcome and kNNDM folding
pdf(file = "code/figures/case_study/predictors_samples.pdf", width = 14, height = 9)
plot_grid(plot_predictors, plot_samples, ncol = 2, nrow = 1, rel_widths = c(735,265), labels = c("(a)", "(b)"))
invisible(dev.off())

# cv results
cv_results = rbind(global_validation(rfmodel_ffs)) |>
  as.data.frame() |> mutate("CV" = c("knndm"),
                            "predictors" = c(ncol(rfmodel_ffs$trainingData)-1))
knitr::kable(cv_results)

# plot variable importance
plotVarImp = plot(varImp(rfmodel_ffs,scale = F), col="black")

pdf(file = "code/figures/case_study/varImp.pdf", width = 10, height = 4)
plotVarImp
invisible(dev.off())

# model prediction
prediction_ffs <- predict(predictors, rfmodel_ffs, na.rm = TRUE)



###############################
## Calculate DI, LPD and AOA ##
###############################
# parameters indices = TRUE 
AOA = CAST::aoa(newdata = predictors, model = rfmodel_ffs, method = "L2", LPD = TRUE)


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
train_parameters = get_parameters_train(AOA) 
train_parameters = train_parameters |> mutate("Training (CV)" = c(""), .before = 1)

pred_parameters = get_parameters_pred(AOA)
pred_parameters = pred_parameters |> mutate("Prediction" = c(""), .before = 1)

parameters = cbind(train_parameters, pred_parameters)

knitr::kable(parameters)


# plot DI layer
plot_DI = ggplot() +
  geom_spatraster(data = AOA$DI) +
  scale_fill_viridis_c(name = "DI", na.value = "transparent", option = "viridis", begin = 0.1) +
  labs(title = "DI") +
  theme +
  theme(legend.position=c(.8,.2),
        legend.text=element_text(size = 10),
        legend.title=element_text(size = 14),
        legend.background=element_blank(),
        legend.key.height = unit(0.05, "npc"),
        legend.key.width = unit(0.04, "npc"))

# plot LPD layer
plot_LPD = ggplot() +
  geom_spatraster(data = AOA$LPD) +
  scale_fill_viridis_c(name = "LPD", na.value = "transparent", option = "viridis", begin = 0.1) +
  labs(title = "LPD") +
  theme +
  theme(legend.position=c(.8,.2),
        legend.text=element_text(size = 10),
        legend.title=element_text(size = 14),
        legend.background=element_blank(),
        legend.key.height = unit(0.05, "npc"),
        legend.key.width = unit(0.04, "npc"))


# plot predictions inside AOA

# workaround to mask outside AOA
colors <- as.character(values(AOA$AOA))
colors[colors==0] <- "violetred"
colors[colors==1] <- "transparent"

# workaround for the legend
p = st_point(c(-50,-55.5))
p = st_sfc(p, crs = "epsg:4326")
p = st_sf(p)
p$label = "Outside AOA"


plot_prediction_AOA = ggplot() +
  geom_spatraster(data = prediction_ffs) +
  scale_fill_viridis_c(name = "Predicted\nSpecies Richness", na.value = "transparent", option = "viridis", begin = 0.1) +
  geom_spatraster(data = AOA$AOA , fill = colors, na.rm = TRUE, show.legend = T) +
  geom_sf(data = p, color = c("violetred"), shape = 15, size = 5) +
  geom_sf_text(data = p, aes(label = label), nudge_x = 7) +
  labs(title = "Prediction inside AOA") +
  theme +
  theme(legend.position=c(.78,.22),
        legend.text=element_text(size = 10),
        legend.title=element_text(size = 12),
        legend.background=element_blank(),
        legend.key.height = unit(0.05, "npc"),
        legend.key.width = unit(0.04, "npc"))

# generate PDF with DI, LPD and prediction_AOA plot
pdf(file = "code/figures/case_study/DI_LPD_predAOA.pdf", width = 14, height = 8)
ggarrange(plot_DI,
          plot_LPD,
          plot_prediction_AOA,
          ncol=3,
          labels = c("(a)", "(b)", "(c)"))
invisible(dev.off())


#######################################
## DI and LPD as uncertainty measure ##
#######################################

# calculate error models
DI_errormodel <- errorProfiles(model = rfmodel_ffs, trainDI = AOA$parameters, variable = "DI", calib = "scam")
LPD_errormodel <- errorProfiles(model = rfmodel_ffs, trainDI = AOA$parameters, variable = "LPD", calib = "scam")


# plot errormodels
plot_DI_errormodel = plot(DI_errormodel) + 
                        scale_linetype_manual(name = "", values = c("model" = "solid"), guide = "none") +
                        scale_color_manual(name = "", values = c("cross-validation" = "black",
                                                                 "model" = "blue")) +
                        scale_shape_manual(name = "", values = c("cross-validation" = 16)) +
                        scale_y_continuous(breaks = seq(0, 120, 30)) +
                        ylim(c(0,115)) +
                        theme_bw() +
                        theme(legend.title = element_blank(), legend.position = "none", panel.grid = element_blank())
plot_LPD_errormodel = plot(LPD_errormodel) +
                        scale_linetype_manual(name = "", values = c("model" = "solid"), guide = "none") +
                        scale_color_manual(name = "", values = c("cross-validation" = "black",
                                                                 "model" = "blue")) +
                        scale_shape_manual(name = "", values = c("cross-validation" = 16)) +
                        scale_y_continuous(breaks = seq(0, 120, 30)) +
                        ylim(c(0,115)) +
                        ylab("") +
                        theme_bw() +
                        theme(legend.title = element_blank(), legend.position = "none", panel.grid = element_blank())

# genrate pdf plot for single error models
pdf(file = "code/figures/case_study/errormodel_DI_LPD_complete.pdf", width = 10, height = 4)
ggarrange(plot_DI_errormodel, 
          plot_LPD_errormodel, 
          ncol = 2, 
          labels = c("(a)", "(b)"),
          hjust = 0,
          common.legend = T,
          legend = "bottom"
          )
invisible(dev.off())


# estimate model performance in the target area with the error models
# reliable predictions can only be made inside the AOA as models were fitted there
DI_error_prediction <- predict(AOA$DI, DI_errormodel)
LPD_error_prediction <- predict(AOA$LPD, LPD_errormodel)

# mask predictions outside AOA
DI_error_prediction[AOA$AOA == 0] <- NA
LPD_error_prediction[AOA$AOA == 0] <- NA

# color mask for outside AOA
colors <- as.character(values(AOA$AOA))
colors[colors==0] <- "violetred"
colors[colors==1] <- "transparent"

# work around for outside AOA legend
p = st_point(c(-50,-45.5))
p = st_sfc(p, crs = "epsg:4326")
p = st_sf(p)
p$label = "Outside AOA"

# plot DI error prediction
plot_DI_error_prediction = ggplot() +
  geom_spatraster(data = DI_error_prediction) +
  scale_fill_viridis_c(name = "RMSE", na.value = "transparent", option = "D", limits = c(0,60), breaks = seq(0, 60, 20)) +
  geom_spatraster(data = AOA$AOA , fill = colors, na.rm = TRUE, show.legend = T) +
  # geom_sf(data = p, color = c("violetred"), shape = 15, size = 7) +
  # geom_sf_text(data = p, aes(label = label), nudge_x = 7) +
  labs(title = "RMSE Prediction (DI)") +
  theme +
  theme(legend.position=c(.683,.2),
        legend.text=element_text(size = 10),
        legend.title=element_text(size = 12),
        legend.background=element_blank(),
        legend.key.height = unit(0.05, "npc"),
        legend.key.width = unit(0.04, "npc"))

# plot LPD error prediction
plot_LPD_error_prediction = ggplot() +
  geom_spatraster(data = LPD_error_prediction) +
  scale_fill_viridis_c(name = "RMSE", na.value = "transparent", option = "D", limits = c(0,60), breaks = seq(0, 60, 20)) +
  geom_spatraster(data = AOA$AOA , fill = colors, na.rm = TRUE, show.legend = T) +
  # geom_sf(data = p, color = c("violetred"), shape = 15, size = 7) +
  # geom_sf_text(data = p, aes(label = label), nudge_x = 9) +
  labs(title = "RMSE Prediction (LPD)") +
  theme +
  theme(legend.position=c(.683,.2),
        legend.text=element_text(size = 10),
        legend.title=element_text(size = 12),
        legend.background = element_blank(),
        legend.key.height = unit(0.05, "npc"),
        legend.key.width = unit(0.04, "npc"))


pdf(file = "code/figures/case_study/error_predictions.pdf", width = 12, height = 6)
ggarrange(plot_DI_error_prediction,
          plot_LPD_error_prediction,
          ncol=2,
          labels = c("(a)", "(b)"),
          common.legend = T,
          legend = "right")
invisible(dev.off())

