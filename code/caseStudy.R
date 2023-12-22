################
## Case study ##
################
# This case study is in parts based on the code of the scripts found on
# https://github.com/LOEK-RS/CAST4ecology?tab=readme-ov-file

# for (spatial) data handling
library(terra)
library(sf)
library(tidyverse)

# model training and AOA calculation
remotes::install_github("fab-scm/CAST")
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

# uncomment if you want to inspect the calcaulted layers with the explore AOA functions
# packages are not installed with the CAST packages, so make sure you install them first
# # exploreAOA
# library(leaflet)
# library(shiny)
# library(shinycssloaders)
# library(rlist)
# library(bslib)
# library(plotly)
# library(shinyWidgets)


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
# rfmodel_ffs <- readRDS("data/model/rfmodel_ffs.rds") # if uncommented the model training can be skipped

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
  theme(legend.position=c(.75,.16),
        legend.text=element_text(size = 8),
        legend.title=element_text(size = 10),
        legend.background=element_blank(),
        legend.key.height = unit(0.03, "npc"),
        legend.key.width = unit(0.015, "npc"))
plot_samples_with_outcome

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
plot_sample_cv_design

# arrange plots
plot_samples = arrangeGrob(plot_samples_with_outcome,plot_sample_cv_design, nrow=2)

# plot geodistance
gd_knndm <- geodist(splotdata,
                    modeldomain,
                    cvfolds = splotdata$fold)

plot(gd_knndm, stat = "density") +
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0),
                     trans = "log10",
                     labels = trans_format("log10", math_format(10^.x)))+
  theme_light()+
  theme(legend.position = "bottom")

plot(gd_knndm, stat = "ecdf")+
  scale_x_continuous(limits = c(0,1500000))+
  theme_light()+
  theme(legend.position = "bottom")


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
selected_predictors <- c("bio_1", "bio_8", "bio_12", "bio_14", "elev")

# plot predictors
plot_predictors = ggplot() +
  geom_spatraster(data = stretch(predictors[[selected_predictors]],0,1)) +
  facet_wrap(~lyr, ncol = 3) +
  scale_fill_viridis(na.value = "transparent") +
  theme +
  theme(legend.key.height = unit(0.17, "npc"),
        legend.title = element_blank())

# genrate pdf plot of predictors, sampling locations wit outcome and kNNDM folding
pdf(file = "analysis/figures/case_study/predictors_samples.pdf", width = 14, height = 9)
plot_grid(plot_predictors, plot_samples, ncol = 2, nrow = 1, rel_widths = c(735,265), labels = c("(a)", "(b)"))
invisible(dev.off())

# cv results
cv_results = rbind(global_validation(rfmodel_ffs)) |>
  as.data.frame() |> mutate("CV" = c("knndm"),
                            "predictors" = c(ncol(rfmodel_ffs$trainingData)-1))
knitr::kable(cv_results)

# plot ffs results and variable importance
plot_ffs = plot(rfmodel_ffs, plotType = "selected") +
  theme_light()+
  theme(panel.grid.major.y = element_blank(),
        axis.text = element_text(color = "black", size = 10))

plotVarImp = plot(varImp(rfmodel_ffs,scale = F), col="black")

pdf(file = "analysis/figures/case_study/FFS_varImp.pdf", width = 10, height = 5)
plot_grid(plot_ffs, plotVarImp, ncol = 2, labels = c("(a)", "(b)"))
invisible(dev.off())

# plot model prediction
prediction_ffs <- predict(predictors, rfmodel_ffs, na.rm = TRUE)
plot_prediction = ggplot() +
  geom_spatraster(data = prediction_ffs) +
  scale_fill_viridis(name = "Species richness", na.value = "transparent", option = "mako", begin = 0.2) +
  labs(title = "Prediction") +
  theme +
  theme(legend.key.height = unit(0.16, "npc")
  )
plot_prediction



###############################
## Calculate DI, LPD and AOA ##
###############################
AOA = CAST::aoa(newdata = predictors, model = rfmodel_ffs, method = "L2", LPD = TRUE, maxLPD = 1)
plot(AOA)[1]
plot(AOA)[2]

# find sampling locations under data/samples/training_samples: to add them to the map
exploreAOA(AOA)



# plot relationship between DI and LPD
df = data.frame(LPD = values(AOA$LPD, na.rm = TRUE), DI = round(values(AOA$DI, na.rm = TRUE),digits = 3))
df[df$LPD == 0 & df$DI <= 0.34, "DI"] = 0.341 # shift values for clear bins
plotDILPD <- ggplot(df, aes(x = LPD, y = DI)) +
  stat_bin_2d(breaks = list(x = seq(-0.5, max(df$LPD)+1, 1), y = seq(0, ceiling(max(df$DI)), 0.01))) +
  scale_fill_viridis(begin = 0.1) +
  geom_hline(aes(yintercept = AOA$parameters$threshold, linetype = "AOA_threshold")) +
  scale_linetype_manual(name = "", values = c(AOA_threshold = "dashed")) +
  labs(title = "LPD ~ DI (random)") +
  theme_bw() +
  geom_smooth(data = df, aes(LPD, DI, color = "LPD~DI"), method = "gam", se = FALSE) +
  scale_color_manual(name = "",values = c("LPD~DI" = "red"))
plotDILPD # takes a moment due to gam model fitting

pdf(file = "analysis/figures/case_study/DI_LPD_relationship.pdf", width = 10, height = 5)
plotDILPD # takes a moment due to gam model fitting
invisible(dev.off())

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
plot_DI

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
plot_LPD


# plot predictions inside AOA

# workaround to mask outside AOA
colors <- as.character(values(AOA$AOA))
colors[colors==0] <- "darkgoldenrod1"
colors[colors==1] <- "transparent"

# workaround for the legend
p = st_point(c(-50,-55.5))
p = st_sfc(p, crs = "epsg:4326")
p = st_sf(p)
p$label = "Outside AOA"


plot_prediction_AOA = ggplot() +
  geom_spatraster(data = prediction_ffs) +
  scale_fill_viridis_c(name = "Predicted\nSpecies Richness", na.value = "transparent", option = "mako", begin = 0.1) +
  geom_spatraster(data = AOA$AOA , fill = colors, na.rm = TRUE, show.legend = T) +
  geom_sf(data = p, color = c("darkgoldenrod1"), shape = 15, size = 5) +
  geom_sf_text(data = p, aes(label = label), nudge_x = 7) +
  labs(title = "Prediction inside AOA") +
  theme +
  theme(legend.position=c(.78,.22),
        legend.text=element_text(size = 10),
        legend.title=element_text(size = 12),
        legend.background=element_blank(),
        legend.key.height = unit(0.05, "npc"),
        legend.key.width = unit(0.04, "npc"))
plot_prediction_AOA

# generate PDF with DI, LPD and prediction_AOA plot
pdf(file = "analysis/figures/case_study/DI_LPD_predAOA.pdf", width = 14, height = 8)
plot_grid(plot_DI,
          plot_LPD,
          plot_prediction_AOA,
          ncol=3,
          labels = c("(a)", "(b)", "(c)"))
invisible(dev.off())


#######################################
## DI and LPD as uncertainty measure ##
#######################################

# calculate error models
DI_errormodel <- DItoErrormetric(model = rfmodel_ffs, trainDI = AOA$parameters, calib = "scam")
LPD_errormodel <- LPDtoErrormetric(model = rfmodel_ffs, trainDI = AOA$parameters, calib = "exp")
DI_LPD_errormodel <- DI_LPDtoErrormetric(model = rfmodel_ffs, trainDI = AOA$parameters, calib = "scam_exp")

# plot errormodels
plot_DI_errormodel = plot(DI_errormodel) + labs(title = "DI ~ metric (RMSE)")
plot_DI_errormodel
plot_LPD_errormodel = plot(LPD_errormodel) + labs(title = "LPD ~ metric (RMSE)")
plot_LPD_errormodel
plot_DI_LPD_errormodel = plot(DI_LPD_errormodel)
plot_DI_LPD_errormodel %>%
  layout(title = list(text = "DI + LPD ~ metric (RMSE)", x = 0.2, y = 0.8))

# genrate pdf plot for single error models
pdf(file = "analysis/figures/case_study/errormodel_DI_LPD_complete.pdf", width = 10, height = 5)
plot_grid(plot_DI_errormodel, plot_LPD_errormodel, ncol = 2, labels = c("(a)", "(b)"))
invisible(dev.off())

# png plot of combined error model was made screenshotted

# estimate model performance in the target area with the error models
# reliable predictions can only be made inside the AOA as models were fitted there
DI_error_prediction <- predict(AOA$DI, DI_errormodel)
LPD_error_prediction <- predict(AOA$LPD, LPD_errormodel)
DI_LPD_error_prediction <- predict(rast(list(AOA$DI, AOA$LPD)), DI_LPD_errormodel)

# mask predictions outside AOA
DI_error_prediction[AOA$AOA == 0] <- NA
LPD_error_prediction[AOA$AOA == 0] <- NA
DI_LPD_error_prediction[AOA$AOA == 0] <- NA

# color mask for outside AOA
colors <- as.character(values(AOA$AOA))
colors[colors==0] <- "violetred"
colors[colors==1] <- "transparent"

# work around for outside AOA legend
p = st_point(c(-50,-55.5))
p = st_sfc(p, crs = "epsg:4326")
p = st_sf(p)
p$label = "Outside AOA"

# plot DI error prediction
plot_DI_error_prediction = ggplot() +
  geom_spatraster(data = DI_error_prediction) +
  scale_fill_viridis_c(name = "RMSE", na.value = "transparent", option = "D") +
  geom_spatraster(data = AOA$AOA , fill = colors, na.rm = TRUE, show.legend = T) +
  geom_sf(data = p, color = c("violetred"), shape = 15, size = 7) +
  geom_sf_text(data = p, aes(label = label), nudge_x = 7) +
  labs(title = "DI predicted RMSE") +
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
  scale_fill_viridis_c(name = "RMSE", na.value = "transparent", option = "D") +
  geom_spatraster(data = AOA$AOA , fill = colors, na.rm = TRUE, show.legend = T) +
  geom_sf(data = p, color = c("violetred"), shape = 15, size = 7) +
  geom_sf_text(data = p, aes(label = label), nudge_x = 7) +
  labs(title = "LPD predicted RMSE") +
  theme +
  theme(legend.position=c(.683,.2),
        legend.text=element_text(size = 10),
        legend.title=element_text(size = 12),
        legend.background=element_blank(),
        legend.key.height = unit(0.05, "npc"),
        legend.key.width = unit(0.04, "npc"))

# plot DI+LPD error prediction
plot_DI_LPD_error_prediction = ggplot() +
  geom_spatraster(data = DI_LPD_error_prediction) +
  scale_fill_viridis_c(name = "RMSE", na.value = "transparent", option = "D") +
  geom_spatraster(data = AOA$AOA , fill = colors, na.rm = TRUE, show.legend = T) +
  geom_sf(data = p, color = c("violetred"), shape = 15, size = 7) +
  geom_sf_text(data = p, aes(label = label), nudge_x = 7) +
  labs(title = "DI and LPD predicted RMSE") +
  theme +
  theme(legend.position=c(.683,.2),
        legend.text=element_text(size = 10),
        legend.title=element_text(size = 12),
        legend.background=element_blank(),
        legend.key.height = unit(0.05, "npc"),
        legend.key.width = unit(0.04, "npc"))


pdf(file = "analysis/figures/case_study/error_predictions.pdf", width = 14, height = 8)
plot_grid(plot_DI_error_prediction,
          plot_LPD_error_prediction,
          plot_DI_LPD_error_prediction,
          ncol=3,
          labels = c("(a)", "(b)", "(c)"))
invisible(dev.off())

