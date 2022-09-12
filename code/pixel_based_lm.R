###########################################################
####################gpp seasonlaity regression#############
###########################################################
library('plyr');library('tidyverse');library('caret');library('leaps');library('raster')
source('R:/cssp_brazil/cssp_brazil_R/code/all_extract_gpp_seasonality_functions.R')
pred_var <- c('LAI','CICA','GSDSR','APAR','ET')
extractfun <- function(m) {
  tinfo <- summary(m)$sigma
  r2 <- summary(m)$adj.r.squared
  data.frame(Rsq = r2,Resid.se = tinfo)
}
amazon_nw_poly <- shapefile("./data/amazon_nw.shp");amazon_sw_poly <- shapefile("./data/amazon_sw.shp");amazon_ec_poly <- shapefile("./data/amazon_ec.shp");amazon_bs_poly <- shapefile("./data/amazon_bs.shp");amazon_gs_poly <- shapefile("./data/amazon_gs.shp")
####extract data from raster into dataframe
amazonia_regres_data <- reg_card_var_regression(reg_data[[6]],reg_names[6],var_names)
summary(amazonia_regres_data)
# all_amazonia_region_regres_data <- na.omit(all_reg_var_regres_run(reg_data,reg_names,var_names))
# summary(all_amazonia_region_regres_data)
# rownames(amazonia_regres_data)<-1:nrow(amazonia_regres_data)
# rownames(all_amazonia_region_regres_data)<-1:nrow(all_amazonia_region_regres_data)

####run linear regression on all cells
# split_pixels_lm <- lapply(split(amazonia_regres_data,amazonia_regres_data$cell_no),function(s){lm(gpp~lai, data=s)})
# split_pixels_lm <- lm(gpp~lai*cell_no-1,data=amazonia_regres_data)
#split_pixels_lm <- dlply(amazonia_regres_data, .(cell_no), function(s) lm(gpp~lai+cica+gsdsr+apar+et, data=s))

split_pixels_lm_gpp_lai <- dlply(na.omit(amazonia_regres_data), .(cell_no,x,y), function(s) lm(gpp~lai, data=s))
summary(split_pixels_lm_gpp_lai)
summary(split_pixels_lm$C0592);names(summary(split_pixels_lm$C0592))
summary(split_pixels_lm$C0592)$sigma
cv(split_pixels_lm$C0592)

fit <- lm(gpp~lai, data=na.omit(amazonia_regres_data))
cv(fit)
var_name <- 'test'
assign(var_name,1:3)
test

# run_all_var<- function(data,variables){
#   for (i in variables){
#     split_data_lm <- dlply(na.omit(data), .(cell_no,x,y), function(s) lm(gpp~i, data=s))
#       lm_summary<-ldply(split_data_lm, extractfun)
#       pixels_lm_rsq_raster <- rasterFromXYZ(lm_summary[,2:4])
#       pixels_lm_rse_raster <- rasterFromXYZ(lm_summary[,c(2:3,5)])
#       crs(pixels_lm_rsq_raster) <- '+proj=longlat +datum=WGS84 +no_defs'
#       crs(pixels_lm_rse_raster) <- '+proj=longlat +datum=WGS84 +no_defs'
#       }
#   par(mfrow = c(1,2))
#   plot(pixels_lm_rsq_raster);plot(pixels_lm_rse_raster)
# }

run_all_var<- function(data,variables){
  for (i in variables){
    if (variables=='LAI'){
      split_data_lm <- dlply(na.omit(data), .(cell_no,x,y), function(s) lm(gpp~lai, data=s))
      gpp_lai_lm_summary<-ldply(split_data_lm, extractfun)
      pixels_gpp_lai_lm_rsq_raster <- rasterFromXYZ(gpp_lai_lm_summary[,2:4])
      pixels_gpp_lai_lm_rse_raster <- rasterFromXYZ(gpp_lai_lm_summary[,c(2:3,5)])
      crs(pixels_gpp_lai_lm_rsq_raster) <- '+proj=longlat +datum=WGS84 +no_defs'
      crs(pixels_gpp_lai_lm_rse_raster) <- '+proj=longlat +datum=WGS84 +no_defs'
      par(mfrow = c(1,2),oma = c(0, 0, 2, 0))
      plot(pixels_gpp_lai_lm_rsq_raster,main='R sq');plot(amazon_nw_poly,add=T);plot(amazon_sw_poly,add=T);plot(amazon_ec_poly,add=T);plot(amazon_bs_poly,add=T);plot(amazon_gs_poly,add=T)
      plot(pixels_gpp_lai_lm_rse_raster,main='Residual SE');plot(amazon_nw_poly,add=T);plot(amazon_sw_poly,add=T);plot(amazon_ec_poly,add=T);plot(amazon_bs_poly,add=T);plot(amazon_gs_poly,add=T)
      mtext(paste(i), line=-1, outer=TRUE, cex=1.5)
    }
    else if (variables=='CICA'){
      split_data_lm <- dlply(na.omit(data), .(cell_no,x,y), function(s) lm(gpp~cica, data=s))
      gpp_cica_lm_summary<-ldply(split_data_lm, extractfun)
      pixels_gpp_cica_lm_rsq_raster <- rasterFromXYZ(gpp_cica_lm_summary[,2:4])
      pixels_gpp_cica_lm_rse_raster <- rasterFromXYZ(gpp_cica_lm_summary[,c(2:3,5)])
      crs(pixels_gpp_cica_lm_rsq_raster) <- '+proj=longlat +datum=WGS84 +no_defs'
      crs(pixels_gpp_cica_lm_rse_raster) <- '+proj=longlat +datum=WGS84 +no_defs'
      par(mfrow = c(1,2),oma = c(0, 0, 2, 0))
      plot(pixels_gpp_cica_lm_rsq_raster,main='R sq');plot(amazon_nw_poly,add=T);plot(amazon_sw_poly,add=T);plot(amazon_ec_poly,add=T);plot(amazon_bs_poly,add=T);plot(amazon_gs_poly,add=T)
      plot(pixels_gpp_cica_lm_rse_raster,main='Residual SE');plot(amazon_nw_poly,add=T);plot(amazon_sw_poly,add=T);plot(amazon_ec_poly,add=T);plot(amazon_bs_poly,add=T);plot(amazon_gs_poly,add=T)
      mtext(paste(i), line=-1, outer=TRUE, cex=1.5)
    }
    else if (variables=='GSDSR'){
      split_data_lm <- dlply(na.omit(data), .(cell_no,x,y), function(s) lm(gpp~gsdsr, data=s))
      gpp_gsdsr_lm_summary<-ldply(split_data_lm, extractfun)
      pixels_gpp_gsdsr_lm_rsq_raster <- rasterFromXYZ(gpp_gsdsr_lm_summary[,2:4])
      pixels_gpp_gsdsr_lm_rse_raster <- rasterFromXYZ(gpp_gsdsr_lm_summary[,c(2:3,5)])
      crs(pixels_gpp_gsdsr_lm_rsq_raster) <- '+proj=longlat +datum=WGS84 +no_defs'
      crs(pixels_gpp_gsdsr_lm_rse_raster) <- '+proj=longlat +datum=WGS84 +no_defs'
      par(mfrow = c(1,2),oma = c(0, 0, 2, 0))
      plot(pixels_gpp_gsdsr_lm_rsq_raster,main='R sq');plot(amazon_nw_poly,add=T);plot(amazon_sw_poly,add=T);plot(amazon_ec_poly,add=T);plot(amazon_bs_poly,add=T);plot(amazon_gs_poly,add=T)
      plot(pixels_gpp_gsdsr_lm_rse_raster,main='Residual SE');plot(amazon_nw_poly,add=T);plot(amazon_sw_poly,add=T);plot(amazon_ec_poly,add=T);plot(amazon_bs_poly,add=T);plot(amazon_gs_poly,add=T)
      mtext(paste(i), line=-1, outer=TRUE, cex=1.5)
    }
    else if (variables=='APAR'){
      split_data_lm <- dlply(na.omit(data), .(cell_no,x,y), function(s) lm(gpp~apar, data=s))
      gpp_apar_lm_summary<-ldply(split_data_lm, extractfun)
      pixels_gpp_apar_lm_rsq_raster <- rasterFromXYZ(gpp_apar_lm_summary[,2:4])
      pixels_gpp_apar_lm_rse_raster <- rasterFromXYZ(gpp_apar_lm_summary[,c(2:3,5)])
      crs(pixels_gpp_apar_lm_rsq_raster) <- '+proj=longlat +datum=WGS84 +no_defs'
      crs(pixels_gpp_apar_lm_rse_raster) <- '+proj=longlat +datum=WGS84 +no_defs'
      par(mfrow = c(1,2),oma = c(0, 0, 2, 0))
      plot(pixels_gpp_apar_lm_rsq_raster,main='R sq');plot(amazon_nw_poly,add=T);plot(amazon_sw_poly,add=T);plot(amazon_ec_poly,add=T);plot(amazon_bs_poly,add=T);plot(amazon_gs_poly,add=T)
      plot(pixels_gpp_apar_lm_rse_raster,main='Residual SE');plot(amazon_nw_poly,add=T);plot(amazon_sw_poly,add=T);plot(amazon_ec_poly,add=T);plot(amazon_bs_poly,add=T);plot(amazon_gs_poly,add=T)
      mtext(paste(i), line=-1, outer=TRUE, cex=1.5)
    }
    else if (variables=='ET'){
      split_data_lm <- dlply(na.omit(data), .(cell_no,x,y), function(s) lm(gpp~et, data=s))
      gpp_et_lm_summary<-ldply(split_data_lm, extractfun)
      pixels_gpp_et_lm_rsq_raster <- rasterFromXYZ(gpp_et_lm_summary[,2:4])
      pixels_gpp_et_lm_rse_raster <- rasterFromXYZ(gpp_et_lm_summary[,c(2:3,5)])
      crs(pixels_gpp_et_lm_rsq_raster) <- '+proj=longlat +datum=WGS84 +no_defs'
      crs(pixels_gpp_et_lm_rse_raster) <- '+proj=longlat +datum=WGS84 +no_defs'
      par(mfrow = c(1,2),oma = c(0, 0, 2, 0))
      plot(pixels_gpp_et_lm_rsq_raster,main='R sq');plot(amazon_nw_poly,add=T);plot(amazon_sw_poly,add=T);plot(amazon_ec_poly,add=T);plot(amazon_bs_poly,add=T);plot(amazon_gs_poly,add=T)
      plot(pixels_gpp_et_lm_rse_raster,main='Residual SE');plot(amazon_nw_poly,add=T);plot(amazon_sw_poly,add=T);plot(amazon_ec_poly,add=T);plot(amazon_bs_poly,add=T);plot(amazon_gs_poly,add=T)
      mtext(paste(i), line=-1, outer=TRUE, cex=1.5)
    }
  }
}
for (i in pred_var){
  run_all_var(amazonia_regres_data,i)
}

# split_pixels_lm <- dlply(all_amazonia_region_regres_data, .(cell_no,region_name), function(s) lm(gpp~lai, data=s))
# str(split_pixels_lm)
# split_pixels_lm$C0592

# set.seed(123) 
# train.control <- trainControl(method = "cv", number = 10)
# 
# split_pixels_lm_cv <- dlply(all_amazonia_region_regres_data, .(cell_no,region_name), function(s) train(gpp ~lai+cica+apar, data = s,
#                                                                                                        method = "leapSeq",
#                                                                                                        tuneGrid = data.frame(nvmax = 1:3),
#                                                                                                        trControl = train.control))
# 
# names(split_pixels_lm_cv)
# split_pixels_lm_cv$C0115.amazon_nw$bestTune$nvmax
# split_pixels_lm_cv$C0115.amazon_nw$results[split_pixels_lm_cv$C0115.amazon_nw$results$nvmax==split_pixels_lm_cv$C0115.amazon_nw$bestTune$nvmax,]
# 
# split_pixels_lm_cv$C0115.amazon_nw$finalModel
# split_pixels_lm_cv$C0115.amazon_nw$coefnames
# names(coef(split_pixels_lm_cv$C0115.amazon_nw$finalModel,split_pixels_lm_cv$C0115.amazon_nw$bestTune$nvmax))
# 
# names(split_pixels_lm_cv$C0115.amazon_nw)

# Take a list (of models) as input and output a data frame:
pixels_gpp_lai_lm_summary<-ldply(split_pixels_lm_gpp_lai, extractfun)

# extractfun_cv <- function(m) {
#   best_coef_no <- m$bestTune
#   best_result <- m$result[m$result$nvmax==best_coef_no[1,],]
#   
#   data.frame(no_of_coef = best_coef_no[1,], no_of_coef_factor = as.factor(best_coef_no[1,]), Rsq = best_result[,3],RMSE=best_result[,2])
# }

# # Take a list (of models) as input and output a data frame:
# summary_pixel_lm_cv_step <- ldply(split_pixels_lm_cv, extractfun_cv)
# summary(summary_pixel_lm_cv_step)
# 
# range(summary_pixel_lm_cv_step[summary_pixel_lm_cv_step$region_name=='amazon_nw',]$Rsq)
# range(summary_pixel_lm_cv_step[summary_pixel_lm_cv_step$region_name=='amazon_sw',]$Rsq)
# range(summary_pixel_lm_cv_step[summary_pixel_lm_cv_step$region_name=='amazon_ec',]$Rsq)
# range(summary_pixel_lm_cv_step[summary_pixel_lm_cv_step$region_name=='amazon_gs',]$Rsq)
# range(summary_pixel_lm_cv_step[summary_pixel_lm_cv_step$region_name=='amazon_bs',]$Rsq)
# 
# summary(summary_pixel_lm_cv_step[summary_pixel_lm_cv_step$region_name=='amazon_nw',])
# summary(summary_pixel_lm_cv_step[summary_pixel_lm_cv_step$region_name=='amazon_sw',])
# summary(summary_pixel_lm_cv_step[summary_pixel_lm_cv_step$region_name=='amazon_ec',])
# summary(summary_pixel_lm_cv_step[summary_pixel_lm_cv_step$region_name=='amazon_gs',])
# summary(summary_pixel_lm_cv_step[summary_pixel_lm_cv_step$region_name=='amazon_bs',])
# 
# # Define training control
# set.seed(123) 
# train.control <- trainControl(method = "cv", number = 10)
# # Train the model
# gpp_model_test <- train(gpp ~., data = amazonia_regres_data[,c(3:5,7)], method = "lm",
#                trControl = train.control)
# # Summarize the results
# summary(gpp_model_test)
# names(gpp_model_test)
# gpp_model_test$coefnames
# 
# # Train the model
# gpp_step.model_test <- train(gpp ~., data = amazonia_regres_data[,c(3:8)],
#                              method = "leapSeq",
#                              tuneGrid = data.frame(nvmax = 1:3),
#                              trControl = train.control
#                              )
# summary(gpp_step.model_test)
# names(gpp_step.model_test)
# gpp_step.model_test$results
# gpp_step.model_test$bestTune
# gpp_step.model_test$finalModel
# gpp_step.model_test$method
# gpp_step.model_test$modelInfo
# gpp_step.model_test$modelType
# gpp_step.model_test$results
# gpp_step.model_test$metric
# gpp_step.model_test$coefnames
# gpp_step.model_test$call
# gpp_step.model_test$coefnames
# gpp_step.model_test$coefnames
# 
# 
#   
# summary(gpp_step.model_test$finalModel)

########################################
###########plot lm results##############
########################################
# summary(amazonia_regres_data);str(amazonia_regres_data)
# trial_regres_raster<-amazonia_regres_data[amazonia_regres_data$date=='X2001.01.01',c(1,2,5)]
# dfr <- rasterFromXYZ(na.omit(trial_regres_raster))
# extent(dfr) <- extent(biomass_amazon_gCm2)#Convert first two columns as lon-lat and third as value                
# crs(dfr) <- '+proj=longlat +datum=WGS84 +no_defs'
# par(mfrow = c(1,1))
# plot(dfr)

pixels_gpp_lai_lm_rsq_raster <- rasterFromXYZ(pixels_gpp_lai_lm_summary[,2:4])
pixels_gpp_lai_lm_rse_raster <- rasterFromXYZ(pixels_gpp_lai_lm_summary[,c(2:3,5)])

crs(pixels_gpp_lai_lm_rsq_raster) <- '+proj=longlat +datum=WGS84 +no_defs'
crs(pixels_gpp_lai_lm_rse_raster) <- '+proj=longlat +datum=WGS84 +no_defs'

par(mfrow = c(1,2))

plot(pixels_gpp_lai_lm_rsq_raster);plot(pixels_gpp_lai_lm_rse_raster)
