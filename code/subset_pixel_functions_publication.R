#All codes being Used
#install packages----
library(ncdf4); library(raster); library(dplyr); library(ggplot2);library(ggpubr);library(quantreg);library(ggpp);library(rgeos);library(ggpmisc);library(rgdal);library(Metrics)

#done
#data prep pre-cursors----
# amazonia_subset_initial <- shapefile("R:/cssp_brazil/cssp_brazil_R/data/amazonia_subset.shp")
# amazonia_subset <- shapefile("R:/cssp_brazil/cssp_brazil_R/data/amazonia_subset_ifl.shp")
amazonia_ifl_layer <- shapefile("R:/brazil_leeds_maps/ifl_2000_amazonia.shp")
compare_var_names <- c('WOOD','NPP_wood_flx','OUTPUT_wood_flx','MTT_wood')
plot(amazonia_ifl_layer)
# mod_var <- c('esa_cci_agb','rainfor_biomass_productivity_2005','rainfor_biomass_annual_productivity')
# time_period <- c('2000-2009','2010-2016')
prefix <- 'R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/amazonia_'
suffix <- '_nomngt_2000_updated_2010.nc'
# # suffix_old <- '_1deg_monthly_2001_updated_2019.nc'
# # rt_amazon_00_09_day<-biomass_amazon_gCm2/biommort_00_09_gCm2d
# # rt_amazon_00_09_year<-rt_amazon_00_09_day/365.25
# # rt_amazon_10_16_day<-biomass_amazon_gCm2/biommort_10_16_gCm2d
# # rt_amazon_10_16_year<-rt_amazon_10_16_day/365.25
# 
# # plot(biomass_amazon_gCm2,main='Woody Biomass t/ha');plot(amazon_nw_poly,add=T);plot(amazon_sw_poly,add=T);plot(amazon_ec_poly,add=T);plot(amazon_bs_poly,add=T);plot(amazon_gs_poly,add=T);plot(amazonia_subset,add=T)
# # 
# # plot(woodprod_00_09_gCm2d,main='Woody Productivity 2000-2009 t/ha/year');plot(amazon_nw_poly,add=T);plot(amazon_sw_poly,add=T);plot(amazon_ec_poly,add=T);plot(amazon_bs_poly,add=T);plot(amazon_gs_poly,add=T);plot(amazonia_subset,add=T)
# # plot(woodprod_10_16_gCm2d,main='Woody Productivity 2010-2016 t/ha/year');plot(amazon_nw_poly,add=T);plot(amazon_sw_poly,add=T);plot(amazon_ec_poly,add=T);plot(amazon_bs_poly,add=T);plot(amazon_gs_poly,add=T);plot(amazonia_subset,add=T)
# # plot(biommort_00_09_gCm2d,main='Woody Mortality 2000-2009 t/ha/year');plot(amazon_nw_poly,add=T);plot(amazon_sw_poly,add=T);plot(amazon_ec_poly,add=T);plot(amazon_bs_poly,add=T);plot(amazon_gs_poly,add=T);plot(amazonia_subset,add=T)
# # plot(biommort_10_16_gCm2d,main='Woody Mortality 2010-2016 t/ha/year');plot(amazon_nw_poly,add=T);plot(amazon_sw_poly,add=T);plot(amazon_ec_poly,add=T);plot(amazon_bs_poly,add=T);plot(amazon_gs_poly,add=T);plot(amazonia_subset,add=T)
# # plot(rt_amazon_00_09_year,main='Woody Residence Time 2000-2009 years');plot(amazon_nw_poly,add=T);plot(amazon_sw_poly,add=T);plot(amazon_ec_poly,add=T);plot(amazon_bs_poly,add=T);plot(amazon_gs_poly,add=T);plot(amazonia_subset,add=T)
# # plot(rt_amazon_10_16_year,main='Woody Residence Time 2010-2016 years');plot(amazon_nw_poly,add=T);plot(amazon_sw_poly,add=T);plot(amazon_ec_poly,add=T);plot(amazon_bs_poly,add=T);plot(amazon_gs_poly,add=T);plot(amazonia_subset,add=T)
# 
# 
# reference_data<- c(biomass_amazon_gCm2,woodprod_00_09_gCm2d,woodprod_10_16_gCm2d,biommort_00_09_gCm2d,biommort_10_16_gCm2d,rt_amazon_00_09_year,rt_amazon_10_16_year)
# reference_data<- c(biomass_amazon_gCm2,woodprod_00_09_gCm2d,biommort_00_09_gCm2d,rt_amazon_00_09_year)

########################################################################
####new code to include variable, model and region######################
########################################################################
#extract subset of wood, npp and mortality from rainfor benchmark data
extract_subset <- function (region,reference){
  masked1 <- mask(reference, region)
  masked2 <- masked1 > -Inf
  poly <- rasterToPolygons(masked2, dissolve=TRUE)
  data_region <- mask(reference,poly)
  return(data_region)
}
biomass_ifl_subset<-extract_subset(amazonia_ifl_layer,biomass_amazon_gCm2)
# plot(biomass_ifl_subset)
# writeRaster(biomass_ifl_subset, "./data/biomass_ifl_subset.tif")
            
###################################################
################Woody Biomass######################
###################################################
#extract subset of wood, npp and mortality from models
extract_biomass_subset_all_models <- function (region,cardamom_var,model_variant,reference){
  for (i in model_variant) {
    data <- brick(paste(prefix,model_variant[1],suffix,sep=""),varname=cardamom_var)
    # data<-data[[which(getZ(data) >= as.Date("2000-01-01") & getZ(data) <= as.Date("2010-12-01"))]]
    data_mean <- stackApply(data, indices =  rep(1,nlayers(data)), fun = "mean")
    masked1 <- mask(data_mean, region)
    masked2 <- masked1 > -Inf
    poly <- rasterToPolygons(masked2, dissolve=TRUE)
    data_region <- mask(data_mean,poly)
    }
  return(data_region)
}

collate_all_models_cwood <- function (region,cardamom_var,model_variant,reference){
  for (i in model_variant) {
    wood_biomass_raster <- extract_biomass_subset_all_models(region,cardamom_var,i,reference)
    return(wood_biomass_raster)
  }
}

res_df_merge_plot_cwood <- function (region,cardamom_var,model_variant,reference){
  benchmark <- extract_subset(region,reference)
  df_benchmark <- as.data.frame(benchmark, xy=TRUE)
  for (i in model_variant) {
        if (i=='esa_cci_agb') {
          model <- extract_biomass_subset_all_models(region,cardamom_var,i,reference)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$esa_cci_agb, model_benchmark_df$benchmark, main="",xlab="ESA CCI Biomass only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=2.5, cex.axis=2.5)+
                   abline(coef = c(0,1),col='red', lwd=5)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''), cex=2.5)
                   )
        }
        else if (i=='rainfor_biomass_productivity_2005') {
          model <- extract_biomass_subset_all_models(region,cardamom_var,i,reference)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$rainfor_biomass_productivity_2005, model_benchmark_df$benchmark, main="",xlab="2005 Biomass and Productivity only", ylab="", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=2.5, cex.axis=2.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=2.5))
        }
        else if (i=='rainfor_biomass_annual_productivity') {
          model <- extract_biomass_subset_all_models(region,cardamom_var,i,reference)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$rainfor_biomass_annual_productivity, model_benchmark_df$benchmark, main="",xlab="Biomass annual and Productivity only", ylab="", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=2.5, cex.axis=2.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=2.5))
        }
      }
      par(mfrow = c(1,3),oma = c(0, 0, 3, 0))
      plot_esa_cci_agb
      plot_rainfor_biomass_productivity_2005
      plot_rainfor_biomass_annual_productivity
      mtext(bquote("Woody Biomass 2000-2010 ("~ g~ m^-2~")"), line=-2, side=3, outer=TRUE, cex=2.5)
}

# extract_biomass_subset_all_models(amazonia_subset,compare_var_names[1],mod_var[1],reference_data[[1]])
# collate_all_models_cwood(amazonia_subset,compare_var_names[1],mod_var[1],reference_data[[1]])
# par(mfrow = c(1,3));res_df_merge_plot_cwood(amazonia_subset,compare_var_names[1],mod_var,reference_data[[1]])

#done
#######################################################
################NPP and mortality######################
#######################################################
#extract subset of wood, npp and mortality from models
extract_npp_mort_subset_models <- function (region,cardamom_var,model_variant,reference){
  for (i in model_variant) {
    data <- brick(paste(prefix,model_variant[1],suffix,sep=""),varname=cardamom_var)
    # data<-data[[which(getZ(data) >= as.Date("2000-01-01") & getZ(data) <= as.Date("2010-12-01"))]]
    data_mean <- stackApply(data, indices =  rep(1,nlayers(data)), fun = "mean")
    masked1 <- mask(data_mean, region)
    masked2 <- masked1 > -Inf
    poly <- rasterToPolygons(masked2, dissolve=TRUE)
    data_region <- mask(data_mean,poly)
    }
  return(data_region)
}

collate_all_models_npp_mort <- function (region,cardamom_var,model_variant,reference){
  for (i in model_variant) {
    wood_biomass_raster <- extract_biomass_subset_all_models(region,cardamom_var,i,reference)
    return(wood_biomass_raster)
  }
}

#extract and plot for each variable
res_df_merge_plot_npp_mort <- function (region,cardamom_var,model_variant,reference){
  benchmark <- extract_subset(region,reference)
  df_benchmark <- as.data.frame(benchmark, xy=TRUE)
  for (j in cardamom_var) {
    if (j =='NPP_wood_flx') {
      for (i in model_variant) {
        if (i=='esa_cci_agb') {
          model <- extract_npp_mort_subset_models(region,j,i,reference)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$esa_cci_agb, model_benchmark_df$benchmark, main="",xlab="ESA CCI Biomass only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=2.5, cex.axis=2.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=2.5))
        }
        else if (i=='rainfor_biomass_productivity_2005') {
          model <- extract_npp_mort_subset_models(region,j,i,reference)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$rainfor_biomass_productivity_2005, model_benchmark_df$benchmark, main="",xlab="2005 Biomass and Productivity only", ylab="", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=2.5, cex.axis=2.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=2.5))
        }
        else if (i=='rainfor_biomass_annual_productivity') {
          model <- extract_npp_mort_subset_models(region,j,i,reference)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$rainfor_biomass_annual_productivity, model_benchmark_df$benchmark, main="",xlab="Biomass annual and Productivity only", ylab="", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=2.5, cex.axis=2.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=2.5))
        }
      }
      par(mfrow = c(1,3),oma = c(0, 0, 3, 0))
      plot_esa_cci_agb
      plot_rainfor_biomass_productivity_2005
      plot_rainfor_biomass_annual_productivity
      mtext(bquote("Woody Productivity 2000-2010 ("~ g~ m^-2~day^-1~")"), line=-2, side=3, outer=TRUE, cex=2.5)
    }
    else if (j == 'OUTPUT_wood_flx'){
      for (i in model_variant) {
        if (i=='esa_cci_agb') {
          model <- extract_npp_mort_subset_models(region,j,i,reference)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$esa_cci_agb, model_benchmark_df$benchmark, main="",xlab="ESA CCI Biomass only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=2.5, cex.axis=2.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=2.5))
        }
        else if (i=='rainfor_biomass_productivity_2005') {
          model <- extract_npp_mort_subset_models(region,j,i,reference)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$rainfor_biomass_productivity_2005, model_benchmark_df$benchmark, main="",xlab="2005 Biomass and Productivity only", ylab="", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=2.5, cex.axis=2.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=2.5))
        }
        else if (i=='rainfor_biomass_annual_productivity') {
          model <- extract_npp_mort_subset_models(region,j,i,reference)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$rainfor_biomass_annual_productivity, model_benchmark_df$benchmark, main="",xlab="Biomass annual and Productivity only", ylab="", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=2.5, cex.axis=2.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=2.5))
        }
      }
      par(mfrow = c(1,3),oma = c(0, 0, 3, 0))
      plot_esa_cci_agb
      plot_rainfor_biomass_productivity_2005
      plot_rainfor_biomass_annual_productivity
      mtext(bquote("Woody Mortality 2000-2010 ("~ g~ m^-2~day^-1~")"), line=-2, side=3, outer=TRUE, cex=2.5)
    }
  }
  }
# expression("Woody Productivity 2000-2010 g m"^-2 "day"^-1
# extract_npp_mort_subset_models(amazonia_subset,compare_var_names[2],mod_var[1],woodprod_00_09_gCm2d)
# collate_all_models_npp_mort(amazonia_subset,compare_var_names[2],mod_var[1],reference_data[[3]])
# par(mfrow = c(1,3));res_df_merge_plot_npp_mort(amazonia_subset,compare_var_names[2],mod_var,reference_data[[2]])
# par(mfrow = c(1,3));res_df_merge_plot_npp_mort(amazonia_subset,compare_var_names[3],mod_var,reference_data[[3]])

#done
####################################################
################Residence Time######################
####################################################
res_df_merge_plot_rt <- function (region,cardamom_var,model_variant,reference){
  benchmark <- extract_subset(region,reference)
  df_benchmark <- as.data.frame(benchmark, xy=TRUE)
  for (i in model_variant) {
        if (i=='esa_cci_agb') {
          model <- extract_npp_mort_subset_models(region,cardamom_var,i,reference)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$esa_cci_agb, model_benchmark_df$benchmark, main="",xlab="ESA CCI Biomass only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=2.5, cex.axis=2.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=2.5))
        }
        else if (i=='rainfor_biomass_productivity_2005') {
          model <- extract_npp_mort_subset_models(region,cardamom_var,i,reference)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$rainfor_biomass_productivity_2005, model_benchmark_df$benchmark, main="",xlab="2005 Biomass and Productivity only", ylab="", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=2.5, cex.axis=2.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=2.5))
        }
        else if (i=='rainfor_biomass_annual_productivity') {
          model <- extract_npp_mort_subset_models(region,cardamom_var,i,reference)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$rainfor_biomass_annual_productivity, model_benchmark_df$benchmark, main="",xlab="Biomass annual and Productivity only", ylab="", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=2.5, cex.axis=2.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=2.5))
        }
      }
  par(mfrow = c(1,3),oma = c(0, 0, 3, 0))
  plot_esa_cci_agb
  plot_rainfor_biomass_productivity_2005
  plot_rainfor_biomass_annual_productivity
  mtext(bquote("Woody Mean Residence Time 2000-2010 ("~ years~")"), line=-2, side=3, outer=TRUE, cex=2.5)
  }

# par(mfrow = c(1,3));res_df_merge_plot_rt(amazonia_subset,compare_var_names[4],mod_var,reference_data[[4]])

#done
#######################################################
################Combine all variables##################
#######################################################
one_function_to_plot_them_all <- function(region,cardamom_var,model_variant,reference) {
  for (i in cardamom_var) {
    if (cardamom_var=='WOOD') {
      res_df_merge_plot_cwood(region,i,model_variant,reference[[1]])
      }
    else if (cardamom_var=='NPP_wood_flx') {
      res_df_merge_plot_npp_mort(region,i,model_variant,reference[[2]])
      }
    else if (cardamom_var=='OUTPUT_wood_flx') {
      res_df_merge_plot_npp_mort(region,i,model_variant,reference[[3]])
      }
    else if (cardamom_var=='MTT_wood') {
      res_df_merge_plot_rt(region,i,model_variant,reference[[4]])
      }
  }
}

par(mfrow = c(1,3));res_df_merge_plot_cwood(amazonia_subset,compare_var_names[1],mod_var,reference_data[[1]])
par(mfrow = c(1,3));res_df_merge_plot_npp_mort(amazonia_subset,compare_var_names[2],mod_var,reference_data[[2]])
par(mfrow = c(1,3));res_df_merge_plot_npp_mort(amazonia_subset,compare_var_names[3],mod_var,reference_data[[3]])
par(mfrow = c(1,3));res_df_merge_plot_rt(amazonia_subset,compare_var_names[4],mod_var,reference_data[[4]])
# done

#######################################################
###############extract model parameters################
#######################################################
# prefix_new <- 'M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/amazonia_ifl_'
# midfix <- '_nomngt/RESULTS_PROCESSED/amazonia_ifl_'
# suffix_new <- '_nomngt_stock_flux.RData'
# 
# mod_var_new <- c('esa_cci_agb','rainfor_biomass_annual','rainfor_biomass_productivity_2005','rainfor_biomass_annual_productivity')
# mod_var_abr <- c('ECB','B_a','B_P_2005','B_a_P')

model_pars_list <- list()
# model_pars_list_sd <- list()  

extract_parameters_all_models <- function (model_variant,model_variant_abr){
  for (i in seq_along(model_variant)) {
    load(paste(prefix_new,model_variant[i],midfix,model_variant[i],suffix_new,sep=""))
    data <- grid_output$parameters[,,,7]
    data_par_mean <- colMeans(data, dims = 2,na.rm=T)
    data_par_mean[34] <- model_variant_abr[i]
    model_pars_list[[i]]<-data_par_mean
  }
  return(model_pars_list)
}
plot_all_parameters <- function(data,parameter_names) {
  for (i in 2:ncol(data)){
    print(ggplot(data, aes(x = model_var,y=data[,i])) +
            geom_bar(stat = "identity", width=0.5)+
            labs(title=parameter_names[i], y = parameter_names[i])+
            scale_fill_grey(start = 0, end = 1) +
            theme_bw()+
            theme(text = element_text(size = 16),panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")))
  }
}
round_df <- function(x, digits) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}


# extract_parameters_all_models_sd <- function (model_variant,model_variant_abr){
#   for (i in seq_along(model_variant)) {
#     load(paste(prefix_new,model_variant[i],midfix,model_variant[i],suffix_new,sep=""))
#     data <- grid_output$parameters[,,,4]
#     # for (j in seq_along(data[37,32,])){
#     #   data_par_sd <- sd(data[37,32,j],na.rm=T)
#     # }
#     # data_par_sd[34] <- model_variant_abr[i]
#     # model_pars_list_sd[[i]]<-data_par_sd
#     data_par_sd <- sd(data, dims = 2,na.rm=T)
#     data_par_sd[34] <- model_variant_abr[i]
#     model_pars_list[[i]]<-data_par_mean
#   }
#   return(model_pars_list_sd)
# }

# model_pars_df_sd<-as.data.frame(t(as.data.frame(extract_parameters_all_models_sd(mod_var_new,mod_var_abr))), stringsAsFactors = FALSE)
# model_pars_df_sd[] <- lapply(model_pars_df_sd, type.convert, as.is = TRUE)

# par_names <- c('Lit2SOM (day_1)','GPP%Ra','NPP_fol_frac','NPP_root_frac','Leaf lifespan','TO Wood','TO Roots','Mineralise Lit','SOM2Rh',
#                'Temp fac','Canopy Eff (gC/m2leaf/day)','Max bud burst day','NPP_lab_frac','Labile_release_period','max leaf fall',
#                'leaf fall period','LMA (gC/m2)','Labile_initial','Foliar_initial','Roots_initial','Wood_initial','Lit_initial',
#                'SOM_initial','Soil Water_initial','Cwood_coarseR_frac','CoarseR Biomass 50% Max Depth','Max Root Depth',
#                'Res factor Burned Cstocks','CCF Fol','CCF Wood and FineR','CCF Soil','CCF Fol and FineR Lit','unknown','model_var')
# # summary(model_pars_df)
# 
# param_names <- c('model_var','Lit2SOM (day_1)','GPP%Ra','NPP_fol_frac','NPP_root_frac','Leaf lifespan','TO Wood','TO Roots','Mineralise Lit','SOM2Rh',
#                'Temp fac','Canopy Eff (gC/m2leaf/day)','Max bud burst day','NPP_lab_frac','Labile_release_period','max leaf fall',
#                'leaf fall period','LMA (gC/m2)','Labile_initial','Foliar_initial','Roots_initial','Wood_initial','Lit_initial',
#                'SOM_initial','Soil Water_initial','Cwood_coarseR_frac','CoarseR Biomass 50% Max Depth','Max Root Depth',
#                'Res factor Burned Cstocks','CCF Fol','CCF Wood and FineR','CCF Soil','CCF Fol and FineR Lit')

model_pars_df<-as.data.frame(t(as.data.frame(extract_parameters_all_models(mod_var_new,mod_var_abr))), stringsAsFactors = FALSE)
model_pars_df[] <- lapply(model_pars_df, type.convert, as.is = TRUE)
rownames(model_pars_df)<-  seq_along(model_pars_df[,1])
colnames(model_pars_df)<-  par_names
model_pars_df$model_var<-as.factor(model_pars_df$model_var)
model_pars_df<-model_pars_df[,c(34,1:32)]

model_pars_df$MRT_wood<-(1/model_pars_df[,7])/365.25
model_pars_df$MRT_roots<-(1/model_pars_df[,8])/365.25
model_pars_df$NPP_wood_frac<- 1-rowSums(model_pars_df[,c(4,5,14)])

new_model_pars_df<-round_df(model_pars_df[,2:36], 3)
rownames(new_model_pars_df)<-mod_var_new
print(t(new_model_pars_df))

plot_all_parameters(model_pars_df,param_names)

# write.csv(t(new_model_pars_df),"rainfor_cardamom_analysis_model_paramaters_ifl_nomngt_updated.csv")
# write.csv(t(new_model_pars_df),"rainfor_cardamom_analysis_model_paramaters_ifl_nomngt_updated_2pt5.csv")
write.csv(t(new_model_pars_df),"rainfor_cardamom_analysis_model_paramaters_ifl_nomngt_updated_97pt5.csv")
