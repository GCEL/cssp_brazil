#All codes being Used
#install packages----
# library(ncdf4); library(raster); library(dplyr); library(ggplot2);library(ggpubr);library(quantreg);library(ggpp);library(rgeos);library(ggpmisc);library(rgdal);library(Metrics)

#done
#data prep pre-cursors----
# amazonia_subset <- shapefile("R:/cssp_brazil/cssp_brazil_R/data/amazonia_subset.shp")
# compare_var_names <- c('WOOD','NPP_wood_flx','OUTPUT_wood_flx','MTT_wood')
# mod_var <- c('esa_cci_agb_nbe','esa_cci_biomass_only',
#              'biomass_2010','biomass_initial','biomass_annual',
#              'mortality_only','productivity_only',
#              'biomass_2010_productivity','biomass_2010_mortality',
#              'biomass_initial_productivity','biomass_annual_productivity',
#              'biomass_initial_mortality','biomass_annual_mortality',
#              'productivity_mortality_only',
#              'biomass_2010_productivity_mortality','biomass_initial_productivity_mortality','biomass_annual_productivity_mortality'
#              )
# time_period <- c('2000-2009','2010-2016')
# prefix <- 'R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/DATA/benchmark/Amazon_subset_'
# suffix <- '_2000_updated_2019.nc'
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

###################################################
################Woody Biomass######################
###################################################
#extract subset of wood, npp and mortality from models
extract_biomass_subset_all_models <- function (region,cardamom_var,model_variant,reference,tp){
  for (i in model_variant) {
    data <- brick(paste(prefix,model_variant[1],suffix,sep=""),varname=cardamom_var)
    data_0109<-data[[which(getZ(data) >= as.Date("2000-01-01") & getZ(data) <= as.Date("2009-12-01"))]]
    data_1016<-data[[which(getZ(data) >= as.Date("2010-01-01") & getZ(data) <= as.Date("2016-12-01"))]]
    data_0109_mean <- stackApply(data_0109, indices =  rep(1,nlayers(data_0109)), fun = "mean")
    data_1016_mean <- stackApply(data_1016, indices =  rep(1,nlayers(data_1016)), fun = "mean")
  }
  return(c(data_0109_mean,data_1016_mean))
}

collate_all_models_cwood <- function (region,cardamom_var,model_variant,reference,tp){
  for (i in model_variant) {
    wood_biomass_raster <- extract_biomass_subset_all_models(region,cardamom_var,i,reference,j)
  }
  for (j in tp) {
    if (j=='2000-2009'){
      wood_biomass_raster <-wood_biomass_raster[[1]]
      return(wood_biomass_raster)
    }
    else if (j=='2010-2016') {
      wood_biomass_raster <-wood_biomass_raster[[2]]
      return(wood_biomass_raster)
    }
  }
}

res_df_merge_plot_cwood <- function (region,cardamom_var,model_variant,reference,tp){
  for (j in tp) {
    if (j=='2000-2009') {
      benchmark <- extract_subset(region,reference)
      df_benchmark <- as.data.frame(benchmark, xy=TRUE)
      for (i in model_variant) {
        if (i=='esa_cci_biomass_only') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$esa_cci_biomass_only, model_benchmark_df$benchmark, main="",xlab="ESA CCI Biomass only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''), cex=1.5))
        }
        else if (i=='esa_cci_agb_nbe') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$esa_cci_agb_nbe, model_benchmark_df$benchmark, main="",xlab="ESA CCI Biomass and NBE", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''), cex=1.5))
        }
        else if (i=='biomass_2010') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_2010, model_benchmark_df$benchmark, main="",xlab="2010 Biomass only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_initial') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_initial, model_benchmark_df$benchmark, main="",xlab="Initial Biomass only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_annual') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_annual, model_benchmark_df$benchmark, main="",xlab="Annual Biomass only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='productivity_only') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$productivity_only, model_benchmark_df$benchmark, main="",xlab="Productivity only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='mortality_only') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$mortality_only, model_benchmark_df$benchmark, main="",xlab="Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_2010_mortality') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_2010_mortality, model_benchmark_df$benchmark, main="",xlab="2010 Biomass and Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_2010_productivity') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_2010_productivity, model_benchmark_df$benchmark, main="",xlab="2010 Biomass and Productivity only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_initial_mortality') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_initial_mortality, model_benchmark_df$benchmark, main="",xlab="Biomass Initial and Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_initial_productivity') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_initial_productivity, model_benchmark_df$benchmark, main="",xlab="Biomass Initial and Productivity only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_annual_mortality') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_annual_mortality, model_benchmark_df$benchmark, main="",xlab="Biomass annual and Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_annual_productivity') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_annual_productivity, model_benchmark_df$benchmark, main="",xlab="Biomass annual and Productivity only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_annual_productivity_nbe') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_annual_productivity_nbe, model_benchmark_df$benchmark, main="",xlab="Biomass annual and Productivity only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='productivity_mortality_only') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$productivity_mortality_only, model_benchmark_df$benchmark, main="",xlab="Productivity and Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_2010_productivity_mortality') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_2010_productivity_mortality, model_benchmark_df$benchmark, main="",xlab="2010 Bio Prod and Mort only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_initial_productivity_mortality') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_initial_productivity_mortality, model_benchmark_df$benchmark, main="",xlab="Bio Initial Prod and Mort only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_annual_productivity_mortality') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_annual_productivity_mortality, model_benchmark_df$benchmark, main="",xlab="Bio Annual Prod and Mort only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
      }
      par(mfrow = c(5,5),oma = c(0, 0, 2, 0))
      plot_esa_cci_biomass_only
      plot_esa_cci_agb_nbe
      plot_biomass_2010
      plot_biomass_initial
      plot_biomass_annual
      plot_productivity_only
      plot_mortality_only
      plot_biomass_2010_mortality
      plot_biomass_2010_productivity
      plot_biomass_initial_mortality
      plot_biomass_initial_productivity
      plot_biomass_annual_mortality
      plot_biomass_annual_productivity
      plot_biomass_annual_productivity_nbe
      plot_productivity_mortality_only
      # plot_biomass_2010_productivity_mortality
      plot_biomass_initial_productivity_mortality
      # plot_biomass_annual_productivity_mortality
      mtext(expression("Woody Biomass 2000-2009 g m"^-2), line=-2, side=3, outer=TRUE, cex=1.5)
    }
    else if (j=='2010-2016') {
      benchmark <- extract_subset(region,reference)
      df_benchmark <- as.data.frame(benchmark, xy=TRUE)
      for (i in model_variant) {
        if (i=='esa_cci_biomass_only') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$esa_cci_biomass_only, model_benchmark_df$benchmark, main="",xlab="ESA CCI Biomass only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''), cex=1.5))
        }
        else if (i=='esa_cci_agb_nbe') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$esa_cci_agb_nbe, model_benchmark_df$benchmark, main="",xlab="ESA CCI Biomass and NBE", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''), cex=1.5))
        }
        else if (i=='biomass_2010') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_2010, model_benchmark_df$benchmark, main="",xlab="2010 Biomass only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_initial') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_initial, model_benchmark_df$benchmark, main="",xlab="Initial Biomass only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_annual') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_annual, model_benchmark_df$benchmark, main="",xlab="Annual Biomass only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='productivity_only') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$productivity_only, model_benchmark_df$benchmark, main="",xlab="Productivity only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='mortality_only') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$mortality_only, model_benchmark_df$benchmark, main="",xlab="Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_2010_mortality') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_2010_mortality, model_benchmark_df$benchmark, main="",xlab="2010 Biomass and Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_2010_productivity') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_2010_productivity, model_benchmark_df$benchmark, main="",xlab="2010 Biomass and Productivity only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_initial_mortality') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_initial_mortality, model_benchmark_df$benchmark, main="",xlab="Biomass Initial and Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_initial_productivity') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_initial_productivity, model_benchmark_df$benchmark, main="",xlab="Biomass Initial and Productivity only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_annual_mortality') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_annual_mortality, model_benchmark_df$benchmark, main="",xlab="Biomass annual and Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_annual_productivity') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_annual_productivity, model_benchmark_df$benchmark, main="",xlab="Biomass annual and Productivity only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_annual_productivity_nbe') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_annual_productivity_nbe, model_benchmark_df$benchmark, main="",xlab="Biomass annual and Productivity only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='productivity_mortality_only') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$productivity_mortality_only, model_benchmark_df$benchmark, main="",xlab="Productivity and Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_2010_productivity_mortality') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_2010_productivity_mortality, model_benchmark_df$benchmark, main="",xlab="Bio Prod and Mort only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_2010_productivity_mortality') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_2010_productivity_mortality, model_benchmark_df$benchmark, main="",xlab="2010 Bio Prod and Mort only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_initial_productivity_mortality') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_initial_productivity_mortality, model_benchmark_df$benchmark, main="",xlab="Bio Initial Prod and Mort only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_annual_productivity_mortality') {
          model <- collate_all_models_cwood(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 0)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_annual_productivity_mortality, model_benchmark_df$benchmark, main="",xlab="Bio Annual Prod and Mort only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(20000, 2000, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
      }
      par(mfrow = c(5,5),oma = c(0, 0, 2, 0))
      plot_esa_cci_biomass_only
      plot_esa_cci_agb_nbe
      plot_biomass_2010
      plot_biomass_initial
      plot_biomass_annual
      plot_productivity_only
      plot_mortality_only
      plot_biomass_2010_mortality
      plot_biomass_2010_productivity
      plot_biomass_initial_mortality
      plot_biomass_initial_productivity
      plot_biomass_annual_mortality
      plot_biomass_annual_productivity
      plot_biomass_annual_productivity_nbe
      plot_productivity_mortality_only
      # plot_biomass_2010_productivity_mortality
      plot_biomass_initial_productivity_mortality
      # plot_biomass_annual_productivity_mortality
      mtext(expression("Woody Biomass 2010-2016 g m"^-2), line=-2, side=3, outer=TRUE, cex=1.5)
    }
  }
}
# extract_biomass_subset_all_models(amazonia_subset,compare_var_names[1],mod_var[1],biomass_amazon_gCm2)
# collate_all_models_cwood(amazonia_subset,compare_var_names[1],mod_var[1],biomass_amazon_gCm2)
# par(mfrow = c(2,4));res_df_merge_plot_cwood(amazonia_subset,compare_var_names[1],mod_var,biomass_amazon_gCm2)

#done
#######################################################
################NPP and mortality######################
#######################################################
#extract subset of wood, npp and mortality from models
extract_npp_mort_subset_models <- function (region,cardamom_var,model_variant,reference,tp){
  for (i in model_variant) {
    data <- brick(paste(prefix,model_variant[1],suffix,sep=""),varname=cardamom_var)
    data_0109<-data[[which(getZ(data) >= as.Date("2000-01-01") & getZ(data) <= as.Date("2009-12-01"))]]
    data_1016<-data[[which(getZ(data) >= as.Date("2010-01-01") & getZ(data) <= as.Date("2016-12-01"))]]
    data_0109_mean <- stackApply(data_0109, indices =  rep(1,nlayers(data_0109)), fun = "mean")
    data_1016_mean <- stackApply(data_1016, indices =  rep(1,nlayers(data_1016)), fun = "mean")
  }
  return(c(data_0109_mean,data_1016_mean))
}

collate_all_models_npp_mort <- function (region,cardamom_var,model_variant,reference,tp){
  for (i in model_variant) {
    wood_npp_mort_raster <- extract_npp_mort_subset_models(region,cardamom_var,i,reference,j)
  }
  for (j in tp) {
    if (j=='2000-2009'){
      wood_npp_mort_raster <-wood_npp_mort_raster[[1]]
      return(wood_npp_mort_raster)
    }
    else if (j=='2010-2016') {
      wood_npp_mort_raster <-wood_npp_mort_raster[[2]]
      return(wood_npp_mort_raster)
    }
  }
}

#extract and plot for each variable
res_df_merge_plot_npp_mort <- function (region,cardamom_var,model_variant,reference,tp){
  for (j in tp) {
    if (j=='2000-2009') {
      benchmark <- extract_subset(region,reference)
      df_benchmark <- as.data.frame(benchmark, xy=TRUE)
      for (i in model_variant) {
        if (i=='esa_cci_biomass_only') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$esa_cci_biomass_only, model_benchmark_df$benchmark, main="",xlab="ESA CCI Biomass only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='esa_cci_agb_nbe') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$esa_cci_agb_nbe, model_benchmark_df$benchmark, main="",xlab="ESA CCI Biomass and NBE", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_2010') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_2010, model_benchmark_df$benchmark, main="",xlab="2010 Biomass only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_initial') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_initial, model_benchmark_df$benchmark, main="",xlab="Initial Biomass Only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_annual') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_annual, model_benchmark_df$benchmark, main="",xlab="Annual Biomass Only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='productivity_only') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df[,4], model_benchmark_df[,3], main="",xlab="Productivity only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='mortality_only') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df[,4], model_benchmark_df[,3], main="",xlab="Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_2010_mortality') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_2010_mortality, model_benchmark_df$benchmark, main="",xlab="2010 Biomass and Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_2010_productivity') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_2010_productivity, model_benchmark_df$benchmark, main="",xlab="2010 Biomass and Productivity only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_initial_mortality') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_initial_mortality, model_benchmark_df$benchmark, main="",xlab="Biomass Initial and Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_initial_productivity') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_initial_productivity, model_benchmark_df$benchmark, main="",xlab="Biomass Initial and Productivity only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_annual_mortality') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_annual_mortality, model_benchmark_df$benchmark, main="",xlab="Biomass annual and Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_annual_productivity') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_annual_productivity, model_benchmark_df$benchmark, main="",xlab="Biomass annual and Productivity only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_annual_productivity_nbe') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_annual_productivity_nbe, model_benchmark_df$benchmark, main="",xlab="Biomass annual and Productivity only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='productivity_mortality_only') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df[,4], model_benchmark_df[,3], main="",xlab="Productivity and Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_2010_productivity_mortality') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_2010_productivity_mortality, model_benchmark_df$benchmark, main="",xlab="2010 Bio Prod and Mort only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_initial_productivity_mortality') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_initial_productivity_mortality, model_benchmark_df$benchmark, main="",xlab="Bio Initial Prod and Mort only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_annual_productivity_mortality') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_annual_productivity_mortality, model_benchmark_df$benchmark, main="",xlab="Bio Annual Prod and Mort only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
      }
      par(mfrow = c(5,5),oma = c(0, 0, 2, 0))
      plot_esa_cci_biomass_only
      plot_esa_cci_agb_nbe
      plot_biomass_2010
      plot_biomass_initial
      plot_biomass_annual
      plot_productivity_only
      plot_mortality_only
      plot_biomass_2010_mortality
      plot_biomass_2010_productivity
      plot_biomass_initial_mortality
      plot_biomass_initial_productivity
      plot_biomass_annual_mortality
      plot_biomass_annual_productivity
      plot_biomass_annual_productivity_nbe
      plot_productivity_mortality_only
      # plot_biomass_2010_productivity_mortality
      plot_biomass_initial_productivity_mortality
      # plot_biomass_annual_productivity_mortality
      mtext(bquote(.(cardamom_var)~.(j)~g~ m^-2~day^-1), line=-2, side=3, outer=TRUE, cex=1.5)
    }
    else if (j=='2010-2016') {
      benchmark <- extract_subset(region,reference)
      df_benchmark <- as.data.frame(benchmark, xy=TRUE)
      for (i in model_variant) {
        if (i=='esa_cci_biomass_only') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$esa_cci_biomass_only, model_benchmark_df$benchmark, main="",xlab="ESA CCI Biomass only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='esa_cci_agb_nbe') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$esa_cci_agb_nbe, model_benchmark_df$benchmark, main="",xlab="ESA CCI Biomass and NBE", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_2010') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_2010, model_benchmark_df$benchmark, main="",xlab="2010 Biomass Only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_initial') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_initial, model_benchmark_df$benchmark, main="",xlab="Initial Biomass Only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_annual') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_annual, model_benchmark_df$benchmark, main="",xlab="Annual Biomass Only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='productivity_only') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df[,4], model_benchmark_df[,3], main="",xlab="Productivity only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='mortality_only') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df[,4], model_benchmark_df[,3], main="",xlab="Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_2010_mortality') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_2010_mortality, model_benchmark_df$benchmark, main="",xlab="2010 Biomass and Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_2010_productivity') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_2010_productivity, model_benchmark_df$benchmark, main="",xlab="2010 Biomass and Productivity only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_initial_mortality') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_initial_mortality, model_benchmark_df$benchmark, main="",xlab="Biomass Initial and Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_initial_productivity') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_initial_productivity, model_benchmark_df$benchmark, main="",xlab="Biomass Initial and Productivity only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_annual_mortality') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_annual_mortality, model_benchmark_df$benchmark, main="",xlab="Biomass annual and Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_annual_productivity') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_annual_productivity, model_benchmark_df$benchmark, main="",xlab="Biomass annual and Productivity only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_annual_productivity_nbe') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_annual_productivity_nbe, model_benchmark_df$benchmark, main="",xlab="Biomass annual and Productivity only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='productivity_mortality_only') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df[,4], model_benchmark_df[,3], main="",xlab="Productivity and Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_2010_productivity_mortality') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_2010_productivity_mortality, model_benchmark_df$benchmark, main="",xlab="2010 Bio Prod and Mort only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_initial_productivity_mortality') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_initial_productivity_mortality, model_benchmark_df$benchmark, main="",xlab="Bio Initial Prod and Mort only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_annual_productivity_mortality') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_annual_productivity_mortality, model_benchmark_df$benchmark, main="",xlab="Bio Annual Prod and Mort only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
      }
      par(mfrow = c(5,5),oma = c(0, 0, 2, 0))
      plot_esa_cci_biomass_only
      plot_esa_cci_agb_nbe
      plot_biomass_2010
      plot_biomass_initial
      plot_biomass_annual
      plot_productivity_only
      plot_mortality_only
      plot_biomass_2010_mortality
      plot_biomass_2010_productivity
      plot_biomass_initial_mortality
      plot_biomass_initial_productivity
      plot_biomass_annual_mortality
      plot_biomass_annual_productivity
      plot_biomass_annual_productivity_nbe
      plot_productivity_mortality_only
      # plot_biomass_2010_productivity_mortality
      plot_biomass_initial_productivity_mortality
      # plot_biomass_annual_productivity_mortality
      mtext(bquote(.(cardamom_var)~.(j)~g~ m^-2~day^-1), line=-2, side=3, outer=TRUE, cex=1.5)
    }
  }
}

# extract_npp_mort_subset_models(amazonia_subset,compare_var_names[2],mod_var[1],woodprod_00_09_gCm2d,time_period)
# collate_all_models_npp_mort(amazonia_subset,compare_var_names[2],mod_var[1],reference_data[[3]],time_period[1])
# res_df_merge_plot_npp_mort(amazonia_subset,compare_var_names[3],mod_var,reference_data[[3]],time_period[2])
# res_df_merge_plot_npp_mort(amazonia_subset,compare_var_names[4],mod_var,reference_data[[6]],time_period[1])

#done
####################################################
################Residence Time######################
####################################################
res_df_merge_plot_rt <- function (region,cardamom_var,model_variant,reference,tp){
  for (j in tp) {
    if (j=='2000-2009') {
      benchmark <- extract_subset(region,reference)
      df_benchmark <- as.data.frame(benchmark, xy=TRUE)
      for (i in model_variant) {
        if (i=='esa_cci_biomass_only') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$esa_cci_biomass_only, model_benchmark_df$benchmark, main="",xlab="ESA CCI Biomass only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='esa_cci_agb_nbe') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$esa_cci_agb_nbe, model_benchmark_df$benchmark, main="",xlab="ESA CCI Biomass and NBE", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_2010') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_2010, model_benchmark_df$benchmark, main="",xlab="2010 Biomass only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_initial') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_initial, model_benchmark_df$benchmark, main="",xlab="Initial Biomass Only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_annual') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_annual, model_benchmark_df$benchmark, main="",xlab="Annual Biomass Only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='productivity_only') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df[,4], model_benchmark_df[,3], main="",xlab="Productivity only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='mortality_only') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df[,4], model_benchmark_df[,3], main="",xlab="Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_2010_mortality') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_2010_mortality, model_benchmark_df$benchmark, main="",xlab="2010 Biomass and Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_2010_productivity') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_2010_productivity, model_benchmark_df$benchmark, main="",xlab="2010 Biomass and Productivity only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_initial_mortality') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_initial_mortality, model_benchmark_df$benchmark, main="",xlab="Biomass Initial and Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_initial_productivity') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_initial_productivity, model_benchmark_df$benchmark, main="",xlab="Biomass Initial and Productivity only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_annual_mortality') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_annual_mortality, model_benchmark_df$benchmark, main="",xlab="Biomass annual and Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_annual_productivity') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_annual_productivity, model_benchmark_df$benchmark, main="",xlab="Biomass annual and Productivity only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_annual_productivity_nbe') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_annual_productivity_nbe, model_benchmark_df$benchmark, main="",xlab="Biomass annual and Productivity only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='productivity_mortality_only') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df[,4], model_benchmark_df[,3], main="",xlab="Productivity and Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_2010_productivity_mortality') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_2010_productivity_mortality, model_benchmark_df$benchmark, main="",xlab="2010 Bio Prod and Mort only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_initial_productivity_mortality') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_initial_productivity_mortality, model_benchmark_df$benchmark, main="",xlab="Bio Initial Prod and Mort only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_annual_productivity_mortality') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_annual_productivity_mortality, model_benchmark_df$benchmark, main="",xlab="Bio Annual Prod and Mort only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
      }
      par(mfrow = c(5,5),oma = c(0, 0, 2, 0))
      plot_esa_cci_biomass_only
      plot_esa_cci_agb_nbe
      plot_biomass_2010
      plot_biomass_initial
      plot_biomass_annual
      plot_productivity_only
      plot_mortality_only
      plot_biomass_2010_mortality
      plot_biomass_2010_productivity
      plot_biomass_initial_mortality
      plot_biomass_initial_productivity
      plot_biomass_annual_mortality
      plot_biomass_annual_productivity
      plot_biomass_annual_productivity_nbe
      plot_productivity_mortality_only
      # plot_biomass_2010_productivity_mortality
      plot_biomass_initial_productivity_mortality
      # plot_biomass_annual_productivity_mortality
      mtext(bquote(.(cardamom_var)~.(j)~years), line=-2, side=3, outer=TRUE, cex=1.5)
    }
    else if (j=='2010-2016') {
      benchmark <- extract_subset(region,reference)
      df_benchmark <- as.data.frame(benchmark, xy=TRUE)
      for (i in model_variant) {
        if (i=='esa_cci_biomass_only') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$esa_cci_biomass_only, model_benchmark_df$benchmark, main="",xlab="ESA CCI Biomass only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='esa_cci_agb_nbe') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$esa_cci_agb_nbe, model_benchmark_df$benchmark, main="",xlab="ESA CCI Biomass and NBE", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_2010') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_2010, model_benchmark_df$benchmark, main="",xlab="2010 Biomass Only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_initial') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_initial, model_benchmark_df$benchmark, main="",xlab="Initial Biomass Only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_annual') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_annual, model_benchmark_df$benchmark, main="",xlab="Annual Biomass Only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='productivity_only') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df[,4], model_benchmark_df[,3], main="",xlab="Productivity only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='mortality_only') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df[,4], model_benchmark_df[,3], main="",xlab="Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_2010_mortality') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_2010_mortality, model_benchmark_df$benchmark, main="",xlab="2010 Biomass and Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_2010_productivity') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_2010_productivity, model_benchmark_df$benchmark, main="",xlab="2010 Biomass and Productivity only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_initial_mortality') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_initial_mortality, model_benchmark_df$benchmark, main="",xlab="Biomass Initial and Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_initial_productivity') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_initial_productivity, model_benchmark_df$benchmark, main="",xlab="Biomass Initial and Productivity only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_annual_mortality') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_annual_mortality, model_benchmark_df$benchmark, main="",xlab="Biomass annual and Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_annual_productivity') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_annual_productivity, model_benchmark_df$benchmark, main="",xlab="Biomass annual and Productivity only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_annual_productivity_nbe') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_annual_productivity_nbe, model_benchmark_df$benchmark, main="",xlab="Biomass annual and Productivity only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='productivity_mortality_only') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df[,4], model_benchmark_df[,3], main="",xlab="Productivity and Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_2010_productivity_mortality') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_2010_productivity_mortality, model_benchmark_df$benchmark, main="",xlab="2010 Bio Prod and Mort only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_initial_productivity_mortality') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_initial_productivity_mortality, model_benchmark_df$benchmark, main="",xlab="Bio Initial Prod and Mort only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
        else if (i=='biomass_annual_productivity_mortality') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]), digits = 2)
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_annual_productivity_mortality, model_benchmark_df$benchmark, main="",xlab="Bio Annual Prod and Mort only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,80), ylim=c(0,80), 
                      cex.lab=1.5, cex.axis=1.5)+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(60,10, paste('RMSE = ',model_rmse, sep=''),cex=1.5))
        }
      }
      par(mfrow = c(5,5),oma = c(0, 0, 2, 0))
      plot_esa_cci_biomass_only
      plot_esa_cci_agb_nbe
      plot_biomass_2010
      plot_biomass_initial
      plot_biomass_annual
      plot_productivity_only
      plot_mortality_only
      plot_biomass_2010_mortality
      plot_biomass_2010_productivity
      plot_biomass_initial_mortality
      plot_biomass_initial_productivity
      plot_biomass_annual_mortality
      plot_biomass_annual_productivity
      plot_biomass_annual_productivity_nbe
      plot_productivity_mortality_only
      # plot_biomass_2010_productivity_mortality
      plot_biomass_initial_productivity_mortality
      # plot_biomass_annual_productivity_mortality
      mtext(bquote(.(cardamom_var)~.(j)~years), line=-2, side=3, outer=TRUE, cex=1.5)
    }
  }
}
# res_df_merge_plot_rt(amazonia_subset,compare_var_names[4],mod_var,reference_data[[6]],time_period[2])

#done
#######################################################
################Combine all variables##################
#######################################################
one_function_to_plot_them_all <- function(region,cardamom_var,model_variant,reference,tp) {
  for (i in cardamom_var) {
    if (cardamom_var=='WOOD') {
      for (j in tp) {
        if (tp =='2000-2009'){res_df_merge_plot_cwood(region,i,model_variant,reference[[1]],j)
        }
        else if (tp =='2010-2016'){res_df_merge_plot_cwood(region,i,model_variant,reference[[1]],j)
        }
      }
    }
    else if (cardamom_var=='NPP_wood_flx') {
      for (j in tp) {
        if (tp =='2000-2009'){res_df_merge_plot_npp_mort(region,i,model_variant,reference[[2]],j)
        }
        else if (tp =='2010-2016'){res_df_merge_plot_npp_mort(region,i,model_variant,reference[[3]],j)
        }
      }
      }
    else if (cardamom_var=='OUTPUT_wood_flx') {
      for (j in tp) {
        if (tp =='2000-2009'){res_df_merge_plot_npp_mort(region,i,model_variant,reference[[4]],j)
        }
        else if (tp =='2010-2016'){res_df_merge_plot_npp_mort(region,i,model_variant,reference[[5]],j)
        }
      }
    }
    else if (cardamom_var=='MTT_wood') {
      for (j in tp) {
        if (tp =='2000-2009'){res_df_merge_plot_rt(region,i,model_variant,reference[[6]],j)
        }
        else if (tp =='2010-2016'){res_df_merge_plot_rt(region,i,model_variant,reference[[7]],j)
        }
      }
    }
  }
}

# par(mfrow = c(3,6));one_function_to_plot_them_all(amazonia_subset,compare_var_names[1],mod_var,reference_data,time_period[1])
# par(mfrow = c(3,6));one_function_to_plot_them_all(amazonia_subset,compare_var_names[1],mod_var,reference_data,time_period[2])
# par(mfrow = c(3,6));one_function_to_plot_them_all(amazonia_subset,compare_var_names[2],mod_var,reference_data,time_period[1])
# par(mfrow = c(3,6));one_function_to_plot_them_all(amazonia_subset,compare_var_names[2],mod_var,reference_data,time_period[2])
# par(mfrow = c(3,6));one_function_to_plot_them_all(amazonia_subset,compare_var_names[3],mod_var,reference_data,time_period[1])
# par(mfrow = c(3,6));one_function_to_plot_them_all(amazonia_subset,compare_var_names[3],mod_var,reference_data,time_period[2])
# par(mfrow = c(3,6));one_function_to_plot_them_all(amazonia_subset,compare_var_names[4],mod_var,reference_data,time_period[1])
# par(mfrow = c(3,6));one_function_to_plot_them_all(amazonia_subset,compare_var_names[4],mod_var,reference_data,time_period[2])
#done

#######################################################
###############extract model parameters################
#######################################################
# prefix_new <- 'M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/'
# midfix <- '_subset/RESULTS_PROCESSED/'
# suffix_new <- '_subset_stock_flux.RData'
# 
# mod_var_new <- c('esa_cci_agb_nbe','esa_cci_agb_only',
#                  'Rainfor_woody_biomass_2010','Rainfor_woody_biomass_initial','Rainfor_woody_biomass_annual',
#                  'Rainfor_woody_mortality_only','Rainfor_woody_productivity_only',
#                  'Rainfor_woody_biomass_2010_productivity','Rainfor_woody_biomass_2010_and_mortality',
#                  'Rainfor_woody_biomass_initial_productivity','Rainfor_woody_biomass_annual_productivity',
#                  'Rainfor_woody_biomass_initial_mortality','Rainfor_woody_biomass_annual_mortality',
#                  'Rainfor_woody_productivity_and_mortality_only','Rainfor_woody_biomass_2010_productivity_mortality',
#                  'Rainfor_woody_biomass_initial_productivity_mortality','Rainfor_woody_biomass_annual_productivity_mortality'
# )
# mod_var_abr <- c('NWD','ECB','B_2010','B_i','B_a','M','P','B_2010_P','B_2010_M',
#                  'B_i_P','B_a_P','B_i_M','B_a_M','P_M','B_2010_P_M','B_i_P_M','B_a_P_M')

model_pars_list <- list()
# model_pars_list_sd <- list()  

extract_parameters_all_models <- function (model_variant,model_variant_abr){
  for (i in seq_along(model_variant)) {
    load(paste(prefix_new,model_variant[i],midfix,model_variant[i],suffix_new,sep=""))
    data <- grid_output$parameters[,,,4]
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

# model_pars_df<-as.data.frame(t(as.data.frame(extract_parameters_all_models(mod_var_new,mod_var_abr))), stringsAsFactors = FALSE)
# model_pars_df[] <- lapply(model_pars_df, type.convert, as.is = TRUE)
# rownames(model_pars_df)<-  seq_along(model_pars_df[,1])
# colnames(model_pars_df)<-  par_names
# model_pars_df$model_var<-as.factor(model_pars_df$model_var)
# model_pars_df<-model_pars_df[,c(34,1:32)]
# new_model_pars_df<-round_df(model_pars_df[,2:33], 3)
# rownames(new_model_pars_df)<-mod_var_new
# print(new_model_pars_df)
# plot_all_parameters(model_pars_df,param_names)
