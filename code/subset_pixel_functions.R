#All codes being Used
#import libraries
#install packages----
library(ncdf4); library(raster); library(dplyr); library(ggplot2);library(ggpubr);library(quantreg);library(ggpp);library(rgeos);library(ggpmisc);library(rgdal);library(Metrics)

#done
#data prep pre-cursors----
compare_var_names <- c('WOOD','NPP_wood_flx','OUTPUT_wood_flx')
amazonia_subset <- shapefile("R:/cssp_brazil/cssp_brazil_R/data/amazonia_subset.shp")
mod_var <- c('norainfor','biomass','productivity','mortality','productivity_mortality','biomass_productivity_mortality')
time_period <- c('01/09','10/16')
prefix <- 'R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/DATA/benchmark/Amazon_subset_'
suffix <- '_1deg_monthly_2001_updated_2019.nc'
reference_data<- c(biomass_amazon_gCm2,woodprod_00_09_gCm2d,woodprod_10_16_gCm2d,biommort_00_09_gCm2d,biommort_10_16_gCm2d)

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
extract_biomass_subset_all_models <- function (region,cardamom_var,model_variant,reference){
  for (i in model_variant) {
    if (i=='norainfor') {
      data <- brick(paste(prefix,i,suffix,sep=""),varname=cardamom_var)
      data <- crop(data, extent(reference))
      a <- data$X2001.01.01
      cropped <- crop(a, extent(region))
      masked1 <- mask(cropped, region)
      masked2 <- masked1 > -Inf
      poly <- rasterToPolygons(masked2, dissolve=TRUE)
      data_region <- mask(data,poly)
      data_2010<-data_region[[which(getZ(data_region) >= as.Date("2010-01-01") & getZ(data_region) <= as.Date("2010-12-01"))]]
      data_2010_mean <- stackApply(data_2010, indices =  rep(1,nlayers(data_2010)), fun = "mean")
    }
    else {
      data <- brick(paste(prefix,i,suffix,sep=""),varname=cardamom_var)
      data_2010<-data[[which(getZ(data) >= as.Date("2010-01-01") & getZ(data) <= as.Date("2010-12-01"))]]
      data_2010_mean <- stackApply(data_2010, indices =  rep(1,nlayers(data_2010)), fun = "mean")
    }
  }
  return(data_2010_mean)
}
# extract_biomass_subset_all_models(amazonia_subset,compare_var_names[1],mod_var[1],biomass_amazon_gCm2)

collate_all_models_cwood <- function (region,cardamom_var,model_variant,reference){
  for (i in model_variant) {
    # if (i=='norainfor') {
    #   cwood_raster <- extract_biomass_subset_all_models(region,cardamom_var,i,reference)
    # }
    # else if (i=='biomass') {
    #   cwood_raster <- extract_biomass_subset_all_models(region,cardamom_var,i,reference)
    # }
    # else if (i=='productivity') {
    #   cwood_raster <- extract_biomass_subset_all_models(region,cardamom_var,i,reference)
    # }
    # else if (i=='mortality') {
    #   cwood_raster <- extract_biomass_subset_all_models(region,cardamom_var,i,reference)
    # }
    # else if (i=='biomass_productivity_mortality') {
    #   cwood_raster <- extract_biomass_subset_all_models(region,cardamom_var,i,reference)
    # }
    cwood_raster <- extract_biomass_subset_all_models(region,cardamom_var,i,reference)
  }
  # cwood_list <-c(cwood_norainfor,cwood_biomass,cwood_productivity,cwood_mortality)
  # names(cwood_list) <- model_variant
  return(cwood_raster)
}
#collate_all_models_cwood(amazonia_subset,compare_var_names[1],mod_var[1],biomass_amazon_gCm2)

res_df_merge_plot_cwood <- function (region,cardamom_var,model_variant,reference){
  benchmark <- extract_subset(region,reference)
  df_benchmark <- as.data.frame(benchmark, xy=TRUE)
  for (i in model_variant) {
    if (i=='norainfor') {
      model <- collate_all_models_cwood(region,cardamom_var,i,reference)
      model_res <- resample(model,benchmark)
      model_res_df <- as.data.frame(model_res, xy=TRUE)
      model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
      names(model_benchmark_df)<-c('x','y','benchmark',i)
      model_benchmark_df<-na.omit(model_benchmark_df)
      model_rmse<- rmse(model_benchmark_df[,3],model_benchmark_df[,4])
      assign(paste("plot", i, sep = "_"),
             plot(model_benchmark_df$norainfor, model_benchmark_df$benchmark, main="",xlab="NORAINFOR ESA CCI Biomass", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000))+
               abline(coef = c(0,1),col='red', lwd=3)+
               text(20000, 2000, paste('RMSE = ',model_rmse, sep='')))
    }
    else if (i=='biomass') {
      model <- collate_all_models_cwood(region,cardamom_var,i,reference)
      model_res <- resample(model,benchmark)
      model_res_df <- as.data.frame(model_res, xy=TRUE)
      model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
      names(model_benchmark_df)<-c('x','y','benchmark',i)
      model_benchmark_df<-na.omit(model_benchmark_df)
      model_rmse<- rmse(model_benchmark_df[,3],model_benchmark_df[,4])
      assign(paste("plot", i, sep = "_"),
             plot(model_benchmark_df$biomass, model_benchmark_df$benchmark, main="",xlab="Biomass Only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000))+
               abline(coef = c(0,1),col='red', lwd=3)+
               text(20000, 2000, paste('RMSE = ',model_rmse, sep='')))
    }
    else if (i=='productivity') {
      model <- collate_all_models_cwood(region,cardamom_var,i,reference)
      model_res <- resample(model,benchmark)
      model_res_df <- as.data.frame(model_res, xy=TRUE)
      model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
      names(model_benchmark_df)<-c('x','y','benchmark',i)
      model_benchmark_df<-na.omit(model_benchmark_df)
      model_rmse<- rmse(model_benchmark_df[,3],model_benchmark_df[,4])
      assign(paste("plot", i, sep = "_"),
             plot(model_benchmark_df$productivity, model_benchmark_df$benchmark, main="",xlab="Productivity only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000))+
               abline(coef = c(0,1),col='red', lwd=3)+
               text(20000, 2000, paste('RMSE = ',model_rmse, sep='')))
    }
    else if (i=='mortality') {
      model <- collate_all_models_cwood(region,cardamom_var,i,reference)
      model_res <- resample(model,benchmark)
      model_res_df <- as.data.frame(model_res, xy=TRUE)
      model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
      names(model_benchmark_df)<-c('x','y','benchmark',i)
      model_benchmark_df<-na.omit(model_benchmark_df)
      model_rmse<- rmse(model_benchmark_df[,3],model_benchmark_df[,4])
      assign(paste("plot", i, sep = "_"),
             plot(model_benchmark_df$mortality, model_benchmark_df$benchmark, main="",xlab="Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000))+
             abline(coef = c(0,1),col='red', lwd=3)+
               text(20000, 2000, paste('RMSE = ',model_rmse, sep='')))
    }
    else if (i=='productivity_mortality') {
      model <- collate_all_models_cwood(region,cardamom_var,i,reference)
      model_res <- resample(model,benchmark)
      model_res_df <- as.data.frame(model_res, xy=TRUE)
      model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
      names(model_benchmark_df)<-c('x','y','benchmark',i)
      model_benchmark_df<-na.omit(model_benchmark_df)
      model_rmse<- rmse(model_benchmark_df[,3],model_benchmark_df[,4])
      assign(paste("plot", i, sep = "_"),
             plot(model_benchmark_df$productivity, model_benchmark_df$benchmark, main="",xlab="Productivity and Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000))+
               abline(coef = c(0,1),col='red', lwd=3)+
               text(20000, 2000, paste('RMSE = ',model_rmse, sep='')))
    }
    else if (i=='biomass_productivity_mortality') {
      model <- collate_all_models_cwood(region,cardamom_var,i,reference)
      model_res <- resample(model,benchmark)
      model_res_df <- as.data.frame(model_res, xy=TRUE)
      model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
      names(model_benchmark_df)<-c('x','y','benchmark',i)
      model_benchmark_df<-na.omit(model_benchmark_df)
      model_rmse<- rmse(model_benchmark_df[,3],model_benchmark_df[,4])
      assign(paste("plot", i, sep = "_"),
             plot(model_benchmark_df$biomass_productivity_mortality, model_benchmark_df$benchmark, main="",xlab="Bio Prod and Mort only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,25000), ylim=c(0,25000))+
               abline(coef = c(0,1),col='red', lwd=3)+
               text(20000, 2000, paste('RMSE = ',model_rmse, sep='')))
    }
  }
  par(mfrow = c(2,3))
  plot_norainfor
  plot_biomass
  plot_productivity
  plot_mortality
  plot_productivity_mortality
  plot_biomass_productivity_mortality
  mtext("Aboveground Biomass 2010", line=-2, side=3, outer=TRUE, cex=1.5)
  }
#par(mfrow = c(2,3));res_df_merge_plot_cwood(amazonia_subset,compare_var_names[1],mod_var,biomass_amazon_gCm2)

#done
#######################################################
################NPP and mortality######################
#######################################################
#extract subset of wood, npp and mortality from models
extract_npp_mort_subset_models <- function (region,cardamom_var,model_variant,reference,tp){
  for (i in model_variant) {
    if (i=='norainfor') {
      data <- brick(paste(prefix,i,suffix,sep=""),varname=cardamom_var)
      data <- crop(data, extent(reference))
      a <- data$X2001.01.01
      cropped <- crop(a, extent(region))
      masked1 <- mask(cropped, region)
      masked2 <- masked1 > -Inf
      poly <- rasterToPolygons(masked2, dissolve=TRUE)
      data_region <- mask(data,poly)
      data_0109<-data_region[[which(getZ(data_region) >= as.Date("2001-01-01") & getZ(data_region) <= as.Date("2009-12-01"))]]
      data_1016<-data_region[[which(getZ(data_region) >= as.Date("2010-01-01") & getZ(data_region) <= as.Date("2016-12-01"))]]
      data_0109_mean <- stackApply(data_0109, indices =  rep(1,nlayers(data_0109)), fun = "mean")
      data_1016_mean <- stackApply(data_1016, indices =  rep(1,nlayers(data_1016)), fun = "mean")
    }
    else {
      data <- brick(paste(prefix,model_variant[1],suffix,sep=""),varname=cardamom_var)
      data_0109<-data[[which(getZ(data) >= as.Date("2001-01-01") & getZ(data) <= as.Date("2009-12-01"))]]
      data_1016<-data[[which(getZ(data) >= as.Date("2010-01-01") & getZ(data) <= as.Date("2016-12-01"))]]
      data_0109_mean <- stackApply(data_0109, indices =  rep(1,nlayers(data_0109)), fun = "mean")
      data_1016_mean <- stackApply(data_1016, indices =  rep(1,nlayers(data_1016)), fun = "mean")
    }
  }
  return(c(data_0109_mean,data_1016_mean))
}

# extract_npp_mort_subset_models(amazonia_subset,compare_var_names[2],mod_var[1],woodprod_00_09_gCm2d,time_period)
collate_all_models_npp_mort <- function (region,cardamom_var,model_variant,reference,tp){
  for (i in model_variant) {
    # if (i=='norainfor') {
    #   wood_npp_mort_raster <- extract_npp_mort_subset_models(region,cardamom_var,i,reference,j)
    # }
    # else if (i=='biomass') {
    #   wood_npp_mort_raster <- extract_npp_mort_subset_models(region,cardamom_var,i,reference,j)
    # }
    # else if (i=='productivity') {
    #   wood_npp_mort_raster <- extract_npp_mort_subset_models(region,cardamom_var,i,reference,j)
    # }
    # else if (i=='mortality') {
    #   wood_npp_mort_raster <- extract_npp_mort_subset_models(region,cardamom_var,i,reference,j)
    # }
    # else if (i=='biomass_productivity_mortality') {
    #   wood_npp_mort_raster <- extract_npp_mort_subset_models(region,cardamom_var,i,reference,j)
    # }
    wood_npp_mort_raster <- extract_npp_mort_subset_models(region,cardamom_var,i,reference,j)
  }
  for (j in tp) {
    if (j=='01/09'){
      wood_npp_mort_raster <-wood_npp_mort_raster[[1]]
      return(wood_npp_mort_raster)
    }
    else if (j=='10/16') {
      wood_npp_mort_raster <-wood_npp_mort_raster[[2]]
      return(wood_npp_mort_raster)
    }
  }
}
#collate_all_models_npp_mort(amazonia_subset,compare_var_names[2],mod_var[1],reference_data[[3]],time_period[1])

#extract and plot for each variable
res_df_merge_plot_npp_mort <- function (region,cardamom_var,model_variant,reference,tp){
  for (j in tp) {
    if (j=='01/09') {
      benchmark <- extract_subset(region,reference)
      df_benchmark <- as.data.frame(benchmark, xy=TRUE)
      for (i in model_variant) {
        if (i=='norainfor') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- rmse(model_benchmark_df[,3],model_benchmark_df[,4])
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$norainfor, model_benchmark_df$benchmark, main="",xlab="NORAINFOR ESA CCI Biomass", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4))+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep='')))
        }
        else if (i=='biomass') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- rmse(model_benchmark_df[,3],model_benchmark_df[,4])
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass, model_benchmark_df$benchmark, main="",xlab="Biomass Only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4))+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep='')))
        }
        else if (i=='productivity') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- rmse(model_benchmark_df[,3],model_benchmark_df[,4])
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$productivity, model_benchmark_df$benchmark, main="",xlab="Productivity only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4))+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep='')))
        }
        else if (i=='mortality') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- rmse(model_benchmark_df[,3],model_benchmark_df[,4])
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$mortality, model_benchmark_df$benchmark, main="",xlab="Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4))+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep='')))
        }
        else if (i=='productivity_mortality') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- rmse(model_benchmark_df[,3],model_benchmark_df[,4])
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$productivity, model_benchmark_df$benchmark, main="",xlab="Productivity and Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4))+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep='')))
        }
        else if (i=='biomass_productivity_mortality') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- rmse(model_benchmark_df[,3],model_benchmark_df[,4])
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_productivity_mortality, model_benchmark_df$benchmark, main="",xlab="Bio Prod and Mort only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4))+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep='')))
        }
      }
      par(mfrow = c(2,3))
      plot_norainfor
      plot_biomass
      plot_productivity
      plot_mortality
      plot_productivity_mortality
      plot_biomass_productivity_mortality
      mtext(paste(cardamom_var, j, sep = "_"), line=-2, side=3, outer=TRUE, cex=1.5)
    }
    else if (j=='10/16') {
      benchmark <- extract_subset(region,reference)
      df_benchmark <- as.data.frame(benchmark, xy=TRUE)
      for (i in model_variant) {
        if (i=='norainfor') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- rmse(model_benchmark_df[,3],model_benchmark_df[,4])
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$norainfor, model_benchmark_df$benchmark, main="",xlab="NORAINFOR ESA CCI Biomass", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4))+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep='')))
        }
        else if (i=='biomass') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- rmse(model_benchmark_df[,3],model_benchmark_df[,4])
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass, model_benchmark_df$benchmark, main="",xlab="Biomass Only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4))+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep='')))
        }
        else if (i=='productivity') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- rmse(model_benchmark_df[,3],model_benchmark_df[,4])
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$productivity, model_benchmark_df$benchmark, main="",xlab="Productivity only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4))+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep='')))
        }
        else if (i=='mortality') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- rmse(model_benchmark_df[,3],model_benchmark_df[,4])
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$mortality, model_benchmark_df$benchmark, main="",xlab="Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4))+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep='')))
        }
        else if (i=='productivity_mortality') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- rmse(model_benchmark_df[,3],model_benchmark_df[,4])
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$productivity, model_benchmark_df$benchmark, main="",xlab="Productivity and Mortality only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4))+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep='')))
        }
        else if (i=='biomass_productivity_mortality') {
          model <- collate_all_models_npp_mort(region,cardamom_var,i,reference,j)
          model_res <- resample(model,benchmark)
          model_res_df <- as.data.frame(model_res, xy=TRUE)
          model_benchmark_df <- merge(df_benchmark, model_res_df, by=c("x","y"))
          names(model_benchmark_df)<-c('x','y','benchmark',i)
          model_benchmark_df<-na.omit(model_benchmark_df)
          model_rmse<- rmse(model_benchmark_df[,3],model_benchmark_df[,4])
          assign(paste("plot", i, sep = "_"),
                 plot(model_benchmark_df$biomass_productivity_mortality, model_benchmark_df$benchmark, main="",xlab="Bio Prod and Mort only", ylab="RAINFOR benchmark", pch=19, xlim=c(0,4), ylim=c(0,4))+
                   abline(coef = c(0,1),col='red', lwd=3)+
                   text(1, 3.5, paste('RMSE = ',model_rmse, sep='')))
        }
      }
      par(mfrow = c(2,3))
      plot_norainfor
      plot_biomass
      plot_productivity
      plot_mortality
      plot_productivity_mortality
      plot_biomass_productivity_mortality
      mtext(paste(cardamom_var, j, sep = "_"), line=-2, side=3, outer=TRUE, cex=1.5)
    }
  }
}
#res_df_merge_plot_npp_mort(amazonia_subset,compare_var_names[3],mod_var,reference_data[[3]],time_period[2])

#done

#######################################################
################Combine all variables##################
#######################################################
one_function_to_plot_them_all <- function(region,cardamom_var,model_variant,reference,tp) {
  for (i in cardamom_var) {
    if (cardamom_var=='WOOD') {
      res_df_merge_plot_cwood(region,i,model_variant,reference[[1]])
    }
    else if (cardamom_var=='NPP_wood_flx') {
      for (j in tp) {
        if (tp =='01/09'){res_df_merge_plot_npp_mort(region,i,model_variant,reference[[2]],j)
        }
        else if (tp =='10/16'){res_df_merge_plot_npp_mort(region,i,model_variant,reference[[3]],j)
        }
      }
      }
    else if (cardamom_var=='OUTPUT_wood_flx') {
      for (j in tp) {
        if (tp =='01/09'){res_df_merge_plot_npp_mort(region,i,model_variant,reference[[4]],j)
        }
        else if (tp =='10/16'){res_df_merge_plot_npp_mort(region,i,model_variant,reference[[5]],j)
        }
      }
    }
  }
}
par(mfrow = c(2,3))
one_function_to_plot_them_all(amazonia_subset,compare_var_names[1],mod_var,reference_data,time_period)
one_function_to_plot_them_all(amazonia_subset,compare_var_names[2],mod_var,reference_data,time_period[1])
one_function_to_plot_them_all(amazonia_subset,compare_var_names[2],mod_var,reference_data,time_period[2])
one_function_to_plot_them_all(amazonia_subset,compare_var_names[3],mod_var,reference_data,time_period[1])
one_function_to_plot_them_all(amazonia_subset,compare_var_names[3],mod_var,reference_data,time_period[2])
#done