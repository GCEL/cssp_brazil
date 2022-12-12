############Pre processing###############
library(ncdf4); library(raster); library(dplyr);library(ggplot2);library(ggpubr);library(quantreg);library(ggpp);library(rgeos);library(ggpmisc);library(rgdal);library(tidyverse);library(lubridate)

#end
# for (f in 1:length(nbe_list)) {
#   nbe_list [[f]] <- flip(nbe_list[[f]],2)
# }
# for (f in 1:length(nbe_unc_list)) {
#   nbe_unc_list [[f]] <- flip(nbe_unc_list[[f]],2)
# }

nbe_stack<- stack(nbe_list)
nbe_unc_stack<- stack(nbe_unc_list)


layer_dates_2001_2019 <- getZ(cardamom_nppwood)
layer_dates_2001_2019 <- paste("X",layer_dates_2001_2019,sep="")
layer_dates_2015_2019 <- layer_dates_2001_2019[169:228]
names(lai_b) <- layer_dates_2001_2019;names(lai_unc_b) <- layer_dates_2001_2019

layer_date_format_2001_2019 <- getZ(cardamom_nppwood)
layer_date_format_2015_2019 <- layer_date_format_2001_2019[169:228]

names(nbe_stack) <- layer_dates_2015_2019;names(nbe_unc_stack) <- layer_dates_2015_2019;names(lai_b)<-layer_dates_2015_2019
nbe_stack_15_19_monthly <- setZ(nbe_stack,layer_date_format_2015_2019,"Date");nbe_unc_stack_15_19_monthly <- setZ(nbe_unc_stack,layer_date_format_2015_2019,"Date")
lai_stack_01_19_monthly <- setZ(lai_b,layer_date_format_2001_2019,"Date");lai_unc_stack_01_19_monthly <- setZ(lai_unc_b,layer_date_format_2001_2019,"Date")

nbe_prefix <- 'R://ILAMB_beta_devel/RAINFOR_leeds_run/nbe_analysis/'
nbe_suffix <- '_Amazon_1deg_monthly_2001_updated_2019.nc'
nbe_model_variants <- c('raw_GEOSCHEM','default_CARDAMOM','geoschem_CARDAMOM')
lai_model_variants <- c('raw_COPERNICUS','default_CARDAMOM','geoschem_CARDAMOM')

nbe_time_period <- c('2015-2019')
amazonia_nbe_subset <- shapefile("R:/cssp_brazil/cssp_brazil_R/data/amazon_site_nbe_pixels.shp")
reference_nbe_data<- brick('R://brazil_leeds_maps/AbovegroundBiomass_Mg_perHa_111km.tif')

amazon_nw_poly <- shapefile("R:/cssp_brazil/cssp_brazil_R/data/amazon_nw.shp")
amazon_sw_poly <- shapefile("R:/cssp_brazil/cssp_brazil_R/data/amazon_sw.shp")
amazon_ec_poly <- shapefile("R:/cssp_brazil/cssp_brazil_R/data/amazon_ec.shp")
amazon_bs_poly <- shapefile("R:/cssp_brazil/cssp_brazil_R/data/amazon_bs.shp")
amazon_gs_poly <- shapefile("R:/cssp_brazil/cssp_brazil_R/data/amazon_gs.shp")

amazon_nw_pixel <- shapefile("R:/brazil_leeds_maps/amazon_nw.shp")
amazon_sw_pixel <- shapefile("R:/brazil_leeds_maps/amazon_sw.shp")
amazon_ec_pixel <- shapefile("R:/brazil_leeds_maps/amazon_ec.shp")
amazon_bs_pixel <- shapefile("R:/brazil_leeds_maps/amazon_bs.shp")
amazon_gs_pixel <- shapefile("R:/brazil_leeds_maps/amazon_gs.shp")

nbe_var <- c('NBE','NBE_2pt5pc','NBE_97pt5pc')
lai_var <- c('LAI','LAI_2pt5pc','LAI_97pt5pc')

pre_nbe_reg_data <- list(amazon_nw_poly,amazon_sw_poly,amazon_ec_poly,amazon_bs_poly,amazon_gs_poly)
pre_nbe_pixel_data <- list(amazon_nw_pixel,amazon_sw_pixel,amazon_ec_pixel,amazon_bs_pixel,amazon_gs_pixel)

nbe_reg_names <- c('amazon_nw','amazon_sw','amazon_ec','amazon_bs','amazon_gs','amazonia')
nbe_pixel_names <- c('amazon_nw','amazon_sw','amazon_ec','amazon_bs','amazon_gs')

extract_amazonia_regions <- function (region,reference){
  masked1 <- mask(reference, region)
  masked2 <- masked1 > -Inf
  poly <- rasterToPolygons(masked2, dissolve=TRUE)
  return(poly)
}
amazonia_nbe_subset_poly <- extract_amazonia_regions(amazonia_nbe_subset,reference_nbe_data)

nbe_reg_data <- list()
nbe_pixel_data <- list()

for (r in 1:length(pre_nbe_reg_data)){
  nbe_reg_data[r]<-extract_amazonia_regions(pre_nbe_reg_data[[r]],reference_nbe_data)
  }

for (r in 1:length(pre_nbe_pixel_data)){
  nbe_pixel_data[r]<-extract_amazonia_regions(pre_nbe_pixel_data[[r]],reference_nbe_data)
}

reference_nbe_data_mask <- reference_nbe_data > -Inf
nbe_reg_data[6] <- rasterToPolygons(reference_nbe_data_mask, dissolve=TRUE)

##########NBE##############
####extract nbe data frame
extract_nbe_all_models <- function (model_variant,region,region_name,variable_name,reference,tp){
  for (i in model_variant) {
    if (i=='raw_GEOSCHEM'){
      data <- nbe_stack_15_19_monthly
      } else {
      data <- brick(paste(nbe_prefix,model_variant[1],nbe_suffix,sep=""),varname=variable_name)
      }
    a <- data$X2015.01.01
    cropped <- crop(a, extent(region))
    masked1 <- mask(cropped, region)
    masked2 <- masked1 > -Inf
    poly <- rasterToPolygons(masked2, dissolve=TRUE)
    data_region <- mask(data,poly)
    data_region_1519<-data_region[[which(getZ(data_region) >= as.Date("2015-01-01") & getZ(data_region) <= as.Date("2019-12-01"))]]
    df_data_region_1519 <- as.data.frame(data_region_1519,xy=F)
    df_data_region_1519<-as.data.frame(t(df_data_region_1519))
    df_data_region_1519 <- data.frame(y=unlist(df_data_region_1519))
    names(df_data_region_1519) <- paste(tolower(variable_name))
    regional_name <- as.factor(rep(region_name,nrow(df_data_region_1519)))
    model_variant_name <- as.factor(rep(model_variant,nrow(df_data_region_1519)))
    no_dp<-nrow(df_data_region_1519)/60
    date<- as.Date(rep(seq(from = as.Date("2015-01-01"), to = as.Date("2019-12-01"), by = 'month'),no_dp))
    month<-format(date,format="%b")
    year <- format(date,format="%Y")
  }
  return(na.omit(cbind(df_data_region_1519,regional_name,date,month,year,model_variant_name)))
}
# nw_nbe <- extract_nbe_all_models(nbe_model_variants[2],nbe_reg_data[[1]],nbe_reg_names[1],nbe_var[1],reference_nbe_data,nbe_time_period)

#extract nbe pixel
extract_unc_nbe_all_models <- function (model_variant,region,region_name,variable_name,reference,tp){
  for (i in model_variant) {
    if (i=='raw_GEOSCHEM'){
      data <- nbe_unc_stack_15_19_monthly
    } else {
      data <- brick(paste(nbe_prefix,model_variant[1],nbe_suffix,sep=""),varname= variable_name)
    }
    a <- data$X2015.01.01
    cropped <- crop(a, extent(region))
    masked1 <- mask(cropped, region)
    masked2 <- masked1 > -Inf
    poly <- rasterToPolygons(masked2, dissolve=TRUE)
    data_region <- mask(data,poly)
    data_region_1519<-data_region[[which(getZ(data_region) >= as.Date("2015-01-01") & getZ(data_region) <= as.Date("2019-12-01"))]]
    df_data_region_1519 <- as.data.frame(data_region_1519,xy=F)
    df_data_region_1519<-as.data.frame(t(df_data_region_1519))
    df_data_region_1519 <- data.frame(y=unlist(df_data_region_1519))
    names(df_data_region_1519) <- paste(tolower(variable_name))
    regional_name <- as.factor(rep(region_name,nrow(df_data_region_1519)))
    model_variant_name <- as.factor(rep(model_variant,nrow(df_data_region_1519)))
    no_dp<-nrow(df_data_region_1519)/60
    date<- as.Date(rep(seq(from = as.Date("2015-01-01"), to = as.Date("2019-12-01"), by = 'month'),no_dp))
    month<-format(date,format="%b")
    year <- format(date,format="%Y")
  }
  return(na.omit(cbind(df_data_region_1519,regional_name,date,month,year,model_variant_name)))
}
# nw_pixel_nbe <- extract_pixel_nbe_all_models(nbe_model_variants[1],nbe_pixel_data[[1]],nbe_pixel_names[1],nbe_var[1],reference_nbe_data,nbe_time_period)
# nw_pixel_nbe_unc <- extract_unc_nbe_all_models(nbe_model_variants[1],nbe_pixel_data[[1]],nbe_pixel_names[1],nbe_var[2],reference_nbe_data,nbe_time_period)

####extract for each region
all_nbe_reg_run <- function(model_variant,region,region_name,variable_name,reference,tp) {
  for (i in region_name) {
    if(i=='amazon_nw'){
      nw_nbe <- extract_nbe_all_models(model_variant,region[[1]],i,variable_name,reference,tp)
      row.names(nw_nbe)<-seq(1:nrow(nw_nbe))
    }
    else if (i=='amazon_sw') {
      sw_nbe <- extract_nbe_all_models(model_variant,region[[2]],i,variable_name,reference,tp)
      row.names(sw_nbe)<-seq(1:nrow(sw_nbe))
    }
    else if (i=='amazon_ec') {
      ec_nbe <- extract_nbe_all_models(model_variant,region[[3]],i,variable_name,reference,tp)
      row.names(ec_nbe)<-seq(1:nrow(ec_nbe))
    }
    else if (i=='amazon_bs') {
      bs_nbe <- extract_nbe_all_models(model_variant,region[[4]],i,variable_name,reference,tp)
      row.names(bs_nbe)<-seq(1:nrow(bs_nbe))
    }
    else if (i=='amazon_gs') {
      gs_nbe <- extract_nbe_all_models(model_variant,region[[5]],i,variable_name,reference,tp)
      row.names(gs_nbe)<-seq(1:nrow(gs_nbe))
    }
    else{
      amazon_nbe <- extract_nbe_all_models(model_variant,region[[6]],i,variable_name,reference,tp)
      row.names(amazon_nbe)<-seq(1:nrow(amazon_nbe))
    }
  }
  return(rbind(nw_nbe,sw_nbe,ec_nbe,bs_nbe,gs_nbe,amazon_nbe))
}
# all_regional_nbe <- all_nbe_reg_run(nbe_model_variants[1],nbe_reg_data,nbe_reg_names,nbe_var[1],reference_nbe_data,nbe_time_period)

all_pixel_reg_run <- function(model_variant,region,region_name,variable_name,reference,tp) {
  for (j in variable_name) {
    if (j=='NBE') {
      for (i in region_name) {
        if(i=='amazon_nw'){
          nw_nbe <- extract_nbe_all_models(model_variant,region[[1]],i,j,reference,tp)
          row.names(nw_nbe)<-seq(1:nrow(nw_nbe))
        }
        else if (i=='amazon_sw') {
          sw_nbe <- extract_nbe_all_models(model_variant,region[[2]],i,j,reference,tp)
          row.names(sw_nbe)<-seq(1:nrow(sw_nbe))
        }
        else if (i=='amazon_ec') {
          ec_nbe <- extract_nbe_all_models(model_variant,region[[3]],i,j,reference,tp)
          row.names(ec_nbe)<-seq(1:nrow(ec_nbe))
        }
        else if (i=='amazon_bs') {
          bs_nbe <- extract_nbe_all_models(model_variant,region[[4]],i,j,reference,tp)
          row.names(bs_nbe)<-seq(1:nrow(bs_nbe))
        }
        else if (i=='amazon_gs') {
          gs_nbe <- extract_nbe_all_models(model_variant,region[[5]],i,j,reference,tp)
          row.names(gs_nbe)<-seq(1:nrow(gs_nbe))
        }
      }
      nbe<- rbind(nw_nbe,sw_nbe,ec_nbe,bs_nbe,gs_nbe)
    }
    else if (j=='NBE_2pt5pc') {
      for (i in region_name) {
        if(i=='amazon_nw'){
          nw_nbe <- extract_unc_nbe_all_models(model_variant,region[[1]],i,j,reference,tp)
          row.names(nw_nbe)<-seq(1:nrow(nw_nbe))
        }
        else if (i=='amazon_sw') {
          sw_nbe <- extract_unc_nbe_all_models(model_variant,region[[2]],i,j,reference,tp)
          row.names(sw_nbe)<-seq(1:nrow(sw_nbe))
        }
        else if (i=='amazon_ec') {
          ec_nbe <- extract_unc_nbe_all_models(model_variant,region[[3]],i,j,reference,tp)
          row.names(ec_nbe)<-seq(1:nrow(ec_nbe))
        }
        else if (i=='amazon_bs') {
          bs_nbe <- extract_unc_nbe_all_models(model_variant,region[[4]],i,j,reference,tp)
          row.names(bs_nbe)<-seq(1:nrow(bs_nbe))
        }
        else if (i=='amazon_gs') {
          gs_nbe <- extract_unc_nbe_all_models(model_variant,region[[5]],i,j,reference,tp)
          row.names(gs_nbe)<-seq(1:nrow(gs_nbe))
        }
      }
      nbe_lower<- rbind(nw_nbe,sw_nbe,ec_nbe,bs_nbe,gs_nbe)
    }
    else {
      for (i in region_name) {
        if(i=='amazon_nw'){
          nw_nbe <- extract_unc_nbe_all_models(model_variant,region[[1]],i,j,reference,tp)
          row.names(nw_nbe)<-seq(1:nrow(nw_nbe))
        }
        else if (i=='amazon_sw') {
          sw_nbe <- extract_unc_nbe_all_models(model_variant,region[[2]],i,j,reference,tp)
          row.names(sw_nbe)<-seq(1:nrow(sw_nbe))
        }
        else if (i=='amazon_ec') {
          ec_nbe <- extract_unc_nbe_all_models(model_variant,region[[3]],i,j,reference,tp)
          row.names(ec_nbe)<-seq(1:nrow(ec_nbe))
        }
        else if (i=='amazon_bs') {
          bs_nbe <- extract_unc_nbe_all_models(model_variant,region[[4]],i,j,reference,tp)
          row.names(bs_nbe)<-seq(1:nrow(bs_nbe))
        }
        else if (i=='amazon_gs') {
          gs_nbe <- extract_unc_nbe_all_models(model_variant,region[[5]],i,j,reference,tp)
          row.names(gs_nbe)<-seq(1:nrow(gs_nbe))
        }
      }
      nbe_upper<- rbind(nw_nbe,sw_nbe,ec_nbe,bs_nbe,gs_nbe)
    }
  }
  region_pixel_model_data <- cbind(nbe,nbe_lower,nbe_upper)
  return(region_pixel_model_data[,c(2:6,1,7,13)])
  # return(region_pixel_model_data)
}
# all_pixel_nbe <- all_pixel_reg_run(nbe_model_variants[1],nbe_pixel_data,nbe_pixel_names,nbe_var[1],reference_nbe_data,nbe_time_period)
# all_pixel_nbe_unc <- all_pixel_reg_run(nbe_model_variants[2],nbe_pixel_data,nbe_pixel_names,nbe_var,reference_nbe_data,nbe_time_period)

####extract for all model variations
all_nbe_reg_model_var_run <- function(model_variant,region,region_name,variable_name,reference,tp) {
  for (i in model_variant) {
    if(i=='raw_GEOSCHEM'){
      raw_nbe <- all_nbe_reg_run(model_variant[1],region,region_name,variable_name,reference,tp)
    }
    else if(i=='default_CARDAMOM'){
      default_nbe <- all_nbe_reg_run(model_variant[2],region,region_name,variable_name,reference,tp)
    }
    else if (i=='geoschem_CARDAMOM') {
      geoschem_nbe <- all_nbe_reg_run(model_variant[3],region,region_name,variable_name,reference,tp)
    }
  }
  return(rbind(raw_nbe,default_nbe,geoschem_nbe))
}

all_nbe_pixel_model_var_run <- function(model_variant,region,region_name,variable_name,reference,tp) {
  for (i in model_variant) {
    if(i=='raw_GEOSCHEM'){
      raw_nbe <- all_pixel_reg_run(model_variant[1],region,region_name,variable_name,reference,tp)
      raw_nbe$nbe_2pt5pc<-raw_nbe$nbe-abs(raw_nbe$nbe_2pt5pc)
      raw_nbe$nbe_97pt5pc<-raw_nbe$nbe+abs(raw_nbe$nbe_97pt5pc)
    }
    else if(i=='default_CARDAMOM'){
      default_nbe <- all_pixel_reg_run(model_variant[2],region,region_name,variable_name,reference,tp)
    }
    else if (i=='geoschem_CARDAMOM') {
      geoschem_nbe <- all_pixel_reg_run(model_variant[3],region,region_name,variable_name,reference,tp)
    }
  }
  return(rbind(raw_nbe,default_nbe,geoschem_nbe))
}

##########LAI##############
####extract lai data frame
extract_lai_all_models <- function (model_variant,region,region_name,variable_name,reference,tp){
  for (i in model_variant) {
    if (i=='raw_COPERNICUS'){
      data <- lai_stack_01_19_monthly
    } else {
      data <- brick(paste(nbe_prefix,model_variant[1],nbe_suffix,sep=""),varname=variable_name)
    }
    a <- data$X2015.01.01
    cropped <- crop(a, extent(region))
    masked1 <- mask(cropped, region)
    masked2 <- masked1 > -Inf
    poly <- rasterToPolygons(masked2, dissolve=TRUE)
    data_region <- mask(data,poly)
    data_region_1519<-data_region[[which(getZ(data_region) >= as.Date("2015-01-01") & getZ(data_region) <= as.Date("2019-12-01"))]]
    df_data_region_1519 <- as.data.frame(data_region_1519,xy=F)
    df_data_region_1519<-as.data.frame(t(df_data_region_1519))
    df_data_region_1519 <- data.frame(y=unlist(df_data_region_1519))
    names(df_data_region_1519) <- paste(tolower(variable_name))
    regional_name <- as.factor(rep(region_name,nrow(df_data_region_1519)))
    model_variant_name <- as.factor(rep(model_variant,nrow(df_data_region_1519)))
    no_dp<-nrow(df_data_region_1519)/60
    date<- as.Date(rep(seq(from = as.Date("2015-01-01"), to = as.Date("2019-12-01"), by = 'month'),no_dp))
    month<-format(date,format="%b")
    year <- format(date,format="%Y")
  }
  return(na.omit(cbind(df_data_region_1519,regional_name,date,month,year,model_variant_name)))
}

# nw_lai <- extract_lai_all_models(lai_model_variants[1],nbe_reg_data[[1]],nbe_reg_names[1],lai_var[1],reference_nbe_data,nbe_time_period)

#extract lai pixel
extract_unc_lai_all_models <- function (model_variant,region,region_name,variable_name,reference,tp){
  for (i in model_variant) {
    if (i=='raw_COPERNICUS'){
      data <- lai_unc_stack_01_19_monthly
    } else {
      data <- brick(paste(nbe_prefix,model_variant[1],nbe_suffix,sep=""),varname= variable_name)
    }
    a <- data$X2015.01.01
    cropped <- crop(a, extent(region))
    masked1 <- mask(cropped, region)
    masked2 <- masked1 > -Inf
    poly <- rasterToPolygons(masked2, dissolve=TRUE)
    data_region <- mask(data,poly)
    data_region_1519<-data_region[[which(getZ(data_region) >= as.Date("2015-01-01") & getZ(data_region) <= as.Date("2019-12-01"))]]
    df_data_region_1519 <- as.data.frame(data_region_1519,xy=F)
    df_data_region_1519<-as.data.frame(t(df_data_region_1519))
    df_data_region_1519 <- data.frame(y=unlist(df_data_region_1519))
    names(df_data_region_1519) <- paste(tolower(variable_name))
    regional_name <- as.factor(rep(region_name,nrow(df_data_region_1519)))
    model_variant_name <- as.factor(rep(model_variant,nrow(df_data_region_1519)))
    no_dp<-nrow(df_data_region_1519)/60
    date<- as.Date(rep(seq(from = as.Date("2015-01-01"), to = as.Date("2019-12-01"), by = 'month'),no_dp))
    month<-format(date,format="%b")
    year <- format(date,format="%Y")
  }
  return(na.omit(cbind(df_data_region_1519,regional_name,date,month,year,model_variant_name)))
}
# nw_pixel_lai <- extract_lai_all_models(lai_model_variants[1],nbe_pixel_data[[1]],nbe_pixel_names[1],lai_var[1],reference_nbe_data,nbe_time_period)
nw_pixel_lai_unc <- extract_unc_lai_all_models(lai_model_variants[1],nbe_pixel_data[[1]],nbe_pixel_names[1],lai_var[2],reference_nbe_data,nbe_time_period)

####extract for each region
all_lai_reg_run <- function(model_variant,region,region_name,variable_name,reference,tp) {
  for (i in region_name) {
    if(i=='amazon_nw'){
      nw_lai <- extract_lai_all_models(model_variant,region[[1]],i,variable_name,reference,tp)
      row.names(nw_lai)<-seq(1:nrow(nw_lai))
    }
    else if (i=='amazon_sw') {
      sw_lai <- extract_lai_all_models(model_variant,region[[2]],i,variable_name,reference,tp)
      row.names(sw_lai)<-seq(1:nrow(sw_lai))
    }
    else if (i=='amazon_ec') {
      ec_lai <- extract_lai_all_models(model_variant,region[[3]],i,variable_name,reference,tp)
      row.names(ec_lai)<-seq(1:nrow(ec_lai))
    }
    else if (i=='amazon_bs') {
      bs_lai <- extract_lai_all_models(model_variant,region[[4]],i,variable_name,reference,tp)
      row.names(bs_lai)<-seq(1:nrow(bs_lai))
    }
    else if (i=='amazon_gs') {
      gs_lai <- extract_lai_all_models(model_variant,region[[5]],i,variable_name,reference,tp)
      row.names(gs_lai)<-seq(1:nrow(gs_lai))
    }
    else{
      amazon_lai <- extract_lai_all_models(model_variant,region[[6]],i,variable_name,reference,tp)
      row.names(amazon_lai)<-seq(1:nrow(amazon_lai))
    }
  }
  return(rbind(nw_lai,sw_lai,ec_lai,bs_lai,gs_lai,amazon_lai))
}
# all_regional_lai <- all_lai_reg_run(lai_model_variants[1],nbe_reg_data,nbe_reg_names,lai_var[1],reference_nbe_data,nbe_time_period)

all_pixel_reg_run <- function(model_variant,region,region_name,variable_name,reference,tp) {
  for (j in variable_name) {
  if (j=='LAI') {
    for (i in region_name) {
      if(i=='amazon_nw'){
        nw_lai <- extract_lai_all_models(model_variant,region[[1]],i,j,reference,tp)
        row.names(nw_lai)<-seq(1:nrow(nw_lai))
      }
      else if (i=='amazon_sw') {
        sw_lai <- extract_lai_all_models(model_variant,region[[2]],i,j,reference,tp)
        row.names(sw_lai)<-seq(1:nrow(sw_lai))
      }
      else if (i=='amazon_ec') {
        ec_lai <- extract_lai_all_models(model_variant,region[[3]],i,j,reference,tp)
        row.names(ec_lai)<-seq(1:nrow(ec_lai))
      }
      else if (i=='amazon_bs') {
        bs_lai <- extract_lai_all_models(model_variant,region[[4]],i,j,reference,tp)
        row.names(bs_lai)<-seq(1:nrow(bs_lai))
      }
      else if (i=='amazon_gs') {
        gs_lai <- extract_lai_all_models(model_variant,region[[5]],i,j,reference,tp)
        row.names(gs_lai)<-seq(1:nrow(gs_lai))
      }
    }
    lai<- rbind(nw_lai,sw_lai,ec_lai,bs_lai,gs_lai)
  }
  else if (j=='LAI_2pt5pc') {
    for (i in region_name) {
      if(i=='amazon_nw'){
        nw_lai <- extract_unc_lai_all_models(model_variant,region[[1]],i,j,reference,tp)
        row.names(nw_lai)<-seq(1:nrow(nw_lai))
      }
      else if (i=='amazon_sw') {
        sw_lai <- extract_unc_lai_all_models(model_variant,region[[2]],i,j,reference,tp)
        row.names(sw_lai)<-seq(1:nrow(sw_lai))
      }
      else if (i=='amazon_ec') {
        ec_lai <- extract_unc_lai_all_models(model_variant,region[[3]],i,j,reference,tp)
        row.names(ec_lai)<-seq(1:nrow(ec_lai))
      }
      else if (i=='amazon_bs') {
        bs_lai <- extract_unc_lai_all_models(model_variant,region[[4]],i,j,reference,tp)
        row.names(bs_lai)<-seq(1:nrow(bs_lai))
      }
      else if (i=='amazon_gs') {
        gs_lai <- extract_unc_lai_all_models(model_variant,region[[5]],i,j,reference,tp)
        row.names(gs_lai)<-seq(1:nrow(gs_lai))
      }
    }
    lai_lower<- rbind(nw_lai,sw_lai,ec_lai,bs_lai,gs_lai)
  }
  else {
    for (i in region_name) {
      if(i=='amazon_nw'){
        nw_lai <- extract_unc_lai_all_models(model_variant,region[[1]],i,j,reference,tp)
        row.names(nw_lai)<-seq(1:nrow(nw_lai))
      }
      else if (i=='amazon_sw') {
        sw_lai <- extract_unc_lai_all_models(model_variant,region[[2]],i,j,reference,tp)
        row.names(sw_lai)<-seq(1:nrow(sw_lai))
      }
      else if (i=='amazon_ec') {
        ec_lai <- extract_unc_lai_all_models(model_variant,region[[3]],i,j,reference,tp)
        row.names(ec_lai)<-seq(1:nrow(ec_lai))
      }
      else if (i=='amazon_bs') {
        bs_lai <- extract_unc_lai_all_models(model_variant,region[[4]],i,j,reference,tp)
        row.names(bs_lai)<-seq(1:nrow(bs_lai))
      }
      else if (i=='amazon_gs') {
        gs_lai <- extract_unc_lai_all_models(model_variant,region[[5]],i,j,reference,tp)
        row.names(gs_lai)<-seq(1:nrow(gs_lai))
      }
    }
    lai_upper<- rbind(nw_lai,sw_lai,ec_lai,bs_lai,gs_lai)
  }
  }
  region_pixel_model_data <- cbind(lai,lai_lower,lai_upper)
  return(region_pixel_model_data[,c(2:6,1,7,13)])
  # return(region_pixel_model_data)
}

# all_pixel_lai_unc <- all_pixel_reg_run(lai_model_variants[1],nbe_pixel_data,nbe_pixel_names,lai_var,reference_nbe_data,nbe_time_period)

####extract for all model variations
all_lai_reg_model_var_run <- function(model_variant,region,region_name,variable_name,reference,tp) {
  for (i in model_variant) {
    if(i=='raw_COPERNICUS'){
      raw_lai <- all_lai_reg_run(model_variant[1],region,region_name,variable_name,reference,tp)
    }
    else if(i=='default_CARDAMOM'){
      default_lai <- all_lai_reg_run(model_variant[2],region,region_name,variable_name,reference,tp)
    }
    else if (i=='geoschem_CARDAMOM') {
      geoschem_lai <- all_lai_reg_run(model_variant[3],region,region_name,variable_name,reference,tp)
    }
  }
  return(rbind(raw_lai,default_lai,geoschem_lai))
}

all_lai_pixel_model_var_run <- function(model_variant,region,region_name,variable_name,reference,tp) {
  for (i in model_variant) {
    if(i=='raw_COPERNICUS'){
      raw_lai <- all_pixel_reg_run(model_variant[1],region,region_name,variable_name,reference,tp)
      raw_lai$lai_2pt5pc<-raw_lai$lai-abs(raw_lai$lai_2pt5pc)
      raw_lai$lai_97pt5pc<-raw_lai$lai+abs(raw_lai$lai_97pt5pc)
    }
    else if(i=='default_CARDAMOM'){
      default_lai <- all_pixel_reg_run(model_variant[2],region,region_name,variable_name,reference,tp)
    }
    else if (i=='geoschem_CARDAMOM') {
      geoschem_lai <- all_pixel_reg_run(model_variant[3],region,region_name,variable_name,reference,tp)
    }
  }
  return(rbind(raw_lai,default_lai,geoschem_lai))
}

#### Run all analysis#### 
#NBE
all_regional_model_var_nbe <- all_nbe_reg_model_var_run(nbe_model_variants,nbe_reg_data,nbe_reg_names,nbe_var[1],reference_nbe_data,nbe_time_period)
str(all_regional_model_var_nbe)
all_regional_model_var_nbe$month<-factor(all_regional_model_var_nbe$month,levels=c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'))

all_pixels_model_var_nbe_nd_unc <- all_nbe_pixel_model_var_run(nbe_model_variants,nbe_pixel_data,nbe_pixel_names,nbe_var,reference_nbe_data,nbe_time_period)

#LAI
all_regional_model_var_lai <- all_lai_reg_model_var_run(lai_model_variants,nbe_reg_data,nbe_reg_names,lai_var[1],reference_nbe_data,nbe_time_period)
str(all_regional_model_var_lai)
all_regional_model_var_lai$month<-factor(all_regional_model_var_lai$month,levels=c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'))

all_pixels_model_var_lai_nd_unc <- all_lai_pixel_model_var_run(lai_model_variants,nbe_pixel_data,nbe_pixel_names,lai_var,reference_nbe_data,nbe_time_period)
str(all_pixels_model_var_lai_nd_unc)


#### Pixel Data manipulation and Plots #### 
#### NBE
# add date as new_date
all_pixels_model_var_nbe_nd_unc<-all_pixels_model_var_nbe_nd_unc %>% 
  mutate(new_date = ymd(date))

# plot them
ggplot(all_pixels_model_var_nbe_nd_unc[all_pixels_model_var_nbe_nd_unc$regional_name !='amazon_nw',],aes(x = new_date, y = nbe))+
  geom_line(data=all_pixels_model_var_nbe_nd_unc[all_pixels_model_var_nbe_nd_unc$model_variant_name=='default_CARDAMOM'&all_pixels_model_var_nbe_nd_unc$regional_name !='amazon_nw',],position=position_dodge(0.1)) +
  geom_ribbon(data=all_pixels_model_var_nbe_nd_unc[all_pixels_model_var_nbe_nd_unc$model_variant_name=='default_CARDAMOM'&all_pixels_model_var_nbe_nd_unc$regional_name !='amazon_nw',],aes(ymin=nbe_2pt5pc, ymax=nbe_97pt5pc), alpha = 0.3, fill = "red", color = "black", linetype = "dotted") + 
  geom_point(data=all_pixels_model_var_nbe_nd_unc[all_pixels_model_var_nbe_nd_unc$model_variant_name=='raw_GEOSCHEM'&all_pixels_model_var_nbe_nd_unc$regional_name !='amazon_nw',],position=position_dodge(0.1)) +
  geom_errorbar(data=all_pixels_model_var_nbe_nd_unc[all_pixels_model_var_nbe_nd_unc$model_variant_name=='raw_GEOSCHEM'&all_pixels_model_var_nbe_nd_unc$regional_name !='amazon_nw',],aes(ymin=nbe_2pt5pc, ymax=nbe_97pt5pc), width=.1, position=position_dodge(0.1)) +
  facet_wrap(.~regional_name) +
  theme_test()

ggplot(all_pixels_model_var_nbe_nd_unc[all_pixels_model_var_nbe_nd_unc$regional_name !='amazon_nw',],aes(x = new_date, y = nbe))+
  geom_line(data=all_pixels_model_var_nbe_nd_unc[all_pixels_model_var_nbe_nd_unc$model_variant_name=='geoschem_CARDAMOM'&all_pixels_model_var_nbe_nd_unc$regional_name !='amazon_nw',],position=position_dodge(0.1)) +
  geom_ribbon(data=all_pixels_model_var_nbe_nd_unc[all_pixels_model_var_nbe_nd_unc$model_variant_name=='geoschem_CARDAMOM'&all_pixels_model_var_nbe_nd_unc$regional_name !='amazon_nw',],aes(ymin=nbe_2pt5pc, ymax=nbe_97pt5pc), alpha = 0.3, fill = "red", color = "black", linetype = "dotted") + 
  geom_point(data=all_pixels_model_var_nbe_nd_unc[all_pixels_model_var_nbe_nd_unc$model_variant_name=='raw_GEOSCHEM'&all_pixels_model_var_nbe_nd_unc$regional_name !='amazon_nw',],position=position_dodge(0.1)) +
  geom_errorbar(data=all_pixels_model_var_nbe_nd_unc[all_pixels_model_var_nbe_nd_unc$model_variant_name=='raw_GEOSCHEM'&all_pixels_model_var_nbe_nd_unc$regional_name !='amazon_nw',],aes(ymin=nbe_2pt5pc, ymax=nbe_97pt5pc), width=.1, position=position_dodge(0.1)) +
  facet_wrap(.~regional_name) +
  theme_test()

all_pixels_model_var_lai_nd_unc<-all_pixels_model_var_lai_nd_unc %>% 
  mutate(new_date = ymd(date))

#### LAI
# plot them
ggplot(all_pixels_model_var_lai_nd_unc,aes(x = new_date, y = lai))+
  geom_line(data=all_pixels_model_var_lai_nd_unc[all_pixels_model_var_lai_nd_unc$model_variant_name=='default_CARDAMOM',],position=position_dodge(0.1)) +
  geom_ribbon(data=all_pixels_model_var_lai_nd_unc[all_pixels_model_var_lai_nd_unc$model_variant_name=='default_CARDAMOM',],aes(ymin=lai_2pt5pc, ymax=lai_97pt5pc), alpha = 0.3, fill = "red", color = "black", linetype = "dotted") + 
  geom_point(data=all_pixels_model_var_lai_nd_unc[all_pixels_model_var_lai_nd_unc$model_variant_name=='raw_COPERNICUS',],position=position_dodge(0.1)) +
  geom_errorbar(data=all_pixels_model_var_lai_nd_unc[all_pixels_model_var_lai_nd_unc$model_variant_name=='raw_COPERNICUS',],aes(ymin=lai_2pt5pc, ymax=lai_97pt5pc), width=.1, position=position_dodge(0.1)) +
  facet_wrap(.~regional_name) +
  theme_test()

ggplot(all_pixels_model_var_lai_nd_unc,aes(x = new_date, y = lai))+
  geom_line(data=all_pixels_model_var_lai_nd_unc[all_pixels_model_var_lai_nd_unc$model_variant_name=='geoschem_CARDAMOM',],position=position_dodge(0.1)) +
  geom_ribbon(data=all_pixels_model_var_lai_nd_unc[all_pixels_model_var_lai_nd_unc$model_variant_name=='geoschem_CARDAMOM',],aes(ymin=lai_2pt5pc, ymax=lai_97pt5pc), alpha = 0.3, fill = "red", color = "black", linetype = "dotted") + 
  geom_point(data=all_pixels_model_var_lai_nd_unc[all_pixels_model_var_lai_nd_unc$model_variant_name=='raw_COPERNICUS',],position=position_dodge(0.1)) +
  geom_errorbar(data=all_pixels_model_var_lai_nd_unc[all_pixels_model_var_lai_nd_unc$model_variant_name=='raw_COPERNICUS',],aes(ymin=lai_2pt5pc, ymax=lai_97pt5pc), width=.1, position=position_dodge(0.1)) +
  facet_wrap(.~regional_name) +
  theme_test()

#### Summarise regional data####
#summarise data
## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  # library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
#NBE
regional_monthly_nbe_summary<-summarySE(all_regional_model_var_nbe, measurevar="nbe", groupvars=c("regional_name","model_variant_name", "month"))
regional_annual_nbe_summary<-summarySE(all_regional_model_var_nbe, measurevar="nbe", groupvars=c("regional_name","model_variant_name", "year"))

regional_monthly_nbe_summary$month<-factor(regional_monthly_nbe_summary$month,levels=c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'))
regional_annual_nbe_summary$year<-factor(regional_annual_nbe_summary$year)

# Use 95% confidence interval instead of SEM
ggplot(regional_monthly_nbe_summary[regional_monthly_nbe_summary$regional_name=='amazonia',], aes(x=month, y=nbe, fill=model_variant_name))+ 
  geom_bar(position=position_dodge(), stat="identity") + 
  geom_errorbar(aes(ymin=nbe-ci, ymax=nbe+ci),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))+ 
  scale_fill_grey(start = 0.8, end = 0.2) +
  theme_bw() + 
  theme(text = element_text(size = 16),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggplot(regional_annual_nbe_summary[regional_annual_nbe_summary$regional_name=='amazonia',], aes(x=year, y=nbe, fill=model_variant_name))+ 
  geom_bar(position="dodge", stat="identity") + 
  geom_errorbar(aes(ymin=nbe-ci, ymax=nbe+ci),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))+ 
  scale_fill_grey(start = 0.8, end = 0.2) +
  theme_bw() + 
  theme(text = element_text(size = 16),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#all
ggplot(regional_monthly_nbe_summary, aes(x=month, y=nbe, fill=model_variant_name))+ 
  geom_bar(position=position_dodge(), stat="identity") + 
  geom_errorbar(aes(ymin=nbe-ci, ymax=nbe+ci),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))+ 
  scale_fill_grey(start = 0.8, end = 0.2)+
  facet_wrap(.~regional_name) +
  theme_bw() + 
  theme(text = element_text(size = 16),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggplot(regional_annual_nbe_summary, aes(x=year, y=nbe, fill=model_variant_name))+ 
  geom_bar(position="dodge", stat="identity") + 
  geom_errorbar(aes(ymin=nbe-ci, ymax=nbe+ci),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))+ 
  scale_fill_grey(start = 0.8, end = 0.2) +
  facet_wrap(.~regional_name)+
  theme_bw() + 
  theme(text = element_text(size = 16),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#LAI
regional_monthly_lai_summary<-summarySE(all_regional_model_var_lai, measurevar="lai", groupvars=c("regional_name","model_variant_name", "month"))
regional_annual_lai_summary<-summarySE(all_regional_model_var_lai, measurevar="lai", groupvars=c("regional_name","model_variant_name", "year"))

regional_monthly_lai_summary$month<-factor(regional_monthly_lai_summary$month,levels=c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'))
regional_annual_lai_summary$year<-factor(regional_annual_lai_summary$year)

# Use 95% confidence interval instead of SEM
ggplot(regional_monthly_lai_summary[regional_monthly_lai_summary$regional_name=='amazonia',], aes(x=month, y=lai, fill=model_variant_name))+ 
  geom_bar(position=position_dodge(), stat="identity") + 
  geom_errorbar(aes(ymin=lai-ci, ymax=lai+ci),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))+ 
  scale_fill_grey(start = 0.8, end = 0.2) +
  theme_bw() + 
  theme(text = element_text(size = 16),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggplot(regional_annual_lai_summary[regional_annual_lai_summary$regional_name=='amazonia',], aes(x=year, y=lai, fill=model_variant_name))+ 
  geom_bar(position="dodge", stat="identity") + 
  geom_errorbar(aes(ymin=lai-ci, ymax=lai+ci),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))+ 
  scale_fill_grey(start = 0.8, end = 0.2) +
  theme_bw() + 
  theme(text = element_text(size = 16),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#all
ggplot(regional_monthly_lai_summary, aes(x=month, y=lai, fill=model_variant_name))+ 
  geom_bar(position=position_dodge(), stat="identity") + 
  geom_errorbar(aes(ymin=lai-ci, ymax=lai+ci),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))+ 
  scale_fill_grey(start = 0.8, end = 0.2)+
  facet_wrap(.~regional_name) +
  theme_bw() + 
  theme(text = element_text(size = 16),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggplot(regional_annual_lai_summary, aes(x=year, y=lai, fill=model_variant_name))+ 
  geom_bar(position="dodge", stat="identity") + 
  geom_errorbar(aes(ymin=lai-ci, ymax=lai+ci),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))+ 
  scale_fill_grey(start = 0.8, end = 0.2) +
  facet_wrap(.~regional_name)+
  theme_bw() + 
  theme(text = element_text(size = 16),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#### Extract model parameters ####
source('R:/cssp_brazil/cssp_brazil_R/code/subset_pixel_functions_only.R')

prefix_new <- 'M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/'
midfix <- '/RESULTS_PROCESSED/'
suffix_new <- '_stock_flux.RData'

nbe_par_model_variants <- c('no_woody_data_copernicus','nbe_data_alone')
nbe_par_model_variants_abr <- c('default_CARDAMOM','geoschem_CARDAMOM')

par_names <- c('Lit2SOM (day_1)','GPP%Ra','NPP_fol_frac','NPP_root_frac','Leaf lifespan','TO Wood','TO Roots','Mineralise Lit','SOM2Rh',
               'Temp fac','Canopy Eff (gC/m2leaf/day)','Max bud burst day','NPP_lab_frac','Labile_release_period','max leaf fall',
               'leaf fall period','LMA (gC/m2)','Labile_initial','Foliar_initial','Roots_initial','Wood_initial','Lit_initial',
               'SOM_initial','Soil Water_initial','Cwood_coarseR_frac','CoarseR Biomass 50% Max Depth','Max Root Depth',
               'Res factor Burned Cstocks','CCF Fol','CCF Wood and FineR','CCF Soil','CCF Fol and FineR Lit','unknown','model_var')

param_names <- c('model_var','Lit2SOM (day_1)','GPP%Ra','NPP_fol_frac','NPP_root_frac','Leaf lifespan','TO Wood','TO Roots','Mineralise Lit','SOM2Rh',
                 'Temp fac','Canopy Eff (gC/m2leaf/day)','Max bud burst day','NPP_lab_frac','Labile_release_period','max leaf fall',
                 'leaf fall period','LMA (gC/m2)','Labile_initial','Foliar_initial','Roots_initial','Wood_initial','Lit_initial',
                 'SOM_initial','Soil Water_initial','Cwood_coarseR_frac','CoarseR Biomass 50% Max Depth','Max Root Depth',
                 'Res factor Burned Cstocks','CCF Fol','CCF Wood and FineR','CCF Soil','CCF Fol and FineR Lit')

model_pars_df<-as.data.frame(t(data.frame(extract_parameters_all_models(nbe_par_model_variants,nbe_par_model_variants_abr))), stringsAsFactors = FALSE)
model_pars_df[] <- lapply(model_pars_df, type.convert, as.is = TRUE)
rownames(model_pars_df)<-  seq_along(model_pars_df[,1])
colnames(model_pars_df)<-  par_names
model_pars_df$model_var<-as.factor(model_pars_df$model_var)
model_pars_df<-model_pars_df[,c(34,1:32)]
model_pars_df$NPP_wood_frac<-1-(model_pars_df[,4]+model_pars_df[,5])
model_pars_df$rt_wood<-1/model_pars_df[,7]
model_pars_df$rt_roots<-1/model_pars_df[,8]
rownames(model_pars_df)<-nbe_par_model_variants_abr
new_model_pars_df<-as.data.frame(t(model_pars_df)[c(2:34),])
new_model_pars_df$default_CARDAMOM<-as.numeric(new_model_pars_df$default_CARDAMOM);new_model_pars_df$geoschem_CARDAMOM<-as.numeric(new_model_pars_df$geoschem_CARDAMOM)
options(scipen = 100, digits = 4)

# write.csv(new_model_pars_df,"nbe_analysis_model_param.csv")
