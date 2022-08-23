##############################################################
################################models########################
##############################################################
#install packages----
library(ncdf4); library(raster); library(dplyr); library(ggplot2);library(ggpubr);library(quantreg);library(ggpp);library(rgeos);library(ggpmisc);library(rgdal)

#done
#functions for extraction----
extract_var_mean_biom <- function (region,data,reference){
  data <- crop(data, extent(reference))
  a <- data$X2001.01.01
  cropped <- crop(a, extent(region))
  masked1 <- mask(cropped, region)
  masked2 <- masked1 > -Inf
  poly <- rasterToPolygons(masked2, dissolve=TRUE)
  data_region <- mask(data,poly)
  data_2010<-data_region[[which(getZ(data_region) >= as.Date("2010-01-01") & getZ(data_region) <= as.Date("2010-12-01"))]]
  data_2010_mean <- stackApply(data_2010, indices =  rep(1,nlayers(data_2010)), fun = "mean")
  return(data_2010_mean)
}
extract_var_mean <- function (region,data,reference){
  data <- crop(data, extent(reference))
  a <- data$X2001.01.01
  cropped <- crop(a, extent(region))
  masked1 <- mask(cropped, region)
  masked2 <- masked1 > -Inf
  poly <- rasterToPolygons(masked2, dissolve=TRUE)
  data_region <- mask(data,poly)
  data_0110<-data_region[[which(getZ(data_region) >= as.Date("2001-01-01") & getZ(data_region) <= as.Date("2009-12-01"))]]
  data_1016<-data_region[[which(getZ(data_region) >= as.Date("2010-01-01") & getZ(data_region) <= as.Date("2016-12-01"))]]
  data_0110_mean <- stackApply(data_0110, indices =  rep(1,nlayers(data_0110)), fun = "mean")
  data_1016_mean <- stackApply(data_1016, indices =  rep(1,nlayers(data_1016)), fun = "mean")
  return(c(data_0110_mean,data_1016_mean))
}
extract_subset <- function (region,reference){
  masked1 <- mask(reference, region)
  masked2 <- masked1 > -Inf
  poly <- rasterToPolygons(masked2, dissolve=TRUE)
  data_region <- mask(reference,poly)
  return(data_region)
}

extract_var_mean_biom_subset <- function (data){
  data_2010<-data[[which(getZ(data) >= as.Date("2010-01-01") & getZ(data) <= as.Date("2010-12-01"))]]
  data_2010_mean <- stackApply(data_2010, indices =  rep(1,nlayers(data_2010)), fun = "mean")
  return(data_2010_mean)
}
extract_var_mean_subset <- function (data){
  data_0109<-data[[which(getZ(data) >= as.Date("2001-01-01") & getZ(data) <= as.Date("2009-12-01"))]]
  data_1016<-data[[which(getZ(data) >= as.Date("2010-01-01") & getZ(data) <= as.Date("2016-12-01"))]]
  data_0109_mean <- stackApply(data_0109, indices =  rep(1,nlayers(data_0109)), fun = "mean")
  data_1016_mean <- stackApply(data_1016, indices =  rep(1,nlayers(data_1016)), fun = "mean")
  return(c(data_0109_mean,data_1016_mean))
}
#done
#CARDAMOM NORAINFOR----
cardamom_cwood <- brick('R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/MODELS/CARDAMOM_NORAINFOR_new/NoRainfor_Amazon_1deg_monthly_2001_updated_2019.nc',varname="WOOD")
cardamom_nppwood <- brick('R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/MODELS/CARDAMOM_NORAINFOR_new/NoRainfor_Amazon_1deg_monthly_2001_updated_2019.nc',varname="NPP_wood_flx")
cardamom_outputwood <- brick('R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/MODELS/CARDAMOM_NORAINFOR_new/NoRainfor_Amazon_1deg_monthly_2001_updated_2019.nc',varname="OUTPUT_wood_flx")

names(cardamom_cwood)
cardamom_cwood_17 <- cardamom_cwood[[which(getZ(cardamom_cwood) >= as.Date("2017-01-01") & getZ(cardamom_cwood) <= as.Date("2017-12-01"))]]
cardamom_cwood_17_mean <- stackApply(cardamom_cwood_17, indices =  rep(1,nlayers(cardamom_cwood_17)), fun = "mean")

cardamom_cwood_18 <- cardamom_cwood[[which(getZ(cardamom_cwood) >= as.Date("2018-01-01") & getZ(cardamom_cwood) <= as.Date("2018-12-01"))]]
cardamom_cwood_18_mean <- stackApply(cardamom_cwood_18, indices =  rep(1,nlayers(cardamom_cwood_18)), fun = "mean")

cardamom_cwood_0110_mean<-extract_var_mean_biom(amazonia,cardamom_cwood,biomass_amazon_gCm2)

cardamom_nppwood_extract<- extract_var_mean(amazonia,cardamom_nppwood,woodprod_00_09_gCm2d)
cardamom_nppwood_0109_mean<-cardamom_nppwood_extract[[1]]
cardamom_nppwood_1016_mean<-cardamom_nppwood_extract[[2]]

cardamom_outputwood_extract<- extract_var_mean(amazonia,cardamom_outputwood,biommort_00_09_gCm2d)
cardamom_outputwood_0109_mean<-cardamom_outputwood_extract[[1]]
cardamom_outputwood_1016_mean<-cardamom_outputwood_extract[[2]]

##subset
amazonia_subset <- shapefile("./data/amazonia_subset.shp")

cardamom_norainfor_cwood_0110_mean_subset<-extract_var_mean_biom(amazonia_subset,cardamom_cwood,biomass_amazon_gCm2)

cardamom_norainfor_nppwood_extract_subset<- extract_var_mean(amazonia_subset,cardamom_nppwood,woodprod_00_09_gCm2d)
cardamom_norainfor_nppwood_0109_mean_subset<-cardamom_norainfor_nppwood_extract_subset[[1]]
cardamom_norainfor_nppwood_1016_mean_subset<-cardamom_norainfor_nppwood_extract_subset[[2]]

cardamom_norainfor_outputwood_extract_subset<- extract_var_mean(amazonia_subset,cardamom_outputwood,biommort_00_09_gCm2d)
cardamom_norainfor_outputwood_0109_mean_subset<-cardamom_norainfor_outputwood_extract_subset[[1]]
cardamom_norainfor_outputwood_1016_mean_subset<-cardamom_norainfor_outputwood_extract_subset[[2]]

#plot(cardamom_cwood_17_mean)
#CARDAMOM ALL RAINFOR----
cardamom_rainfor_cwood <- stack('R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/DATA/benchmark/Amazon_1deg_monthly_2001_2019.nc',varname="WOOD")
names(cardamom_rainfor_cwood)
cardamom_rainfor_cwood_01 <- cardamom_rainfor_cwood[[which(getZ(cardamom_rainfor_cwood) >= as.Date("2001-01-01") & getZ(cardamom_rainfor_cwood) <= as.Date("2001-12-01"))]]
cardamom_rainfor_cwood_01_mean <- stackApply(cardamom_rainfor_cwood_01, indices =  rep(1,nlayers(cardamom_rainfor_cwood_01)), fun = "mean")

cardamom_rainfor_cwood_0119 <- cardamom_rainfor_cwood[[which(getZ(cardamom_rainfor_cwood) >= as.Date("2001-01-01") & getZ(cardamom_rainfor_cwood) <= as.Date("2019-12-01"))]]
cardamom_rainfor_cwood_0119_mean <- stackApply(cardamom_rainfor_cwood_0119, indices =  rep(1,nlayers(cardamom_rainfor_cwood_0119)), fun = "mean")

#CARDAMOM ALL RAINFOR UPDATED----
cardamom_updated_rainfor_cwood <- brick('R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/DATA/benchmark/Amazon_1deg_monthly_2001_updated_2019.nc',varname="WOOD")
cardamom_updated_rainfor_nppwood <- brick('R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/DATA/benchmark/Amazon_1deg_monthly_2001_updated_2019.nc',varname="NPP_wood_flx")
cardamom_updated_rainfor_outputwood <- brick('R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/DATA/benchmark/Amazon_1deg_monthly_2001_updated_2019.nc',varname="OUTPUT_wood_flx")

testing_extract_var_mean <-extract_var_mean(amazonia,cardamom_updated_rainfor_nppwood,woodprod_00_09_gCm2d)

cardamom_updated_rainfor_cwood_0110_mean<-testing_extract_var_mean[[1]]
cardamom_updated_rainfor_nppwood_0110_mean<-testing_extract_var_mean[[1]]
cardamom_updated_rainfor_nppwood_1016_mean<-testing_extract_var_mean[[2]]
cardamom_updated_rainfor_outputwood_0110_mean<-testing_extract_var_mean[[1]]
cardamom_updated_rainfor_outputwood_1016_mean<-testing_extract_var_mean[[2]]


plot(cardamom_updated_rainfor_outputwood_0110_mean)
plot(biommort_00_09_gCm2d)

cardamom_rainfor_cwood_01 <- cardamom_rainfor_cwood[[which(getZ(cardamom_rainfor_cwood) >= as.Date("2001-01-01") & getZ(cardamom_rainfor_cwood) <= as.Date("2001-12-01"))]]
cardamom_rainfor_cwood_01_mean <- stackApply(cardamom_rainfor_cwood_01, indices =  rep(1,nlayers(cardamom_rainfor_cwood_01)), fun = "mean")

cardamom_rainfor_cwood_0119 <- cardamom_rainfor_cwood[[which(getZ(cardamom_rainfor_cwood) >= as.Date("2001-01-01") & getZ(cardamom_rainfor_cwood) <= as.Date("2019-12-01"))]]
cardamom_rainfor_cwood_0119_mean <- stackApply(cardamom_rainfor_cwood_0119, indices =  rep(1,nlayers(cardamom_rainfor_cwood_0119)), fun = "mean")

#done
#Cwood ESA source AGB----
esa_cci_agb_17 <- brick('G://AGB/ESA_CCI_BIOMASS/ESA_CCI_AGB_1deg/AGB_map_MgCha_2017.tif')
#plot(esa_cci_agb_17)
esa_cci_agb_17_amazon <- crop(esa_cci_agb_17, extent(cardamom_cwood_17_mean))
#plot(esa_cci_agb_17_amazon)

check_extent <- extent(esa_cci_agb_17_amazon) #use south america extent from INLAND SA
check_extent_r <- raster(check_extent)
res(check_extent_r) <- res(esa_cci_agb_17_amazon)
values(check_extent_r) <- 1
crs(check_extent_r) <- "+proj=longlat +datum=WGS84 +no_defs"
check_extent_r_amazonia <- mask(check_extent_r, amazonia_poly, updatevalue=0)
check_extent_r_amazonia[check_extent_r_amazonia<1] <- NA

esa_cci_agb_17_amazon_masked <- mask(esa_cci_agb_17_amazon, check_extent_r_amazonia)
esa_cci_agb_17_amazon_masked <- tha_to_gm2_fun(esa_cci_agb_17_amazon_masked)
#plot(esa_cci_agb_17_amazon_masked)

esa_cci_agb_18 <- brick('G://AGB/ESA_CCI_BIOMASS/ESA_CCI_AGB_1deg/AGB_map_MgCha_2018.tif')
#plot(esa_cci_agb_18)
esa_cci_agb_18_amazon <- crop(esa_cci_agb_18, extent(cardamom_cwood_18_mean))
#plot(esa_cci_agb_18_amazon)

check_extent <- extent(esa_cci_agb_18_amazon) #use south america extent from INLAND SA
check_extent_r <- raster(check_extent)
res(check_extent_r) <- res(esa_cci_agb_18_amazon)
values(check_extent_r) <- 1
crs(check_extent_r) <- "+proj=longlat +datum=WGS84 +no_defs"
check_extent_r_amazonia <- mask(check_extent_r, amazonia_poly, updatevalue=0)
check_extent_r_amazonia[check_extent_r_amazonia<1] <- NA

esa_cci_agb_18_amazon_masked <- mask(esa_cci_agb_18_amazon, check_extent_r_amazonia)
esa_cci_agb_18_amazon_masked <- tha_to_gm2_fun(esa_cci_agb_18_amazon_masked)


#cwood rainfor cwood initial agb----
rainfor_agb_01 <- biomass_amazon_gCm2_card

#RAINFOR raw data----
biomass_amazon_gCm2
woodprod_00_09_gCm2d
woodprod_10_16_gCm2d
biommort_00_09_gCm2d
biommort_10_16_gCm2d
#copy RAINFOR data back for CARDAMOM assimilation----
###
## Wood biomass
# Write back out with new name and units etc
writeRaster(biomass_amazon_gCm2, filename = "G://cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/1deg_converted/wood_biomass/wood_biomass_gCm2.tif", format = "GTiff")
# Write an estimate of the uncertainty back out
writeRaster(biomass_amazon_gCm2*0.25, filename = "G://cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/1deg_converted/wood_biomass/unc_wood_biomass_gCm2.tif", format = "GTiff")

###
## Wood productivity
## Write back out with new name and units etc
writeRaster(woodprod_00_09_gCm2d, filename = "G://cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/1deg_converted/wood_productivity/wood_productivity_gCm2_2000_2009.tif", format = "GTiff")
## Write an estimate of the uncertainty back out
writeRaster(woodprod_00_09_gCm2d*0.25, filename = "G://cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/1deg_converted/wood_productivity/unc_wood_productivity_gCm2_2000_2009.tif", format = "GTiff")

## Write back out with new name and units etc
writeRaster(woodprod_10_16_gCm2d, filename = "G://cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/1deg_converted/wood_productivity/wood_productivity_gCm2_2010_2016.tif", format = "GTiff")
## Write an estimate of the uncertainty back out
writeRaster(woodprod_10_16_gCm2d*0.25, filename = "G://cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/1deg_converted/wood_productivity/unc_wood_productivity_gCm2_2010_2016.tif", format = "GTiff")

###
## Wood mortality

# # Write back out with new name and units etc
writeRaster(biommort_00_09_gCm2d, filename = "G://cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/1deg_converted/wood_mortality/wood_mortality_gCm2_2000_2009.tif", format = "GTiff")
# # Write an estimate of the uncertainty back out
writeRaster(biommort_00_09_gCm2d*0.25, filename = "G://cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/1deg_converted/wood_mortality/unc_wood_mortality_gCm2_2000_2009.tif", format = "GTiff")
#
# # Write back out with new name and units etc
writeRaster(biommort_10_16_gCm2d, filename = "G://cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/1deg_converted/wood_mortality/wood_mortality_gCm2_2010_2016.tif", format = "GTiff")
# # Write an estimate of the uncertainty back out
writeRaster(biommort_10_16_gCm2d*0.25, filename = "G://cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/1deg_converted/wood_mortality/unc_wood_mortality_gCm2_2010_2016.tif", format = "GTiff")
#

#done
#done
#resample CARD_RAINFOR and RAINFOR data----
res_df_merge_plot_cwood <- function (benchmark,model){
  model_res <- resample(model,benchmark)
  model_res_df <- as.data.frame(model_res, xy=TRUE)
  benchmark_df <- as.data.frame(benchmark, xy=TRUE)
  model_benchmark_df <- merge(benchmark_df, model_res_df, by=c("x","y"))
  names(model_benchmark_df)<-c('x','y','rainfor','cardamom')
  model_benchmark_df <- model_benchmark_df[!is.na(model_benchmark_df$rainfor)&!is.na(model_benchmark_df$cardamom),]
  plot(model_benchmark_df$cardamom, model_benchmark_df$rainfor, main="AGB comparison 2001-2010",
       xlab="CARDAMOM_RAINFOR Wood biomass g.m-2 ", ylab="RAINFOR Wood biomass g.m-2", pch=19)
  abline(coef = c(0,1),col='red', lwd=3)
}

res_df_merge_plot_nppwood <- function (benchmark,model){
  model_res <- resample(model,benchmark)
  model_res_df <- as.data.frame(model_res, xy=TRUE)
  benchmark_df <- as.data.frame(benchmark, xy=TRUE)
  model_benchmark_df <- merge(benchmark_df, model_res_df, by=c("x","y"))
  names(model_benchmark_df)<-c('x','y','rainfor','cardamom')
  model_benchmark_df <- model_benchmark_df[!is.na(model_benchmark_df$rainfor)&!is.na(model_benchmark_df$cardamom),]
  plot(model_benchmark_df$cardamom, model_benchmark_df$rainfor, main="Woody NPP comparison 2000/1-2009",
       xlab="CARDAMOM_RAINFOR Woody NPP g.m-2.d-1 ", ylab="RAINFOR Wood productivity g.m-2.d-1", pch=19)
  abline(coef = c(0,1),col='red', lwd=3)
}

res_df_merge_plot_nppwood_2 <- function (benchmark,model){
  model_res <- resample(model,benchmark)
  model_res_df <- as.data.frame(model_res, xy=TRUE)
  benchmark_df <- as.data.frame(benchmark, xy=TRUE)
  model_benchmark_df <- merge(benchmark_df, model_res_df, by=c("x","y"))
  names(model_benchmark_df)<-c('x','y','rainfor','cardamom')
  model_benchmark_df <- model_benchmark_df[!is.na(model_benchmark_df$rainfor)&!is.na(model_benchmark_df$cardamom),]
  plot(model_benchmark_df$cardamom, model_benchmark_df$rainfor, main="Woody NPP comparison 2010-2016",
       xlab="CARDAMOM_RAINFOR Woody NPP g.m-2.d-1 ", ylab="RAINFOR Wood productivity g.m-2.d-1", pch=19)
  abline(coef = c(0,1),col='red', lwd=3)
}

res_df_merge_plot_outputwood <- function (benchmark,model){
  model_res <- resample(model,benchmark)
  model_res_df <- as.data.frame(model_res, xy=TRUE)
  benchmark_df <- as.data.frame(benchmark, xy=TRUE)
  model_benchmark_df <- merge(benchmark_df, model_res_df, by=c("x","y"))
  names(model_benchmark_df)<-c('x','y','rainfor','cardamom')
  model_benchmark_df <- model_benchmark_df[!is.na(model_benchmark_df$rainfor)&!is.na(model_benchmark_df$cardamom),]
  plot(model_benchmark_df$cardamom, model_benchmark_df$rainfor, main="Woody Carbon losses comparison 2000/1-2009",
       xlab="CARDAMOM_RAINFOR Output Wood g.m-2.d-1 ", ylab="RAINFOR Wood mortality g.m-2.d-1", pch=19)
  abline(coef = c(0,1),col='red', lwd=3)
}
res_df_merge_plot_outputwood_2 <- function (benchmark,model){
  model_res <- resample(model,benchmark)
  model_res_df <- as.data.frame(model_res, xy=TRUE)
  benchmark_df <- as.data.frame(benchmark, xy=TRUE)
  model_benchmark_df <- merge(benchmark_df, model_res_df, by=c("x","y"))
  names(model_benchmark_df)<-c('x','y','rainfor','cardamom')
  model_benchmark_df <- model_benchmark_df[!is.na(model_benchmark_df$rainfor)&!is.na(model_benchmark_df$cardamom),]
  plot(model_benchmark_df$cardamom, model_benchmark_df$rainfor, main="Woody Carbon losses comparison 2010-2016",
       xlab="CARDAMOM_RAINFOR Output Wood g.m-2.d-1 ", ylab="RAINFOR Wood mortality g.m-2.d-1", pch=19)
  abline(coef = c(0,1),col='red', lwd=3)
}
# plot(model_benchmark_df$cardamom, model_benchmark_df$rainfor, main="AGB comparison 2001-2010",
#      xlab="CARDAMOM_RAINFOR Wood biomass g.m-2 ", ylab="RAINFOR Wood biomass g.m-2", pch=19)
# abline(coef = c(0,1),col='red', lwd=3)
# plot(model_benchmark_df$cardamom, model_benchmark_df$rainfor, main="Woody NPP comparison 2000/1-2010",
#      xlab="CARDAMOM_RAINFOR Woody NPP g.m-2.d-1 ", ylab="RAINFOR Wood productivity g.m-2.d-1", pch=19)
# abline(coef = c(0,1),col='red', lwd=3)

res_df_merge_plot_cwood(biomass_amazon_gCm2,cardamom_updated_rainfor_cwood_0110_mean)
res_df_merge_plot_nppwood(woodprod_00_09_gCm2d,cardamom_updated_rainfor_nppwood_0110_mean)
res_df_merge_plot_nppwood_2(woodprod_10_16_gCm2d,cardamom_updated_rainfor_nppwood_1016_mean)
res_df_merge_plot_outputwood(biommort_00_09_gCm2d,cardamom_updated_rainfor_outputwood_0110_mean)
res_df_merge_plot_outputwood_2(biommort_10_16_gCm2d,cardamom_updated_rainfor_outputwood_1016_mean)

#done
#resample and extract data----
cardamom_cwood_17_mean_b <- resample(cardamom_cwood_17_mean,esa_cci_agb_17_amazon_masked)
#cardamom_cwood_17_mean_b_df<-as.data.frame(cardamom_cwood_17_mean_b, xy=TRUE)
cardamom_cwood_17_mean_b_df<-as.data.frame(cardamom_cwood_17_mean, xy=TRUE)

esa_cci_agb_17_amazon_masked_b <- resample(esa_cci_agb_17_amazon_masked, cardamom_cwood_17_mean)
plot(esa_cci_agb_17_amazon_masked_b)
#esa_cci_agb_17_amazon_masked_df<-as.data.frame(esa_cci_agb_17_amazon_masked, xy=TRUE)
esa_cci_agb_17_amazon_masked_df<-as.data.frame(esa_cci_agb_17_amazon_masked_b, xy=TRUE)

#2018
cardamom_cwood_18_mean_b <- resample(cardamom_cwood_18_mean,esa_cci_agb_18_amazon_masked)
cardamom_cwood_18_mean_b_df<-as.data.frame(cardamom_cwood_18_mean_b, xy=TRUE)
esa_cci_agb_18_amazon_masked_df<-as.data.frame(esa_cci_agb_18_amazon_masked, xy=TRUE)


#RAINFOR
extent(rainfor_agb_01)<-extent(cardamom_rainfor_cwood_01_mean)
rainfor_agb_01_b <- resample(rainfor_agb_01,cardamom_rainfor_cwood_01_mean)
rainfor_agb_01_mmean_df<-as.data.frame(rainfor_agb_01_b, xy=TRUE)

#cardamom_rainfor 2001
cardamom_rainfor_cwood_01_mean_df<-as.data.frame(cardamom_rainfor_cwood_01_mean, xy=TRUE)

#cardamom_rainfor 2001 2019
cardamom_rainfor_cwood_0119_mean_df<-as.data.frame(cardamom_rainfor_cwood_0119_mean, xy=TRUE)

#done
#merge data frame----
esa_cardamom_agb_17<-merge(cardamom_cwood_17_mean_b_df, esa_cci_agb_17_amazon_masked_df, by=c("x","y"))
names(esa_cardamom_agb)<-c('x','y','cardamom_cw','esa_cci_cw')
esa_cardamom_agb<-esa_cardamom_agb[!is.na(esa_cardamom_agb$cardamom_cw)&!is.na(esa_cardamom_agb$esa_cci_cw),]

esa_cardamom_agb_18<-merge(cardamom_cwood_18_mean_b_df, esa_cci_agb_18_amazon_masked_df, by=c("x","y"))
names(esa_cardamom_agb_18)<-c('x','y','cardamom_cw','esa_cci_cw')
esa_cardamom_agb_18<-esa_cardamom_agb_18[!is.na(esa_cardamom_agb_18$cardamom_cw)&!is.na(esa_cardamom_agb_18$esa_cci_cw),]

rainfor_cardamom_agb_01<-cardamom_rainfor_cwood_01_mean_df
names(rainfor_cardamom_agb_01)<-c('x','y','cardamom_rainfor_cw')
rainfor_cardamom_agb_01$rainfor_cw<-rainfor_agb_01_mmean_df$wood_biomass_gCm2

rainfor_cardamom_agb_01<-rainfor_cardamom_agb_01[!is.na(rainfor_cardamom_agb_01$cardamom_rainfor_cw)&!is.na(rainfor_cardamom_agb_01$rainfor_cw),]

rainfor_cardamom_agb_01<-merge(cardamom_rainfor_cwood_01_mean_df, rainfor_agb_01_mmean_df, by=c("x","y"))
names(rainfor_cardamom_agb_01)<-c('x','y','cardamom_rainfor_cw','rainfor_cw')
rainfor_cardamom_agb_01$cardamom_rainfor_cw_0119 <- cardamom_rainfor_cwood_0119_mean_df$index_1
rainfor_cardamom_agb_01<-rainfor_cardamom_agb_01[!is.na(rainfor_cardamom_agb_01$cardamom_rainfor_cw)&!is.na(rainfor_cardamom_agb_01$rainfor_cw)&!is.na(rainfor_cardamom_agb_01$cardamom_rainfor_cw_0119),]

#done
#plot maps----
plot(esa_cardamom_agb$cardamom_cw, esa_cardamom_agb$esa_cci_cw, main="Biomass comparison 2017",
     xlab="CARDAMOM_NORAINFOR Wood biomass g.m-2 ", ylab="ESA CCI Wood biomass g.m-2", pch=19)
abline(coef = c(0,1),col='red', lwd=3)

plot(esa_cardamom_agb_18$cardamom_cw, esa_cardamom_agb_18$esa_cci_cw, main="Biomass comparison 2018",
     xlab="CARDAMOM_NORAINFOR Wood biomass g.m-2 ", ylab="ESA CCI Wood biomass g.m-2", pch=19)
abline(coef = c(0,1),col='red', lwd=3)

plot(rainfor_cardamom_agb_01$cardamom_rainfor_cw, rainfor_cardamom_agb_01$rainfor_cw, ylim = c(0, 20000), main="Biomass comparison 2001",
     xlab="CARDAMOM_RAINFOR Wood biomass g.m-2 ", ylab="RAINFOR Wood biomass g.m-2", pch=19)
abline(coef = c(0,1),col='red', lwd=3)

plot(rainfor_cardamom_agb_01$cardamom_rainfor_cw_0119, rainfor_cardamom_agb_01$rainfor_cw, ylim = c(0, 20000), main="Biomass comparison 2001 2019",
     xlab="CARDAMOM_RAINFOR Wood biomass 2001-19 g.m-2 ", ylab="RAINFOR Wood biomass g.m-2", pch=19)
abline(coef = c(0,1),col='red', lwd=3)

#done
#save relevant data----
save(cardamom_cwood_0110_mean,cardamom_nppwood_0109_mean,cardamom_nppwood_1016_mean,cardamom_outputwood_0109_mean,
     cardamom_outputwood_1016_mean,cardamom_updated_rainfor_cwood_0110_mean,cardamom_updated_rainfor_nppwood_0110_mean,
     cardamom_updated_rainfor_nppwood_1016_mean,cardamom_updated_rainfor_outputwood_0110_mean,
     cardamom_updated_rainfor_outputwood_1016_mean,file='./data/extracted_card.RData')

###################################################
###############Compare Specific pixels#############
###################################################
##extract specifc pixels
cardamom_updated_rainfor_cwood_pixels <- brick('R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/DATA/benchmark/Amazon_subset_biomass_1deg_monthly_2001_updated_2019.nc',varname="WOOD")
cardamom_updated_rainfor_nppwood_pixels <- brick('R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/DATA/benchmark/Amazon_subset_productivity_1deg_monthly_2001_updated_2019.nc',varname="NPP_wood_flx")
cardamom_updated_rainfor_outputwood_pixels <- brick('R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/DATA/benchmark/Amazon_subset_mortality_1deg_monthly_2001_updated_2019.nc',varname="OUTPUT_wood_flx")

biomass_amazon_gCm2_pixels<-extract_subset(amazonia_subset,biomass_amazon_gCm2)
woodprod_00_09_gCm2d_pixels<-extract_subset(amazonia_subset,woodprod_00_09_gCm2d)
woodprod_10_16_gCm2d_pixels<-extract_subset(amazonia_subset,woodprod_10_16_gCm2d)
biommort_00_09_gCm2d_pixels<-extract_subset(amazonia_subset,biommort_00_09_gCm2d)
biommort_10_16_gCm2d_pixels<-extract_subset(amazonia_subset,biommort_10_16_gCm2d)

cardamom_cwood_0110_mean_subset<-extract_var_mean_biom_subset(cardamom_updated_rainfor_cwood_pixels)

cardamom_nppwood_extract_subset<- extract_var_mean_subset(cardamom_updated_rainfor_nppwood_pixels)
cardamom_nppwood_0109_mean_subset<-cardamom_nppwood_extract_subset[[1]]
cardamom_nppwood_1016_mean_subset<-cardamom_nppwood_extract_subset[[2]]

cardamom_outputwood_extract_subset<- extract_var_mean_subset(cardamom_updated_rainfor_outputwood_pixels)
cardamom_outputwood_0109_mean_subset<-cardamom_outputwood_extract_subset[[1]]
cardamom_outputwood_1016_mean_subset<-cardamom_outputwood_extract_subset[[2]]

par(mfrow = c(1, 2))
res_df_merge_plot_cwood(biomass_amazon_gCm2_pixels,cardamom_cwood_0110_mean_subset)
res_df_merge_plot_nppwood(woodprod_00_09_gCm2d_pixels,cardamom_nppwood_0109_mean_subset)
res_df_merge_plot_nppwood_2(woodprod_10_16_gCm2d_pixels,cardamom_nppwood_1016_mean_subset)
res_df_merge_plot_outputwood(biommort_00_09_gCm2d_pixels,cardamom_outputwood_0109_mean_subset)
res_df_merge_plot_outputwood_2(biommort_10_16_gCm2d_pixels,cardamom_outputwood_1016_mean_subset)

res_df_merge_plot_cwood(biomass_amazon_gCm2_pixels,cardamom_norainfor_cwood_0110_mean_subset)
res_df_merge_plot_nppwood(woodprod_00_09_gCm2d_pixels,cardamom_norainfor_nppwood_0109_mean_subset)
res_df_merge_plot_nppwood_2(woodprod_10_16_gCm2d_pixels,cardamom_norainfor_nppwood_1016_mean_subset)
res_df_merge_plot_outputwood(biommort_00_09_gCm2d_pixels,cardamom_norainfor_outputwood_0109_mean_subset)
res_df_merge_plot_outputwood_2(biommort_10_16_gCm2d_pixels,cardamom_norainfor_outputwood_1016_mean_subset)
