#All codes being Used
#install packages----
library(ncdf4); library(raster); library(dplyr); library(ggplot2);library(ggpubr);library(quantreg);library(ggpp);library(rgeos);library(ggpmisc);library(rgdal);library(Metrics);library(gridExtra);library(grid)

#done
######Read in shape files ###########
compare_var_names <- c('WOOD','NPP_wood_flx','OUTPUT_wood_flx','MTT_wood') #parameters for comparisons

## load shape files of Amazonian extent and sub regions [found in shape_files_for_analysis]
amazonia_nw <- shapefile("/exports/csce/datastore/geos/groups/gcel/cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/amazon_nw.shp")
amazonia_sw <- shapefile("/exports/csce/datastore/geos/groups/gcel/cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/amazon_sw.shp")
amazonia_ec <- shapefile("/exports/csce/datastore/geos/groups/gcel/cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/amazon_ec.shp")
amazonia_bs <- shapefile("/exports/csce/datastore/geos/groups/gcel/cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/amazon_bs.shp")
amazonia_gs <- shapefile("/exports/csce/datastore/geos/groups/gcel/cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/amazon_gs.shp")
amazonia_ifl_layer <- shapefile("R:/brazil_leeds_maps/ifl_2000_amazonia.shp") #import Intact forest layers across amazon source: Potapov, P., Yaroshenko, A., Turubanova, S., Dubinin, M., Laestadius, L., Thies, C., Aksenov, D., Egorov, A., Yesipova, Y., Glushkov, I., Karpachevskiy, M., Kostikova, A., Manisha A., Tsybikova, E. & Zhuravleva, I. (2008) Mapping the world’s intact forest landscapes by remote sensing. Ecology and Society, 13(2), 51.
amazonia_extent <- shapefile("./data/amazonia_extent") 

#read in extrapolated RAINFOR data (contact university of leeds for data David Galbraith- D.R.Galbraith@leeds.ac.uk)
#wood productivity t/ha/year
woodprod_00_09 <- brick('R://brazil_leeds_maps/WoodyProductivity20002009_Mg_perHa_perYear_111km.tif')
#biomass mortality t/ha/year
biommort_00_09 <- brick('R://brazil_leeds_maps/BiomassMortality_20002009_Mg_perHa_perYear_111km.tif')
#biomass t/ha
biomass_amazon <- brick('R://brazil_leeds_maps/AbovegroundBiomass_Mg_perHa_111km.tif')

#functions for modification
thayr_to_gCm2day_fun <- function(x) { #convert t/ha/year to gC/m2/day
  x *0.48 * (100/365.25) }
tha_to_gCm2_fun <- function(x) { #convert t/ha to gC/m2
  x *0.48 * 100 }
bgb_inclusion <- function(x){ #include belowground biomass
  x+(0.489 * x ** 0.89)
}

#include BGB to wood dynamics
woodprod_00_09 <- calc (woodprod_00_09,bgb_inclusion)
woodprod_10_16 <- calc (woodprod_10_16,bgb_inclusion)
biommort_00_09 <- calc (biommort_00_09,bgb_inclusion)
biommort_10_16 <- calc (biommort_10_16,bgb_inclusion)
biomass_amazon <- calc (biomass_amazon,bgb_inclusion)
rt_amazon_00_09_year <- biomass_amazon/biommort_00_09

#convert to gcm2
biomass_amazon_gCm2 <- calc(biomass_amazon, tha_to_gCm2_fun)
woodprod_00_09_gCm2d <- calc(woodprod_00_09, thayr_to_gCm2day_fun)
biommort_00_09_gCm2d <- calc(biommort_00_09, thayr_to_gCm2day_fun)

reg_data <- list(amazonia_nw,amazonia_sw,amazonia_ec,amazonia_bs,amazonia_gs,amazonia) #create a list of regional shape files
reg_names <- c('amazon_nw','amazon_sw','amazon_ec','amazon_bs','amazon_gs','amazonia') #names of regions
prefix <- 'R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/amazonia_ifl_' #prefix of file location edit as necessary
suffix <- '_nomngt_2000_updated_2010.nc' #suffix of file name and location
updated_mod_var <- c('esa_cci_agb','rainfor_biomass_annual','rainfor_biomass_productivity_2005','rainfor_biomass_annual_productivity') #CARDAMOM run variants

reference_data<- c(biomass_amazon_gCm2,woodprod_00_09_gCm2d,biommort_00_09_gCm2d,rt_amazon_00_09_year) #create list of benchmark data for comparison

#extract subset of raster data (reference) based on a regional shapefile (region)
extract_subset <- function (region,reference){
  masked1 <- mask(reference, region)
  masked2 <- masked1 > -Inf
  poly <- rasterToPolygons(masked2, dissolve=TRUE)
  data_region <- mask(reference,poly)
  return(data_region)
}

#create subset of RAINFOR data based on Intact forest layers
reference_data_ifl <- c(extract_subset(amazonia_ifl_layer,reference_data[[1]]),
                        extract_subset(amazonia_ifl_layer,reference_data[[2]]),
                        extract_subset(amazonia_ifl_layer,reference_data[[3]]),
                        extract_subset(amazonia_ifl_layer,reference_data[[4]]))

######function to extract RAINFOR amazon regions and correct units#######
extract_benchmark_ecoregion <- function (cardamom_var,reference,ecoregions){
  masked1 <- mask(reference, ecoregions)
  masked2 <- masked1 > -Inf
  poly <- rasterToPolygons(masked2, dissolve=TRUE)
  data_region <- mask(reference,poly)
  if (cardamom_var == 'MTT_wood') {
    data_region <- data_region
  }
  else if (cardamom_var == 'WOOD'){
    data_region <- data_region/100
  }
  else {
    data_region <- data_region * (365.25/100)
  }
  return(data_region)
}

#extract model variants and parameters
#reion: amazonia_ifl_layer (amazonian ifl layer)/cardamom_var: compare_var_names (cardamom paramters for comparison)
#model_variant: updated_mod_var (CARDAMOM run variants), reference: reference_data_ifl (RAINFOR extrapolated data for comparison)
#ecoregions: reg_data(amazonian ecoregions), ecoregions_names: reg_names (names of amazonian ecoregions)
extract_models_ifl <- function (region,cardamom_var,model_variant,reference,ecoregions,ecoregions_names){
  for (i in model_variant) {
    data <- brick(paste(prefix,i,suffix,sep=""),varname=cardamom_var)
    data_mean <- stackApply(data, indices =  rep(1,nlayers(data)), fun = "mean")
    masked1 <- mask(data_mean, region)
    masked2 <- masked1 > -Inf
    poly <- rasterToPolygons(masked2, dissolve=TRUE)
    data_ifl <- mask(data_mean,poly)
  }
  return(data_ifl)
}

extract_model_ecoregion <- function (region,cardamom_var,model_variant,reference,ecoregions,ecoregions_names){
  model_ifl <- extract_models_ifl(region,cardamom_var,model_variant,reference,ecoregions,ecoregions_names)
  model_res <- resample(model_ifl,reference)
  masked1 <- mask(model_res, ecoregions)
  masked2 <- masked1 > -Inf
  poly <- rasterToPolygons(masked2, dissolve=TRUE)
  model_ecoregion <- mask(model_res,poly)
  if (cardamom_var == 'MTT_wood') {
    model_ecoregion <- model_ecoregion
  }
  else if (cardamom_var == 'WOOD'){
    model_ecoregion <- model_ecoregion/100
  }
  else {
    model_ecoregion <- model_ecoregion * (365.25/100)
  }
  return(model_ecoregion)
}

merge_pixel_var_bench_mod <- function (region,cardamom_var,model_variant,reference,ecoregions,ecoregions_names) {
  benchmark <- extract_benchmark_ecoregion(cardamom_var,reference,ecoregions)
  df_benchmark <- as.data.frame(benchmark, xy=TRUE)
  df_benchmark$ecoregion <- as.factor(rep(ecoregions_names,nrow(df_benchmark)))
  rownames(df_benchmark)<-1:nrow(df_benchmark)
  for (i in model_variant) {
    if (i=='esa_cci_agb') {
      model_er <- extract_model_ecoregion(region,cardamom_var,i,reference,ecoregions,ecoregions_names)
      c_res_df <- as.data.frame(model_er, xy=TRUE)
      rownames(c_res_df)<-1:nrow(c_res_df)
      model_c_benchmark_df <- merge(df_benchmark, c_res_df, by=c("x","y"))
      names(model_c_benchmark_df) <- c('x','y','benchmark','ecoregion','cci_esa')
    }
    else if (i=='rainfor_biomass_annual') {
      model_er <- extract_model_ecoregion(region,cardamom_var,i,reference,ecoregions,ecoregions_names)
      r0_res_df <- as.data.frame(model_er, xy=TRUE)
      rownames(r0_res_df)<-1:nrow(r0_res_df)
      model_r0_benchmark_df <- merge(df_benchmark, r0_res_df, by=c("x","y"))
      names(model_r0_benchmark_df) <- c('x','y','benchmark','ecoregion','rainfor_biomass_annual')
    }
    else if (i=='rainfor_biomass_productivity_2005') {
      model_er <- extract_model_ecoregion(region,cardamom_var,i,reference,ecoregions,ecoregions_names)
      r1_res_df <- as.data.frame(model_er, xy=TRUE)
      rownames(r1_res_df)<-1:nrow(r1_res_df)
      model_r1_benchmark_df <- merge(df_benchmark, r1_res_df, by=c("x","y"))
      names(model_r1_benchmark_df) <- c('x','y','benchmark','ecoregion','rainfor_2005')
    }
    else if (i=='rainfor_biomass_annual_productivity') {
      model_er <- extract_model_ecoregion(region,cardamom_var,i,reference,ecoregions,ecoregions_names)
      r2_res_df <- as.data.frame(model_er, xy=TRUE)
      rownames(r2_res_df)<-1:nrow(r2_res_df)
      model_r2_benchmark_df <- merge(df_benchmark, r2_res_df, by=c("x","y"))
      names(model_r2_benchmark_df) <- c('x','y','benchmark','ecoregion','rainfor_annual')
    }
  }
  model_benchmark_df <- cbind(model_c_benchmark_df,model_r0_benchmark_df,model_r1_benchmark_df,model_r2_benchmark_df)
  model_benchmark_df <- model_benchmark_df[,c(1,2,4,3,5,10,15,20)]
  model_benchmark_df <- na.omit(model_benchmark_df)
  rownames(model_benchmark_df)<-1:nrow(model_benchmark_df)
  return(model_benchmark_df)
}

merge_pixel_var_bench_mod_ecoregions <- function (region,cardamom_var,model_variant,reference,ecoregions,ecoregions_names) {
  for (i in ecoregions_names) {
    if (i=='amazon_nw') {
      nw_data <- merge_pixel_var_bench_mod(region,cardamom_var,model_variant,reference,ecoregions[[1]],i)
    }
    else if (i=='amazon_sw') {
      sw_data <- merge_pixel_var_bench_mod(region,cardamom_var,model_variant,reference,ecoregions[[2]],i)
    }
    else if (i=='amazon_ec') {
      ec_data <- merge_pixel_var_bench_mod(region,cardamom_var,model_variant,reference,ecoregions[[3]],i)
    }
    else if (i=='amazon_bs') {
      bs_data <- merge_pixel_var_bench_mod(region,cardamom_var,model_variant,reference,ecoregions[[4]],i)
    }
    else if (i=='amazon_gs') {
      gs_data <- merge_pixel_var_bench_mod(region,cardamom_var,model_variant,reference,ecoregions[[5]],i)
    }
  }
  return(rbind(nw_data,sw_data,ec_data,bs_data,gs_data))
}

new_function_to_plot_them_all <- function(region,cardamom_var,model_variant,reference,ecoregions,ecoregions_names) {
  par(mfrow = c(2,2),mar = c(5, 5, 2, 2),oma = c(0, 2, 3, 0))
  for (i in cardamom_var) {
    if (i == 'WOOD') {
      model_benchmark_er_df<-merge_pixel_var_bench_mod_ecoregions(region,i,model_variant,reference[[1]],ecoregions,ecoregions_names)
      #BIOM
      model_rmse<- round(rmse(model_benchmark_er_df$benchmark,model_benchmark_er_df$cci_esa), digits = 0)
      model_lm<- lm(model_benchmark_er_df$benchmark ~ model_benchmark_er_df$cci_esa)
      print(summary(model_lm))
      plot(model_benchmark_er_df$cci_esa, model_benchmark_er_df$benchmark,pch = c(0,1,2,3,4), col = factor(model_benchmark_er_df$ecoregion), 
           main="",xlab=bquote("CARDAMOM C TCWC Biomass ("~Mg~C~ ha^-1~")"), ylab="", 
           xlim=c(0,220), ylim=c(0,220), 
           cex.lab=1.5, cex.axis=1.5)
      abline(coef = c(0,1),col='black', lwd=3)
      text(100, 25, paste('RMSE = ',model_rmse, sep=''),cex=2)
      # text(100, 5, paste('p-value = ',model_lm.pval, sep=''),cex=1.5)
      legend("bottomright",legend = levels(factor(model_benchmark_er_df$ecoregion)),pch = c(0,1,2,3,4),col = factor(levels(factor(model_benchmark_er_df$ecoregion))),cex=1.5)
      clip(x1 = 10, x2 = 200, y1 = 125, y2 = 220)
      abline(lm(model_benchmark_er_df$benchmark ~ model_benchmark_er_df$cci_esa),col='red', lwd=3,xlim=c(20,200))
      
      model_rmse<- round(rmse(model_benchmark_er_df$benchmark,model_benchmark_er_df$rainfor_biomass_annual), digits = 0)
      model_lm=NULL
      model_lm<- lm(model_benchmark_er_df$benchmark ~ model_benchmark_er_df$rainfor_biomass_annual)
      print(summary(model_lm))
      plot(model_benchmark_er_df$rainfor_biomass_annual,model_benchmark_er_df$benchmark,pch = c(0,1,2,3,4),col = factor(model_benchmark_er_df$ecoregion),
           main="",xlab=bquote("CARDAMOM R0 TCWC Biomass ("~Mg~C~ ha^-1~")"), ylab="",
           xlim=c(75,220), ylim=c(75,220), 
           cex.lab=1.5, cex.axis=1.5)
      abline(coef = c(0,1),col='black', lwd=3)
      text(110, 200, paste('RMSE = ',model_rmse, sep=''),cex=2)
      # text(110, 185, paste('p-value = ',model_lm.pval, sep=''),cex=1.5)
      legend("bottomright",legend = levels(factor(model_benchmark_er_df$ecoregion)),pch = c(0,1,2,3,4),col = factor(levels(factor(model_benchmark_er_df$ecoregion))),cex=1.5)
      clip(x1 = 10, x2 = 200, y1 = 125, y2 = 220)
      abline(lm(model_benchmark_er_df$benchmark ~ model_benchmark_er_df$rainfor_biomass_annual),col='red', lwd=3)
      
      model_rmse<- round(rmse(model_benchmark_er_df$benchmark,model_benchmark_er_df$rainfor_2005), digits = 0)
      model_lm=NULL
      model_lm<- lm(model_benchmark_er_df$benchmark ~ model_benchmark_er_df$rainfor_2005)
      print(summary(model_lm))
      plot(model_benchmark_er_df$rainfor_2005,model_benchmark_er_df$benchmark,pch = c(0,1,2,3,4),col = factor(model_benchmark_er_df$ecoregion),
           main="",xlab=bquote("CARDAMOM R1 TCWC Biomass ("~Mg~C~ ha^-1~")"), ylab="",
           xlim=c(75,220), ylim=c(75,220), 
           cex.lab=1.5, cex.axis=1.5)
      abline(coef = c(0,1),col='black', lwd=3)
      text(110, 200, paste('RMSE = ',model_rmse, sep=''),cex=2)
      # text(110, 185, paste('p-value = ',model_lm.pval, sep=''),cex=1.5)
      legend("bottomright",legend = levels(factor(model_benchmark_er_df$ecoregion)),pch = c(0,1,2,3,4),col = factor(levels(factor(model_benchmark_er_df$ecoregion))),cex=1.5)
      clip(x1 = 75, x2 = 220, y1 = 110, y2 = 220)
      abline(lm(model_benchmark_er_df$benchmark ~ model_benchmark_er_df$rainfor_2005),col='red', lwd=3)
      
      model_rmse<- round(rmse(model_benchmark_er_df$benchmark,model_benchmark_er_df$rainfor_annual), digits = 0)
      model_lm=NULL
      model_lm<- lm(model_benchmark_er_df$benchmark ~ model_benchmark_er_df$rainfor_annual)
      print(summary(model_lm))
      plot(model_benchmark_er_df$rainfor_annual,model_benchmark_er_df$benchmark,pch = c(0,1,2,3,4),col = factor(model_benchmark_er_df$ecoregion),
           main="",xlab=bquote("CARDAMOM R2 TCWC Biomass ("~Mg~C~ ha^-1~")"), ylab="",
           xlim=c(75,220), ylim=c(75,220), 
           cex.lab=1.5, cex.axis=1.5)
      abline(coef = c(0,1),col='black', lwd=3)
      text(110, 200, paste('RMSE = ',model_rmse, sep=''),cex=2)
      # text(110, 185, paste('p-value = ',model_lm.pval, sep=''),cex=1.5)
      legend("bottomright",legend = levels(factor(model_benchmark_er_df$ecoregion)),pch = c(0,1,2,3,4), col = factor(levels(factor(model_benchmark_er_df$ecoregion))),cex=1.5)
      clip(x1 = 75, x2 = 220, y1 = 110, y2 = 220)
      abline(lm(model_benchmark_er_df$benchmark ~ model_benchmark_er_df$rainfor_annual),col='red', lwd=3)
      mtext(bquote("TCWC Biomass 2000-2010 ("~ Mg~C~ ha^-1~")"), line=-2, side=2, outer=TRUE, cex=1.5)
      
    }
    else if (i == 'NPP_wood_flx') {
      model_benchmark_er_df<-merge_pixel_var_bench_mod_ecoregions(region,i,model_variant,reference[[2]],ecoregions,ecoregions_names)
      #PROD
      model_rmse<- round(rmse(model_benchmark_er_df$benchmark,model_benchmark_er_df$cci_esa), digits = 2)
      model_lm=NULL
      model_lm<- lm(model_benchmark_er_df$benchmark ~ model_benchmark_er_df$cci_esa)
      print(summary(model_lm))
      plot(model_benchmark_er_df$cci_esa, model_benchmark_er_df$benchmark, col = factor(model_benchmark_er_df$ecoregion),pch = c(0,1,2,3,4), 
           main="",xlab=bquote("CARDAMOM C TCWC Productivity ("~ Mg~C~ ha^-1~year^-1~")"), ylab="", 
           xlim=c(2,9), ylim=c(2,9),
           cex.lab=1.5, cex.axis=1.5)
      abline(coef = c(0,1),col='black', lwd=3)
      text(7, 2, paste('RMSE = ',model_rmse, sep=''),cex=2)
      # text(4, 2, paste('p-value = ',model_lm.pval, sep=''),cex=1.5)
      legend("topleft",legend = levels(factor(model_benchmark_er_df$ecoregion)),pch = c(0,1,2,3,4),col = factor(levels(factor(model_benchmark_er_df$ecoregion))),cex=1.5)
      clip(x1 = 2, x2 = 9, y1 = 2.5, y2 = 4.5)
      abline(lm(model_benchmark_er_df$benchmark ~ model_benchmark_er_df$cci_esa),col='red', lwd=3)
      
      model_rmse<- round(rmse(model_benchmark_er_df$benchmark,model_benchmark_er_df$rainfor_biomass_annual), digits = 2)
      model_lm=NULL
      model_lm<- lm(model_benchmark_er_df$benchmark ~ model_benchmark_er_df$rainfor_biomass_annual)
      print(summary(model_lm))
      plot(model_benchmark_er_df$rainfor_biomass_annual,model_benchmark_er_df$benchmark,pch = c(0,1,2,3,4),col = factor(model_benchmark_er_df$ecoregion),
           main="",xlab=bquote("CARDAMOM R0 TCWC Productivity ("~ Mg~C~ ha^-1~year^-1~")"), ylab="",
           xlim=c(2,9), ylim=c(2,9), 
           cex.lab=1.5, cex.axis=1.5)
      abline(coef = c(0,1),col='black', lwd=3)
      text(7, 2, paste('RMSE = ',model_rmse, sep=''),cex=2)
      # text(4, 2, paste('p-value = ',model_lm.pval, sep=''),cex=1.5)
      legend("topleft",legend = levels(factor(model_benchmark_er_df$ecoregion)),pch = c(0,1,2,3,4),col = factor(levels(factor(model_benchmark_er_df$ecoregion))),cex=1.5)
      clip(x1 = 2, x2 = 9, y1 = 2.5, y2 = 4.5)
      abline(lm(model_benchmark_er_df$benchmark ~ model_benchmark_er_df$rainfor_biomass_annual),col='red', lwd=3)
      
      model_rmse<- round(rmse(model_benchmark_er_df$benchmark,model_benchmark_er_df$rainfor_2005), digits = 2)
      model_lm=NULL
      model_lm<- lm(model_benchmark_er_df$benchmark ~ model_benchmark_er_df$rainfor_2005)
      print(summary(model_lm))
      plot(model_benchmark_er_df$rainfor_2005,model_benchmark_er_df$benchmark,pch = c(0,1,2,3,4),col = factor(model_benchmark_er_df$ecoregion),
           main="",xlab=bquote("CARDAMOM R1 TCWC Productivity ("~ Mg~C~ ha^-1~year^-1~")"), ylab="",
           xlim=c(2,5), ylim=c(2,5), 
           cex.lab=1.5, cex.axis=1.5)
      abline(coef = c(0,1),col='black', lwd=3)
      text(4, 2, paste('RMSE = ',model_rmse, sep=''),cex=2)
      # text(3, 2, paste('p-value = ',model_lm.pval, sep=''),cex=1.5)
      legend("topleft",legend = levels(factor(model_benchmark_er_df$ecoregion)),pch = c(0,1,2,3,4),col = factor(levels(factor(model_benchmark_er_df$ecoregion))),cex=1.5)
      clip(x1 = 2.7, x2 = 5, y1 = 2.7, y2 = 5)
      abline(lm(model_benchmark_er_df$benchmark ~ model_benchmark_er_df$rainfor_2005),col='red', lwd=3)
      
      model_rmse<- round(rmse(model_benchmark_er_df$benchmark,model_benchmark_er_df$rainfor_annual), digits = 2)
      model_lm=NULL
      model_lm<- lm(model_benchmark_er_df$benchmark ~ model_benchmark_er_df$rainfor_annual)
      print(summary(model_lm))
      plot(model_benchmark_er_df$rainfor_annual,model_benchmark_er_df$benchmark,pch = c(0,1,2,3,4),col = factor(model_benchmark_er_df$ecoregion),
           main="",xlab=bquote("CARDAMOM R2 TCWC Productivity ("~ Mg~C~ ha^-1~year^-1~")"), ylab="",
           xlim=c(2,5), ylim=c(2,5), 
           cex.lab=1.5, cex.axis=1.5)
      abline(coef = c(0,1),col='black', lwd=3)
      text(4, 2, paste('RMSE = ',model_rmse, sep=''),cex=2)
      # text(3, 2, paste('p-value = ',model_lm.pval, sep=''),cex=1.5)
      legend("topleft",legend = levels(factor(model_benchmark_er_df$ecoregion)),pch = c(0,1,2,3,4), col = factor(levels(factor(model_benchmark_er_df$ecoregion))),cex=1.5)
      clip(x1 = 2.7, x2 = 5, y1 = 2.7, y2 = 5)
      abline(lm(model_benchmark_er_df$benchmark ~ model_benchmark_er_df$rainfor_annual),col='red', lwd=3)
      mtext(bquote("TCWC Productivity 2000-2010 ("~ Mg~C~ ha^-1~year^-1~")"), line=-2, side=2, outer=TRUE, cex=1.5)
    }
    else if (i == 'OUTPUT_wood_flx') {
      model_benchmark_er_df<-merge_pixel_var_bench_mod_ecoregions(region,i,model_variant,reference[[3]],ecoregions,ecoregions_names)
      #MORT
      model_rmse<- round(rmse(model_benchmark_er_df$benchmark,model_benchmark_er_df$cci_esa), digits = 2)
      model_lm=NULL
      model_lm<- lm(model_benchmark_er_df$benchmark ~ model_benchmark_er_df$cci_esa)
      print(summary(model_lm))
      plot(model_benchmark_er_df$cci_esa, model_benchmark_er_df$benchmark, col = factor(model_benchmark_er_df$ecoregion),pch = c(0,1,2,3,4), 
           main="",xlab=bquote("CARDAMOM C TCWC Mortality ("~ Mg~C~ ha^-1~year^-1~")"), ylab="", 
           xlim=c(2,9), ylim=c(2,9),
           cex.lab=1.5, cex.axis=1.5)
      abline(coef = c(0,1),col='black', lwd=3)
      text(6, 8, paste('RMSE = ',model_rmse, sep=''),cex=2)
      # text(6, 8.5, paste('p-value = ',model_lm.pval, sep=''),cex=1.5)
      legend("topleft",legend = levels(factor(model_benchmark_er_df$ecoregion)),pch = c(0,1,2,3,4),col = factor(levels(factor(model_benchmark_er_df$ecoregion))),cex=1.5)
      clip(x1 = 2, x2 = 8, y1 = 2, y2 = 4)
      abline(lm(model_benchmark_er_df$benchmark ~ model_benchmark_er_df$cci_esa),col='red', lwd=3)
      
      model_rmse<- round(rmse(model_benchmark_er_df$benchmark,model_benchmark_er_df$rainfor_biomass_annual), digits = 2)
      model_lm=NULL
      model_lm<- lm(model_benchmark_er_df$benchmark ~ model_benchmark_er_df$rainfor_biomass_annual)
      print(summary(model_lm))
      plot(model_benchmark_er_df$rainfor_biomass_annual,model_benchmark_er_df$benchmark,pch = c(0,1,2,3,4),col = factor(model_benchmark_er_df$ecoregion),
           main="",xlab=bquote("CARDAMOM R0 TCWC Mortality ("~ Mg~C~ ha^-1~year^-1~")"), ylab="",
           xlim=c(2,9), ylim=c(2,9), 
           cex.lab=1.5, cex.axis=1.5)
      abline(coef = c(0,1),col='black', lwd=3)
      text(6, 8, paste('RMSE = ',model_rmse, sep=''),cex=2)
      # text(6, 8.5, paste('p-value = ',model_lm.pval, sep=''),cex=1.5)
      legend("topleft",legend = levels(factor(model_benchmark_er_df$ecoregion)),pch = c(0,1,2,3,4),col = factor(levels(factor(model_benchmark_er_df$ecoregion))),cex=1.5)
      clip(x1 = 2, x2 = 8, y1 = 2, y2 = 4)
      abline(lm(model_benchmark_er_df$benchmark ~ model_benchmark_er_df$rainfor_biomass_annual),col='red', lwd=3)
      
      model_rmse<- round(rmse(model_benchmark_er_df$benchmark,model_benchmark_er_df$rainfor_2005), digits = 2)
      model_lm=NULL
      lm(model_benchmark_er_df$benchmark ~ model_benchmark_er_df$rainfor_2005)
      print(summary(model_lm))
      plot(model_benchmark_er_df$rainfor_2005,model_benchmark_er_df$benchmark,pch = c(0,1,2,3,4),col = factor(model_benchmark_er_df$ecoregion),
           main="",xlab=bquote("CARDAMOM R1 TCWC Mortality ("~ Mg~C~ ha^-1~year^-1~")"), ylab="",
           xlim=c(2,5), ylim=c(2,5), 
           cex.lab=1.5, cex.axis=1.5)
      abline(coef = c(0,1),col='black', lwd=3)
      text(3.5, 4.5, paste('RMSE = ',model_rmse, sep=''),cex=2)
      # text(3.5, 4.8, paste('p-value = ',model_lm.pval, sep=''),cex=1.5)
      legend("topleft",legend = levels(factor(model_benchmark_er_df$ecoregion)),pch = c(0,1,2,3,4),col = factor(levels(factor(model_benchmark_er_df$ecoregion))),cex=1.5)
      clip(x1 = 2.75, x2 = 5, y1 = 2, y2 = 4)
      abline(lm(model_benchmark_er_df$benchmark ~ model_benchmark_er_df$rainfor_2005),col='red', lwd=3)
      
      model_rmse<- round(rmse(model_benchmark_er_df$benchmark,model_benchmark_er_df$rainfor_annual), digits = 2)
      model_lm=NULL
      model_lm<- lm(model_benchmark_er_df$benchmark ~ model_benchmark_er_df$rainfor_annual)
      print(summary(model_lm))
      plot(model_benchmark_er_df$rainfor_annual,model_benchmark_er_df$benchmark,pch = c(0,1,2,3,4),col = factor(model_benchmark_er_df$ecoregion),
           main="",xlab=bquote("CARDAMOM R2 TCWC Mortality ("~ Mg~C~ ha^-1~year^-1~")"), ylab="",
           xlim=c(2,5), ylim=c(2,5), 
           cex.lab=1.5, cex.axis=1.5)
      abline(coef = c(0,1),col='black', lwd=3)
      text(3.5, 4.5, paste('RMSE = ',model_rmse, sep=''),cex=2)
      # text(3.5, 4.8, paste('p-value = ',model_lm.pval, sep=''),cex=1.5)
      legend("topleft",legend = levels(factor(model_benchmark_er_df$ecoregion)),pch = c(0,1,2,3,4), col = factor(levels(factor(model_benchmark_er_df$ecoregion))),cex=1.5)
      clip(x1 = 2.75, x2 = 5, y1 = 2, y2 = 4)
      abline(lm(model_benchmark_er_df$benchmark ~ model_benchmark_er_df$rainfor_annual),col='red', lwd=3)
      mtext(bquote("TCWC Mortality 2000-2010 ("~ Mg~C~ ha^-1~year^-1~")"), line=-2, side=2, outer=TRUE, cex=1.5)
    }
    else {
      model_benchmark_er_df<-merge_pixel_var_bench_mod_ecoregions(region,i,model_variant,reference[[4]],ecoregions,ecoregions_names)
      #RT
      model_rmse<- round(rmse(model_benchmark_er_df$benchmark,model_benchmark_er_df$cci_esa), digits = 0)
      model_lm<- lm(model_benchmark_er_df$benchmark ~ model_benchmark_er_df$cci_esa)
      print(summary(model_lm))
      plot(model_benchmark_er_df$cci_esa, model_benchmark_er_df$benchmark, col = factor(model_benchmark_er_df$ecoregion),pch = c(0,1,2,3,4), 
           main="",xlab="CARDAMOM C TCWC RT (years)", ylab="", 
           xlim=c(0,80), ylim=c(0,80),
           cex.lab=1.5, cex.axis=1.5)
      abline(coef = c(0,1),col='black', lwd=3)
      text(40, 5, paste('RMSE = ',model_rmse, sep=''),cex=2)
      legend("bottomright",legend = levels(factor(model_benchmark_er_df$ecoregion)),pch = c(0,1,2,3,4),col = factor(levels(factor(model_benchmark_er_df$ecoregion))),cex=1.5)
      clip(x1 = 0, x2 = 30, y1 = 40, y2 = 80)
      abline(lm(model_benchmark_er_df$benchmark ~ model_benchmark_er_df$cci_esa),col='red', lwd=3)
      
      model_rmse<- round(rmse(model_benchmark_er_df$benchmark,model_benchmark_er_df$rainfor_biomass_annual), digits = 0)
      model_lm<- lm(model_benchmark_er_df$benchmark ~ model_benchmark_er_df$rainfor_biomass_annual)
      print(summary(model_lm))
      plot(model_benchmark_er_df$rainfor_biomass_annual,model_benchmark_er_df$benchmark,pch = c(0,1,2,3,4),col = factor(model_benchmark_er_df$ecoregion),
           main="",xlab="CARDAMOM R0 TCWC RT (years)", ylab="",
           xlim=c(0,80), ylim=c(0,80), 
           cex.lab=1.5, cex.axis=1.5)
      abline(coef = c(0,1),col='black', lwd=3)
      text(40, 5, paste('RMSE = ',model_rmse, sep=''),cex=2)
      legend("bottomright",legend = levels(factor(model_benchmark_er_df$ecoregion)),pch = c(0,1,2,3,4),col = factor(levels(factor(model_benchmark_er_df$ecoregion))),cex=1.5)
      clip(x1 = 10, x2 = 60, y1 = 40, y2 = 80)
      abline(lm(model_benchmark_er_df$benchmark ~ model_benchmark_er_df$rainfor_biomass_annual),col='red', lwd=3)
      
      model_rmse<- round(rmse(model_benchmark_er_df$benchmark,model_benchmark_er_df$rainfor_2005), digits = 0)
      lm(model_benchmark_er_df$benchmark ~ model_benchmark_er_df$rainfor_2005)
      print(summary(model_lm))
      plot(model_benchmark_er_df$rainfor_2005,model_benchmark_er_df$benchmark,pch = c(0,1,2,3,4),col = factor(model_benchmark_er_df$ecoregion),
           main="",xlab="CARDAMOM R1 TCWC RT (years)", ylab="",
           xlim=c(0,80), ylim=c(0,80), 
           cex.lab=1.5, cex.axis=1.5)
      abline(coef = c(0,1),col='black', lwd=3)
      text(40, 5, paste('RMSE = ',model_rmse, sep=''),cex=2)
      legend("bottomright",legend = levels(factor(model_benchmark_er_df$ecoregion)),pch = c(0,1,2,3,4),col = factor(levels(factor(model_benchmark_er_df$ecoregion))),cex=1.5)
      clip(x1 = 20, x2 = 60, y1 = 40, y2 = 80)
      abline(lm(model_benchmark_er_df$benchmark ~ model_benchmark_er_df$rainfor_2005),col='red', lwd=3)
      
      model_rmse<- round(rmse(model_benchmark_er_df$benchmark,model_benchmark_er_df$rainfor_annual), digits = 0)
      model_lm<- lm(model_benchmark_er_df$benchmark ~ model_benchmark_er_df$rainfor_annual)
      print(summary(model_lm))
      plot(model_benchmark_er_df$rainfor_annual,model_benchmark_er_df$benchmark,pch = c(0,1,2,3,4),col = factor(model_benchmark_er_df$ecoregion),
           main="",xlab="CARDAMOM R2 TCWC RT (years)", ylab="",
           xlim=c(0,80), ylim=c(0,80), 
           cex.lab=1.5, cex.axis=1.5)
      abline(coef = c(0,1),col='black', lwd=3)
      text(40, 5, paste('RMSE = ',model_rmse, sep=''),cex=2)
      legend("bottomright",legend = levels(factor(model_benchmark_er_df$ecoregion)),pch = c(0,1,2,3,4), col = factor(levels(factor(model_benchmark_er_df$ecoregion))),cex=1.5)
      clip(x1 = 30, x2 = 70, y1 = 40, y2 = 80)
      abline(lm(model_benchmark_er_df$benchmark ~ model_benchmark_er_df$rainfor_annual),col='red', lwd=3)
      mtext(bquote("TCWC Residence Time 2000-2010 (years)"), line=-2, side=2, outer=TRUE, cex=1.5)
    }
  }
}

#run function to plot all 4 CARDAMOM runs against extrapolated RAINFOR amazonian data
new_function_to_plot_them_all(amazonia_ifl_layer,compare_var_names[4],updated_mod_var,reference_data_ifl,reg_data,reg_names)
