##############################################################
################################models########################
##############################################################
#install packages----
library(ncdf4); library(raster); library(dplyr); library(ggplot2);library(ggpubr);library(quantreg);library(ggpp);library(rgeos);library(ggpmisc);library(rgdal)

#done
#Leeds data----
#wood productivity
woodprod_00_09 <- brick('R://brazil_leeds_maps/WoodyProductivity20002009_Mg_perHa_perYear_111km.tif')
woodprod_10_16 <- brick('R://brazil_leeds_maps/WoodyProductivity20102016_Mg_perHa_perYear_111km.tif')
#biomass mortality
biommort_00_09 <- brick('R://brazil_leeds_maps/BiomassMortality_20002009_Mg_perHa_perYear_111km.tif')
biommort_10_16 <- brick('R://brazil_leeds_maps/BiomassMortality_20102016_Mg_perHa_perYear_111km.tif')
#biomass
biomass_amazon <- brick('R://brazil_leeds_maps/AbovegroundBiomass_Mg_perHa_111km.tif')

#convert from Mg_perHa_perYear to g_perM2_perDay
thayr_to_gm2day_fun <- function(x) { x * (100/365.25) }
thayr_to_gCm2day_fun <- function(x) { (x/2) * (100/365.25) }
tha_to_gCm2_fun <- function(x) { (x/2) * 100 }

biomass_amazon_gCm2 <- calc(biomass_amazon, tha_to_gCm2_fun)

woodprod_00_09_gCm2d <- calc(woodprod_00_09, thayr_to_gCm2day_fun)
woodprod_10_16_gCm2d <- calc(woodprod_10_16, thayr_to_gCm2day_fun)
biommort_00_09_gCm2d <- calc(biommort_00_09, thayr_to_gCm2day_fun)
biommort_10_16_gCm2d <- calc(biommort_10_16, thayr_to_gCm2day_fun)

#done
#copy RAINFOR data back for CARDAMOM assimilation----
###
## Wood biomass
# Write back out with new name and units etc
writeRaster(biomass_amazon_gCm2, filename = "G://cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/1deg_converted/wood_biomass/wood_biomass_gCm2.tif", format = "GTiff",overwrite=TRUE)
# Write an estimate of the uncertainty back out
writeRaster(biomass_amazon_gCm2*0.25, filename = "G://cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/1deg_converted/wood_biomass/unc_wood_biomass_gCm2.tif", format = "GTiff",overwrite=TRUE)
writeRaster(biomass_amazon_gCm2*0.05, filename = "G://cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/1deg_converted/wood_biomass/unc5_wood_biomass_gCm2.tif", format = "GTiff",overwrite=TRUE)

no_of_years <- as.numeric(c(2000:2016))
no_of_years_obj <- c(2001:2016)

output_location_a <- "G://cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/1deg_converted/wood_biomass/wood_biomass_gCm2_"
output_location_b <- "G://cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/1deg_converted/wood_biomass/unc_wood_biomass_gCm2_"

save_annually <- function (a,year_no) {
  for (t in year_no) {
    # Write back out with new name and units etc
    writeRaster(a, filename = paste(output_location_a,t, ".tif",sep=""), format = "GTiff",overwrite=TRUE)
    # Write an estimate of the uncertainty back out
    writeRaster(a*0.25, filename = paste(output_location_b,t, ".tif",sep=""), format = "GTiff",overwrite=TRUE)
  }
}
save_annually(biomass_amazon_gCm2,no_of_years)
###
## Wood productivity
## Write back out with new name and units etc
writeRaster(woodprod_00_09_gCm2d, filename = "G://cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/1deg_converted/wood_productivity/wood_productivity_gCm2_2000_2009.tif", format = "GTiff")
## Write an estimate of the uncertainty back out
writeRaster(woodprod_00_09_gCm2d*0.25, filename = "G://cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/1deg_converted/wood_productivity/unc_wood_productivity_gCm2_2000_2009.tif", format = "GTiff",overwrite=TRUE)
writeRaster(woodprod_00_09_gCm2d*0.05, filename = "G://cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/1deg_converted/wood_productivity/unc5_wood_productivity_gCm2_2000_2009.tif", format = "GTiff",overwrite=TRUE)

## Write back out with new name and units etc
writeRaster(woodprod_10_16_gCm2d, filename = "G://cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/1deg_converted/wood_productivity/wood_productivity_gCm2_2010_2016.tif", format = "GTiff")
## Write an estimate of the uncertainty back out
writeRaster(woodprod_10_16_gCm2d*0.25, filename = "G://cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/1deg_converted/wood_productivity/unc_wood_productivity_gCm2_2010_2016.tif", format = "GTiff",overwrite=TRUE)
writeRaster(woodprod_10_16_gCm2d*0.05, filename = "G://cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/1deg_converted/wood_productivity/unc5_wood_productivity_gCm2_2010_2016.tif", format = "GTiff",overwrite=TRUE)

###
## Wood mortality

# # Write back out with new name and units etc
writeRaster(biommort_00_09_gCm2d, filename = "G://cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/1deg_converted/wood_mortality/wood_mortality_gCm2_2000_2009.tif", format = "GTiff")
# # Write an estimate of the uncertainty back out
writeRaster(biommort_00_09_gCm2d*0.25, filename = "G://cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/1deg_converted/wood_mortality/unc_wood_mortality_gCm2_2000_2009.tif", format = "GTiff",overwrite=TRUE)
writeRaster(biommort_00_09_gCm2d*0.05, filename = "G://cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/1deg_converted/wood_mortality/unc5_wood_mortality_gCm2_2000_2009.tif", format = "GTiff",overwrite=TRUE)
#
# # Write back out with new name and units etc
writeRaster(biommort_10_16_gCm2d, filename = "G://cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/1deg_converted/wood_mortality/wood_mortality_gCm2_2010_2016.tif", format = "GTiff")
# # Write an estimate of the uncertainty back out
writeRaster(biommort_10_16_gCm2d*0.25, filename = "G://cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/1deg_converted/wood_mortality/unc_wood_mortality_gCm2_2010_2016.tif", format = "GTiff",overwrite=TRUE)
writeRaster(biommort_10_16_gCm2d*0.05, filename = "G://cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/1deg_converted/wood_mortality/unc5_wood_mortality_gCm2_2010_2016.tif", format = "GTiff",overwrite=TRUE)
#
####subset######
amazonia_subset <- shapefile("./data/amazonia_subset.shp")

extract_subset <- function (region,reference){
  masked1 <- mask(reference, region)
  masked2 <- masked1 > -Inf
  poly <- rasterToPolygons(masked2, dissolve=TRUE)
  data_region <- mask(reference,poly)
  return(data_region)
}

#amazonia_subset_raster<-extract_subset(amazonia_subset,biomass_amazon_gCm2)
#plot(amazonia_subset_raster)
#writeRaster(amazonia_subset_raster, filename = "G://cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/amazonia_subset.tif", format = "GTiff")

#done
