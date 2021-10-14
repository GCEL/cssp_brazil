#Summary of Leeds data----
#download packages
library(sp)
library(rgdal)
library(raster)
library(ggplot2)
library(viridis)
library(rasterVis)
library(ncdf4)

#upload data
#wood productivity
woodprod_00_09 <- brick('./maps/WoodyProductivity20002009_Mg_perHa_perYear_111km.tif')
woodprod_10_16 <- brick('./maps/WoodyProductivity20102016_Mg_perHa_perYear_111km.tif')
#biomass mortality
biommort_00_09 <- brick('./maps/BiomassMortality_20002009_Mg_perHa_perYear_111km.tif')
biommort_10_16 <- brick('./maps/BiomassMortality_20102016_Mg_perHa_perYear_111km.tif')
#biomass
biomass_amazon <- brick('./maps/AbovegroundBiomass_Mg_perHa_111km.tif')

#summary
summary(woodprod_00_09)
str(woodprod_00_09)

#plots
plot(woodprod_00_09)
plot(woodprod_10_16)
plot(biommort_00_09)
plot(biommort_10_16)
plot(biomass_amazon)

hist(woodprod_00_09)
hist(woodprod_10_16)
hist(biommort_00_09)
hist(biommort_10_16)
hist(biomass_amazon)

#mean(woodprod_00_09)
cellStats(woodprod_00_09, 'mean');cellStats(woodprod_00_09, 'sd')
cellStats(woodprod_10_16, 'mean');cellStats(woodprod_10_16, 'sd')
cellStats(biommort_00_09, 'mean');cellStats(biommort_00_09, 'sd')
cellStats(biommort_10_16, 'mean');cellStats(biommort_10_16, 'sd')
cellStats(biomass_amazon, 'mean');cellStats(biomass_amazon, 'sd')

#models
cardamomnpp <- stack('./maps/CARDAMOM.nc', varname="NPP")
latest_cardamom_nppwood <- stack('./maps/CARDAMOM_Brazil_1x1_2001_2017_v1.0.nc',varname="NPP_wood_flx")

# nc_data <- nc_open('./maps/CARDAMOM_Brazil_1x1_2001_2017_v1.0.nc')
# print(nc_data)
# nc_close(nc_data)

amazon_mask <- brick('./maps/ilamb-mask-AMAZON.nc')
amazon_mask_pol <- rasterToPolygons(amazon_mask, dissolve = TRUE)
cardamomnpp_amazon <- crop(cardamomnpp, amazon_mask_pol)
biomass_amazon_extent <- extent(biomass_amazon)
plot(biomass_amazon_extent)

r <- biomass_amazon > -Inf
pp <- rasterToPolygons(r, dissolve=TRUE)
plot(pp)

cardamomnpp_amazon <- crop(cardamomnpp, pp)

summary(amazon_mask)
plot(cardamomnpp$X2001.01.01)
plot(amazon_mask)
plot(amazon_mask_pol)
plot(cardamomnpp_amazon$X2001.01.01)

names(cardamomnpp_amazon[,1:120,])

cellStats(cardamomnpp_amazon, 'mean');cellStats(cardamomnpp_amazon, 'sd')
cellStats(latest_cardamom_nppwood, 'mean');cellStats(latest_cardamom_nppwood, 'sd')
