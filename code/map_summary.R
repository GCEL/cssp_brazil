#Summary of Leeds data----
#download packages
library(sp)
library(rgdal)
library(raster)
library(dplyr)
library(ggplot2)
library(viridis)
library(rasterVis)
library(ncdf4)

#upload data
#wood productivity
woodprod_00_09 <- brick('R://brazil_leeds_maps/WoodyProductivity20002009_Mg_perHa_perYear_111km.tif')
woodprod_10_16 <- brick('R://brazil_leeds_maps/WoodyProductivity20102016_Mg_perHa_perYear_111km.tif')
#biomass mortality
biommort_00_09 <- brick('R://brazil_leeds_maps/BiomassMortality_20002009_Mg_perHa_perYear_111km.tif')
biommort_10_16 <- brick('R://brazil_leeds_maps/BiomassMortality_20102016_Mg_perHa_perYear_111km.tif')
#biomass
biomass_amazon <- brick('R://brazil_leeds_maps/AbovegroundBiomass_Mg_perHa_111km.tif')

#summary
summary(woodprod_00_09)
str(woodprod_00_09)

#plots
par(mfrow=c(2,2))
plot(woodprod_00_09,main='Woody Productivity 2000-2009')
plot(woodprod_10_16,main='Woody Productivity 2010-2016')
plot(biommort_00_09,main='Biomass Mortality 2000-2009')
plot(biommort_10_16,main='Biomass Mortality 2010-2016')
par(mfrow=c(1,1))
plot(biomass_amazon,main='Aboveground Biomass')

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

ncell(woodprod_00_09)
table(is.na(woodprod_10_16[]))

raster_data <- c("woodprod_00_09","woodprod_10_16","biommort_00_09","biommort_10_16")
cell_n <- c(740,740,740,740)
cell_mean <- c(cellStats(woodprod_00_09, 'mean'),cellStats(woodprod_10_16, 'mean'),cellStats(biommort_00_09, 'mean'),cellStats(biommort_10_16, 'mean'))
cell_sd <- c(cellStats(woodprod_00_09, 'sd'),cellStats(woodprod_10_16, 'sd'),cellStats(biommort_00_09, 'sd'),cellStats(biommort_10_16, 'sd'))

stats_leeds_table <- data.frame(raster_data,cell_n,cell_mean,cell_sd)
str(stats_leeds_table)

stats_leeds_table_error <- stats_leeds_table %>%
  mutate(cell_se =cell_sd/sqrt(cell_n))

ggplot(stats_leeds_table_error) +
  geom_bar( aes(x=raster_data, y=cell_mean), stat="identity", fill="black", alpha=0.5) +
  geom_errorbar( aes(x=raster_data, ymin=cell_mean-cell_se, ymax=cell_mean+cell_se), width=0.4, colour="red", alpha=0.9, size=1.5)+ 
  ylim(0,8)+
  labs(x="Data", y = "Productivity (Mg/ ha/ year)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(text = element_text(size = 20))


#models
nc_card <- nc_open('G://CARDAMOM_ILAMB/CARDAMOM_Brazil_1x1_2001_2017_v1.0.nc')
print(nc_card)
nc_close(nc_card)

cardamom_nppwood <- stack('G://CARDAMOM_ILAMB/CARDAMOM_Brazil_1x1_2001_2017_v1.0.nc',varname="NPP_wood_flx")
class(getZ(cardamom_nppwood))

cardamom_nppwood_01_10 <- cardamom_nppwood[[which(getZ(cardamom_nppwood) >= as.Date("2001-01-01") & getZ(cardamom_nppwood) <= as.Date("2009-12-01"))]]
cardamom_nppwood_10_16 <- cardamom_nppwood[[which(getZ(cardamom_nppwood) >= as.Date("2010-01-01") & getZ(cardamom_nppwood) <= as.Date("2016-12-01"))]]

cardamom_nppwood_01_10_indices <- as.numeric(format(as.Date(names(cardamom_nppwood_01_10), format = "X%Y.%m.%d"), format = "%Y"))
cardamom_nppwood_01_10_mean <- stackApply(cardamom_nppwood_01_10, cardamom_nppwood_01_10_indices, fun = mean)

amazon_mask <- brick('./data/ilamb-mask-AMAZON.nc')
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
cellStats(cardamom_nppwood_01_10_mean, 'mean');cellStats(cardamom_nppwood_01_10_mean, 'sd')
