######################################################
###################Summary of Leeds data##############
######################################################
#download packages----
library(sp)
library(rgdal)
library(raster)
library(dplyr)
library(ggplot2)
library(viridis)
library(rasterVis)
library(ncdf4)

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

woodprod_00_09_gm2d <- calc(woodprod_00_09, thayr_to_gCm2day_fun)
woodprod_10_16_gm2d <- calc(woodprod_10_16, thayr_to_gCm2day_fun)
biommort_00_09_gm2d <- calc(biommort_00_09, thayr_to_gCm2day_fun)
biommort_10_16_gm2d <- calc(biommort_10_16, thayr_to_gCm2day_fun)

biomass_amazon_gCm2 <- calc(biomass_amazon, tha_to_gCm2_fun)
biomass_amazon_gCm2_card <- brick('G://cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/1deg/wood_biomass/wood_biomass_gCm2.tif')
npp_amazon_gCm2_card <- brick('G://cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/1deg/wood_productivity/wood_productivity_gCm2_2000_2009.tif')

#summary
summary(woodprod_00_09)
str(woodprod_00_09)

#plots
par(mfrow=c(1,2))
plot(woodprod_00_09,main='Woody Productivity 2000-2009 Mg/Ha/Year')
plot(woodprod_10_16,main='Woody Productivity 2010-2016 Mg/Ha/Year')
plot(biommort_00_09,main='Biomass Mortality 2000-2009 Mg/Ha/Year')
plot(biommort_10_16,main='Biomass Mortality 2010-2016 Mg/Ha/Year')

par(mfrow=c(1,1))
plot(biomass_amazon,main='Aboveground Biomass Mg/Ha')

par(mfrow=c(2,2))
plot(woodprod_00_09_gm2d,main='Woody Productivity 2000-2009 g/m2/day')
plot(woodprod_10_16_gm2d,main='Woody Productivity 2010-2016 g/m2/day')
plot(biommort_00_09_gm2d,main='Biomass Mortality 2000-2009 g/m2/day')
plot(biommort_10_16_gm2d,main='Biomass Mortality 2010-2016 g/m2/day')

hist(woodprod_00_09)
hist(woodprod_10_16)
hist(biommort_00_09)
hist(biommort_10_16)
hist(biomass_amazon)

#mean(woodprod_00_09)
cellStats(woodprod_00_09_gm2d, 'mean');cellStats(woodprod_00_09_gm2d, 'sd')
cellStats(woodprod_10_16_gm2d, 'mean');cellStats(woodprod_10_16_gm2d, 'sd')
cellStats(biommort_00_09_gm2d, 'mean');cellStats(biommort_00_09_gm2d, 'sd')
cellStats(biommort_10_16_gm2d, 'mean');cellStats(biommort_10_16_gm2d, 'sd')
cellStats(biomass_amazon_gCm2, 'mean');cellStats(biomass_amazon_gCm2, 'sd')

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

#Cwood ESA source AGB----
esa_cci_agb_17 <- brick('G://AGB/ESA_CCI_BIOMASS/ESA_CCI_AGB_1deg/AGB_map_MgCha_2017.tif')
plot(esa_cci_agb_17)
esa_cci_agb_17_amazon <- crop(esa_cci_agb_17, extent(cardamom_benchmarknppwood_1101))
plot(esa_cci_agb_17_amazon)

check_extent <- extent(esa_cci_agb_17_amazon) #use south america extent from INLAND SA
check_extent_r <- raster(check_extent)
res(check_extent_r) <- res(esa_cci_agb_17_amazon)
values(check_extent_r) <- 1
crs(check_extent_r) <- "+proj=longlat +datum=WGS84 +no_defs"
check_extent_r_amazonia <- mask(check_extent_r, amazonia_poly, updatevalue=0)
check_extent_r_amazonia[check_extent_r_amazonia<1] <- NA

esa_cci_agb_17_amazon_masked <- mask(esa_cci_agb_17_amazon, check_extent_r_amazonia)

esa_cci_agb_17_amazon_masked <- tha_to_gm2_fun(esa_cci_agb_17_amazon_masked)



#done
#models----
#CARDAMOM
cardamom_benchmarknppwood <- stack('R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/MODELS/CARDAMOM_NORAINFOR/NoRainfor_Amazon_1deg_monthly_2001_2019.nc',varname="WOOD")
#class(getZ(cardamom_benchmarknppwood))
#plot(cardamom_benchmarknppwood$X2001.01.01)
#names(cardamom_benchmarknppwood)
cardamom_benchmarknppwood_masked_01_16_mean <- stackApply(cardamom_benchmarknppwood, indices =  rep(1,nlayers(cardamom_benchmarknppwood)), fun = "mean")
cellStats(cardamom_benchmarknppwood_masked_01_16_mean, 'mean');cellStats(cardamom_benchmarknppwood_masked_01_16_mean, 'sd')

cardamom_benchmarknppwood_1101<-cardamom_benchmarknppwood$X2001.01.01 #template to get data ignore values (green in plot of previous step)
# cardamom_nppwood_1101[cardamom_nppwood_1101 > 500000] <- NA #ignore values above 500000
# plot(cardamom_nppwood_1101) #view masked out values

#cardamom_nppwood_masked <- mask(cardamom_nppwood, cardamom_nppwood_1101) #mask out data ignore values
#plot(cardamom_nppwood_masked$X2001.01.01) #check successful masking

#extract raster images of two periods of interest
cardamom_benchmarknppwood_masked_01_10 <- cardamom_benchmarknppwood[[which(getZ(cardamom_benchmarknppwood) >= as.Date("2001-01-01") & getZ(cardamom_benchmarknppwood) <= as.Date("2009-12-01"))]]
cardamom_benchmarknppwood_masked_10_16 <- cardamom_benchmarknppwood[[which(getZ(cardamom_benchmarknppwood) >= as.Date("2010-01-01") & getZ(cardamom_benchmarknppwood) <= as.Date("2016-12-01"))]]

#calculate mean over time
cardamom_benchmarknppwood_masked_01_10_mean <- stackApply(cardamom_benchmarknppwood_masked_01_10, indices =  rep(1,nlayers(cardamom_benchmarknppwood_masked_01_10)), fun = "mean")
cardamom_benchmarknppwood_masked_10_16_mean <- stackApply(cardamom_benchmarknppwood_masked_10_16, indices =  rep(1,nlayers(cardamom_benchmarknppwood_masked_10_16)), fun = "mean")
#plot(cardamom_benchmarknppwood_masked_01_10_mean)
#plot(cardamom_benchmarknppwood_masked_10_16_mean)

#visualize mean NPP flux over two periods
#par(mfrow=c(1,2))
#plot(cardamom_benchmarknppwood_masked_01_10_mean,main='CARDAMOM NPP flux for wood 2001-2009 g/m2/day')
#plot(cardamom_benchmarknppwood_masked_10_16_mean,main='CARDAMOM NPP flux for wood 2010-2016')

#calculate spatial statistics (mean, sd) of two periods 
cellStats(cardamom_benchmarknppwood_masked_01_10_mean, 'mean');cellStats(cardamom_benchmarknppwood_masked_01_10_mean, 'sd')
cellStats(cardamom_benchmarknppwood_masked_10_16_mean, 'mean');cellStats(cardamom_benchmarknppwood_masked_10_16_mean, 'sd')

#extract NPP flux variable and convert to rasterer
cardamom_nppwood <- stack('R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/MODELS/CARDAMOM_NORAINFOR/NoRainfor_Amazon_1deg_monthly_2001_2019.nc',varname="OUTPUT_wood_flx")
class(getZ(cardamom_nppwood))
plot(cardamom_nppwood$X2001.01.01)
names(cardamom_nppwood)

cardamom_nppwood_1101<-cardamom_nppwood$X2001.01.01 #template to get data ignore values (green in plot of previous step)
cardamom_nppwood_1101[cardamom_nppwood_1101 > 500000] <- NA #ignore values above 500000
plot(cardamom_nppwood_1101) #view masked out values

cardamom_nppwood_masked <- mask(cardamom_nppwood, cardamom_nppwood_1101) #mask out data ignore values
plot(cardamom_nppwood_masked$X2001.01.01) #check successful masking

#extract raster images of two periods of interest
cardamom_nppwood_masked_01_10 <- cardamom_nppwood[[which(getZ(cardamom_nppwood) >= as.Date("2001-01-01") & getZ(cardamom_nppwood) <= as.Date("2009-12-01"))]]
cardamom_nppwood_masked_10_16 <- cardamom_nppwood[[which(getZ(cardamom_nppwood) >= as.Date("2010-01-01") & getZ(cardamom_nppwood) <= as.Date("2016-12-01"))]]

#calculate mean over time
cardamom_nppwood_masked_01_10_mean <- stackApply(cardamom_nppwood_masked_01_10, indices =  rep(1,nlayers(cardamom_nppwood_masked_01_10)), fun = "mean")
cardamom_nppwood_masked_10_16_mean <- stackApply(cardamom_nppwood_masked_10_16, indices =  rep(1,nlayers(cardamom_nppwood_masked_10_16)), fun = "mean")
plot(cardamom_nppwood_masked_01_10_mean)
plot(cardamom_nppwood_masked_10_16_mean)

#visualize mean NPP flux over two periods
par(mfrow=c(1,2))
plot(cardamom_nppwood_masked_01_10_mean,main='CARDAMOM NPP flux for wood 2001-2009 g/m2/day')
plot(cardamom_nppwood_masked_10_16_mean,main='CARDAMOM NPP flux for wood 2010-2016')

#calculate spatial statistics (mean, sd) of two periods 
cellStats(cardamom_nppwood_masked_01_10_mean, 'mean');cellStats(cardamom_nppwood_masked_01_10_mean, 'sd')
cellStats(cardamom_nppwood_masked_10_16_mean, 'mean');cellStats(cardamom_nppwood_masked_10_16_mean, 'sd')

# cardamom_nppwood_01_10 <- cardamom_nppwood[[which(getZ(cardamom_nppwood) >= as.Date("2001-01-01") & getZ(cardamom_nppwood) <= as.Date("2009-12-01"))]]
# cardamom_nppwood_10_16 <- cardamom_nppwood[[which(getZ(cardamom_nppwood) >= as.Date("2010-01-01") & getZ(cardamom_nppwood) <= as.Date("2016-12-01"))]]
# 
# cardamom_nppwood_01_10_indices <- as.numeric(format(as.Date(names(cardamom_nppwood_01_10), format = "X%Y.%m.%d"), format = "%Y"))
# cardamom_nppwood_01_10_mean <- stackApply(cardamom_nppwood_01_10, cardamom_nppwood_01_10_indices, fun = mean)
# plot(cardamom_nppwood_01_10_mean$index_2001)

#JULES
jules_nppwood <- stack('R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/MODELS/JULES/JULES_monthly_1x1_Amazon_combined_2001_2019.nc',varname="OUTPUT_wood_flx")
class(getZ(jules_nppwood))
plot(jules_nppwood$X2001.01.01)
names(jules_nppwood)

jules_nppwood_1101<-jules_nppwood$X2001.01.01 #template to get data ignore values (green in plot of previous step)
jules_nppwood_1101[jules_nppwood_1101 > 500000] <- NA #ignore values above 500000
plot(jules_nppwood_1101) #view masked out values

extent(jules_nppwood)<-extent(cardamom_nppwood_1101)
jules_nppwood_masked <- mask(jules_nppwood, cardamom_nppwood_1101) #mask out data ignore values
plot(jules_nppwood_masked$X2001.01.01) #check successful masking

#extract raster images of two periods of interest
jules_nppwood_masked_01_10 <- jules_nppwood_masked[[which(getZ(jules_nppwood_masked) >= as.Date("2001-01-01") & getZ(jules_nppwood_masked) <= as.Date("2009-12-01"))]]
jules_nppwood_masked_10_16 <- jules_nppwood_masked[[which(getZ(jules_nppwood_masked) >= as.Date("2010-01-01") & getZ(jules_nppwood_masked) <= as.Date("2016-12-01"))]]

#calculate mean over time
jules_nppwood_masked_01_10_mean <- stackApply(jules_nppwood_masked_01_10, indices =  rep(1,nlayers(jules_nppwood_masked_01_10)), fun = "mean")
jules_nppwood_masked_10_16_mean <- stackApply(jules_nppwood_masked_10_16, indices =  rep(1,nlayers(jules_nppwood_masked_10_16)), fun = "mean")
plot(jules_nppwood_masked_01_10_mean)
plot(jules_nppwood_masked_10_16_mean)

#visualize mean NPP flux over two periods
par(mfrow=c(1,2))
plot(jules_nppwood_masked_01_10_mean,main='jules NPP flux for wood 2001-2009 g/m2/day')
plot(jules_nppwood_masked_10_16_mean,main='jules NPP flux for wood 2010-2016')

#calculate spatial statistics (mean, sd) of two periods 
cellStats(jules_nppwood_masked_01_10_mean, 'mean');cellStats(jules_nppwood_masked_01_10_mean, 'sd')
cellStats(jules_nppwood_masked_10_16_mean, 'mean');cellStats(jules_nppwood_masked_10_16_mean, 'sd')

#extract MTT wood variable and convert to rasterer
cardamom_mttwood <- stack('G://CARDAMOM_ILAMB/CARDAMOM_Brazil_1x1_2001_2017_v1.0.nc',varname="MTT_wood")
class(getZ(cardamom_mttwood))
plot(cardamom_mttwood$X2001.01.01)
cardamom_mttwood_indices <- as.numeric(format(as.Date(names(cardamom_mttwood), format = "X%Y.%m.%d"), format = "%Y"))

cardamom_mttwood_1101<-cardamom_mttwood$X2001.01.01 #template to get data ignore values (green in plot of previous step)
cardamom_mttwood_1101[cardamom_mttwood_1101 > 500000] <- NA #ignore values above 500000
plot(cardamom_mttwood_1101) #view masked out values

cardamom_mttwood_masked <- mask(cardamom_mttwood, cardamom_mttwood_1101) #mask out data ignore values
plot(cardamom_mttwood_masked$X2001.01.01) #check successful masking

#extract MTT wood variable and convert to rasterer
cardamom_cwood <- stack('G://CARDAMOM_ILAMB/CARDAMOM_Brazil_1x1_2001_2017_v1.0.nc',varname="WOOD")
class(getZ(cardamom_cwood))
plot(cardamom_cwood$X2001.01.01)

cardamom_cwood_1101<-cardamom_cwood$X2001.01.01 #template to get data ignore values (green in plot of previous step)
cardamom_cwood_1101[cardamom_cwood_1101 > 500000] <- NA #ignore values above 500000
plot(cardamom_cwood_1101) #view masked out values

cardamom_cwood_masked <- mask(cardamom_cwood, cardamom_cwood_1101) #mask out data ignore values
plot(cardamom_cwood_masked$X2001.01.01) #check successful masking

cardamom_outputwood <- stack(cardamom_cwood_masked/cardamom_mttwood_masked)
names(cardamom_outputwood)
getZ(cardamom_outputwood)

cardamom_outputwood <- overlay(cardamom_cwood_masked, cardamom_mttwood_masked, fun=function(x, y) { sqrt(x/y) } )
names(cardamom_outputwood) <- names(cardamom_nppwood)
cardamom_outputwood <- setZ(cardamom_outputwood, getZ(cardamom_nppwood), name='Date')
plot(cardamom_outputwood$X2001.01.01)

gm2yr_to_gm2day_fun <- function(x) { x /365.25 }
cardamom_outputwood_gm2d <- stack(calc(cardamom_outputwood, gm2yr_to_gm2day_fun))
names(cardamom_outputwood_gm2d) <- names(cardamom_nppwood)
cardamom_outputwood_gm2d <- setZ(cardamom_outputwood_gm2d, getZ(cardamom_nppwood), name='Date')
plot(cardamom_outputwood_gm2d$X2001.01.01)

#extract raster images of two periods of interest
cardamom_outputwood_gm2d_masked_01_10 <- cardamom_outputwood_gm2d[[which(getZ(cardamom_outputwood_gm2d) >= as.Date("2001-01-01") & getZ(cardamom_outputwood_gm2d) <= as.Date("2009-12-01"))]]
cardamom_outputwood_gm2d_masked_10_16 <- cardamom_outputwood_gm2d[[which(getZ(cardamom_outputwood_gm2d) >= as.Date("2010-01-01") & getZ(cardamom_outputwood_gm2d) <= as.Date("2016-12-01"))]]

#calculate mean over time
cardamom_outputwood_gm2d_masked_01_10_mean <- stackApply(cardamom_outputwood_gm2d_masked_01_10, indices =  rep(1,nlayers(cardamom_outputwood_gm2d_masked_01_10)), fun = "mean")
cardamom_outputwood_gm2d_masked_10_16_mean <- stackApply(cardamom_outputwood_gm2d_masked_10_16, indices =  rep(1,nlayers(cardamom_outputwood_gm2d_masked_10_16)), fun = "mean")

#visualize mean NPP flux over two periods
par(mfrow=c(1,2))
plot(cardamom_outputwood_gm2d_masked_01_10_mean,main='CARDAMOM Output wood 2000-2009 g/m2/day')
plot(cardamom_outputwood_gm2d_masked_10_16_mean,main='CARDAMOM Output wood 2010-2016 g/m2/day')

#done
#intersect of Brazil and amazon----
nc_amazonmask <- nc_open('R://ILAMB_tutorial/ilamb-mask-AMAZON.nc')
print(nc_amazonmask) #print properties
nc_close(nc_amazonmask) #close nc file

amazon_mask <- raster('R://ILAMB_tutorial/ilamb-mask-AMAZON.nc')
amazon_mask
summary(amazon_mask)
print(amazon_mask)
plot(amazon_mask)

# write the raster layer (tmpin)
outfile <- "amazon_mask_raster.nc"
#crs(amazon_mask) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
writeRaster(amazon_mask, outfile, overwrite=TRUE, format="CDF", varname="ids", varunit="", xname="lon", yname="lat", zname='ids')

nc_amazon_mask <- nc_open('./amazon_mask_raster.nc')
nc_amazon_mask
print(nc_amazon_mask) #print properties
nc_close(nc_amazon_mask) #close nc file

amazon_mask_new <- raster('./amazon_mask_raster.nc')
amazon_mask_new
summary(amazon_mask_new)
print(amazon_mask_new)
plot(amazon_mask_new)

amazon_mask[amazon_mask < 0] <- NA
amazon_mask <- amazon_mask> -Inf
plot(amazon_mask)
amazon_mask_pol <- rasterToPolygons(amazon_mask, dissolve = TRUE)
plot(amazon_mask_pol,add=T)

biomass_amazon_mask <- biomass_amazon > -Inf
biomass_amazon_mask_pol <- rasterToPolygons(biomass_amazon_mask, dissolve=TRUE)
plot(biomass_amazon_mask_pol)
plot(biomass_amazon_mask_pol,add=T)

# cardamom_brazil_mask <- cardamom_nppwood_1101> -Inf
# cardamom_brazil_mask_pol <- rasterToPolygons(cardamom_brazil_mask, dissolve = TRUE)
# plot(cardamom_brazil_mask_pol)
# 
# cardamomnpp_amazon <- crop(cardamomnpp, amazon_mask_pol)
# biomass_amazon_extent <- extent(biomass_amazon)
# plot(biomass_amazon_extent)
# 
# cardamomnpp_amazon <- crop(cardamomnpp, pp)
# 
# summary(amazon_mask)
# plot(cardamomnpp$X2001.01.01)
# plot(amazon_mask)
# plot(amazon_mask_pol)
# plot(cardamomnpp_amazon$X2001.01.01)
# 
# names(cardamomnpp_amazon[,1:120,])
# 
# brazil_amazon_intersect <- intersect(biomass_amazon,cardamom_nppwood_1101)
# brazil_amazon_intersect <- brazil_amazon_intersect> -Inf
# brazil_amazon_intersect_pol <- rasterToPolygons(brazil_amazon_intersect, dissolve = TRUE)
# plot(brazil_amazon_intersect_pol)
# plot(brazil_amazon_intersect)

#crop and mask----
leeds_agb_amazon_brazil_crop <- crop(biomass_amazon, extent(amazon_mask_pol))
leeds_agb_amazon_brazil <- mask(leeds_agb_amazon_brazil_crop, amazon_mask_pol)
par(mfrow=c(1,1))
plot(leeds_agb_amazon_brazil,main='Leeds Brazilian Amazon Aboveground Biomass Mg/ha')

leeds_woodprod_00_09_gm2d_amazon_brazil_crop <- crop(woodprod_00_09_gm2d, extent(amazon_mask_pol))
leeds_woodprod_00_09_gm2d_amazon_brazil <- mask(leeds_woodprod_00_09_gm2d_amazon_brazil_crop, amazon_mask_pol)

leeds_woodprod_10_16_gm2d_amazon_brazil_crop <- crop(woodprod_10_16_gm2d, extent(amazon_mask_pol))
leeds_woodprod_10_16_gm2d_amazon_brazil <- mask(leeds_woodprod_10_16_gm2d_amazon_brazil_crop, amazon_mask_pol)

leeds_biommort_00_09_gm2d_amazon_brazil_crop <- crop(biommort_00_09_gm2d, extent(amazon_mask_pol))
leeds_biommort_00_09_gm2d_amazon_brazil <- mask(leeds_biommort_00_09_gm2d_amazon_brazil_crop, amazon_mask_pol)

leeds_biommort_10_16_gm2d_amazon_brazil_crop <- crop(biommort_10_16_gm2d, extent(amazon_mask_pol))
leeds_biommort_10_16_gm2d_amazon_brazil <- mask(leeds_biommort_10_16_gm2d_amazon_brazil_crop, amazon_mask_pol)

par(mfrow=c(2,2))
plot(leeds_woodprod_00_09_gm2d_amazon_brazil,main='Leeds Brazilian Amazon Woody Productivity 2000-2009 g/m2/day')
plot(leeds_woodprod_10_16_gm2d_amazon_brazil,main='Leeds Brazilian Amazon Woody Productivity 2010-2016 g/m2/day')
plot(leeds_biommort_00_09_gm2d_amazon_brazil,main='Leeds Brazilian Amazon Biomass Mortality 2000-2009 g/m2/day')
plot(leeds_biommort_10_16_gm2d_amazon_brazil,main='Leeds Brazilian Amazon Biomass Mortality 2010-2016 g/m2/day')

hist(leeds_woodprod_00_09_gm2d_amazon_brazil)
hist(leeds_woodprod_10_16_gm2d_amazon_brazil)
hist(leeds_biommort_00_09_gm2d_amazon_brazil)
hist(leeds_biommort_10_16_gm2d_amazon_brazil)

cellStats(leeds_woodprod_00_09_gm2d_amazon_brazil, 'mean');cellStats(leeds_woodprod_00_09_gm2d_amazon_brazil, 'sd')
cellStats(leeds_woodprod_10_16_gm2d_amazon_brazil, 'mean');cellStats(leeds_woodprod_10_16_gm2d_amazon_brazil, 'sd')
cellStats(leeds_biommort_00_09_gm2d_amazon_brazil, 'mean');cellStats(leeds_biommort_00_09_gm2d_amazon_brazil, 'sd')
cellStats(leeds_biommort_10_16_gm2d_amazon_brazil, 'mean');cellStats(leeds_biommort_10_16_gm2d_amazon_brazil, 'sd')

#visualize mean NPP flux over two periods
cardamom_nppwood_01_09_gm2d_amazon_brazil_crop <- crop(cardamom_nppwood_masked_01_10_mean, extent(amazon_mask_pol))
cardamom_nppwood_01_09_gm2d_amazon_brazil <- mask(cardamom_nppwood_01_09_gm2d_amazon_brazil_crop, amazon_mask_pol)
cardamom_outputwood_01_09_gm2d_amazon_brazil_crop <- crop(cardamom_outputwood_gm2d_masked_01_10_mean, extent(amazon_mask_pol))
cardamom_outputwood_01_09_gm2d_amazon_brazil <- mask(cardamom_outputwood_01_09_gm2d_amazon_brazil_crop, amazon_mask_pol)

cardamom_nppwood_10_16_gm2d_amazon_brazil_crop <- crop(cardamom_nppwood_masked_10_16_mean, extent(amazon_mask_pol))
cardamom_nppwood_10_16_gm2d_amazon_brazil <- mask(cardamom_nppwood_10_16_gm2d_amazon_brazil_crop, amazon_mask_pol)
cardamom_outputwood_10_16_gm2d_amazon_brazil_crop <- crop(cardamom_outputwood_gm2d_masked_10_16_mean, extent(amazon_mask_pol))
cardamom_outputwood_10_16_gm2d_amazon_brazil <- mask(cardamom_outputwood_10_16_gm2d_amazon_brazil_crop, amazon_mask_pol)

par(mfrow=c(2,2))
plot(cardamom_nppwood_01_09_gm2d_amazon_brazil,main='Brazilian Amazon CARDAMOM NPP flux for wood 2001-2009 g/m2/day')
plot(cardamom_nppwood_10_16_gm2d_amazon_brazil,main='Brazilian Amazon CARDAMOM NPP flux for wood 2010-2016 g/m2/day')
plot(cardamom_outputwood_01_09_gm2d_amazon_brazil,main='Brazilian Amazon CARDAMOM Output wood  2001-2009 g/m2/day')
plot(cardamom_outputwood_10_16_gm2d_amazon_brazil,main='Brazilian Amazon CARDAMOM Output wood 2010-2016 g/m2/day')

hist(cardamom_nppwood_01_09_gm2d_amazon_brazil)
hist(cardamom_nppwood_10_16_gm2d_amazon_brazil)
hist(cardamom_outputwood_01_09_gm2d_amazon_brazil)
hist(cardamom_outputwood_10_16_gm2d_amazon_brazil)

cellStats(cardamom_nppwood_01_09_gm2d_amazon_brazil, 'mean');cellStats(cardamom_nppwood_01_09_gm2d_amazon_brazil, 'sd')
cellStats(cardamom_nppwood_10_16_gm2d_amazon_brazil, 'mean');cellStats(cardamom_nppwood_10_16_gm2d_amazon_brazil, 'sd')
cellStats(cardamom_outputwood_01_09_gm2d_amazon_brazil, 'mean');cellStats(cardamom_outputwood_01_09_gm2d_amazon_brazil, 'sd')
cellStats(cardamom_outputwood_10_16_gm2d_amazon_brazil, 'mean');cellStats(cardamom_outputwood_10_16_gm2d_amazon_brazil, 'sd')


#INLAND----
nc_inland <- nc_open('G://ILAMB_runs_output/CSSP_stippling/INLAND/INLAND_monthly_1x1_SAmerica.nc')
print(nc_inland) #print properties
nc_close(nc_inland) #close nc file

#G://ILAMB_runs_output/CSSP_stippling/DATA/benchmark/CARDAMOM_monthly_1x1_SAmerica.nc
#G://ILAMB_runs_output/JULES_ORIG/JULES_modeldata_monthly.nc
#G://ILAMB_runs_output/INLAND_ORIG/INLAND_modeldata.nc
#G://CARDAMOM_Brazil/processed/20200924/Brazil_1deg_monthly_nopotAGB_2001_2017.nc
#G://ILAMBbeta/ILAMB_beta_tutorial/MODELS/INLAND/INLAND.nc
#G://CARDAMOM_ILAMB/CARDAMOM_Brazil_1x1_2001_2017_v1.0.nc

#extract NPP flux variable and convert to rasterer
inland_nppwood <- stack('G://ILAMB_runs_output/CSSP_stippling/INLAND/INLAND_monthly_1x1_SAmerica.nc',varname="alwood")
class(getZ(inland_nppwood))
plot(inland_nppwood$X2001.01.15)

inland_outputwood <- stack('G://ILAMB_runs_output/CSSP_stippling/INLAND/INLAND_monthly_1x1_SAmerica.nc',varname="fallw")
class(getZ(inland_outputwood))
plot(inland_outputwood$X2001.01.15)

#extract MTT wood variable and convert to rasterer
inland_cwood <- stack('G://ILAMB_runs_output/CSSP_stippling/INLAND/INLAND_monthly_1x1_SAmerica.nc',varname="woodbio")
class(getZ(inland_cwood))
plot(inland_cwood$X2001.01.15)

#extract raster images of two periods of interest
inland_nppwood_01_09 <- inland_nppwood[[which(getZ(inland_nppwood) >= as.Date("2001-01-15") & getZ(inland_nppwood) <= as.Date("2009-12-15"))]]
inland_outputwood_01_09 <- inland_outputwood[[which(getZ(inland_outputwood) >= as.Date("2001-01-15") & getZ(inland_outputwood) <= as.Date("2009-12-15"))]]
inland_cwood_01_09 <- inland_cwood[[which(getZ(inland_cwood) >= as.Date("2001-01-15") & getZ(inland_cwood) <= as.Date("2016-12-15"))]]

#calculate mean over time
inland_nppwood_01_09_mean <- stackApply(inland_nppwood_01_09, indices =  rep(1,nlayers(inland_nppwood_01_09)), fun = "mean")
inland_outputwood_01_09_mean <- stackApply(inland_outputwood_01_09, indices =  rep(1,nlayers(inland_outputwood_01_09)), fun = "mean")
inland_cwood_01_09_mean <- stackApply(inland_cwood_01_09, indices =  rep(1,nlayers(inland_cwood_01_09)), fun = "mean")
plot(inland_nppwood_01_09_mean,main='INLAND Wood Productivity 2001-2009 g/m2/day')

#crop and mask brzilian amazon
inland_nppwood_01_09_gm2d_amazon_brazil_crop <- crop(inland_nppwood_01_09_mean, extent(amazon_mask_pol))
inland_nppwood_01_09_gm2d_amazon_brazil <- mask(inland_nppwood_01_09_gm2d_amazon_brazil_crop, amazon_mask_pol)

inland_outputwood_01_09_gm2d_amazon_brazil_crop <- crop(inland_outputwood_01_09_mean, extent(amazon_mask_pol))
inland_outputwood_01_09_gm2d_amazon_brazil <- mask(inland_outputwood_01_09_gm2d_amazon_brazil_crop, amazon_mask_pol)

inland_cwood_01_09_gm2d_amazon_brazil_crop <- crop(inland_cwood_01_09_mean, extent(amazon_mask_pol))
inland_cwood_01_09_gm2d_amazon_brazil <- mask(inland_cwood_01_09_gm2d_amazon_brazil_crop, amazon_mask_pol)

par(mfrow=c(2,2))
plot(inland_nppwood_01_09_gm2d_amazon_brazil,main='INLAND Brazilian Amazon Wood Productivity 2001-2009 g/m2/day')
plot(inland_outputwood_01_09_gm2d_amazon_brazil,main='INLAND Brazilian Amazon Wood Mortality 2001-2009 g/m2/day')
plot(inland_cwood_01_09_gm2d_amazon_brazil,main='INLAND Brazilian Amazon Wood carbon 2001-2009 g/m2')

hist(inland_nppwood_01_09_gm2d_amazon_brazil)
hist(inland_outputwood_01_09_gm2d_amazon_brazil)
hist(inland_cwood_01_09_gm2d_amazon_brazil)

cellStats(inland_nppwood_01_09_gm2d_amazon_brazil, 'mean');cellStats(inland_nppwood_01_09_gm2d_amazon_brazil, 'sd')
cellStats(inland_outputwood_01_09_gm2d_amazon_brazil, 'mean');cellStats(inland_outputwood_01_09_gm2d_amazon_brazil, 'sd')
cellStats(inland_cwood_01_09_gm2d_amazon_brazil, 'mean');cellStats(inland_cwood_01_09_gm2d_amazon_brazil, 'sd')

#done
#JULES----
nc_jules <- nc_open('G://ILAMB_runs_output/CSSP2_proposal/MODELS/JULES/JULES_modeldata.nc')
print(nc_jules) #print properties
nc_close(nc_jules) #close nc file

#G://ILAMB_runs_output/CSSP_stippling/DATA/benchmark/CARDAMOM_monthly_1x1_SAmerica.nc
#G://ILAMB_runs_output/JULES_ORIG/JULES_modeldata_monthly.nc
#G://ILAMB_runs_output/INLAND_ORIG/INLAND_modeldata.nc
#G://CARDAMOM_Brazil/processed/20200924/Brazil_1deg_monthly_nopotAGB_2001_2017.nc
#G://ILAMBbeta/ILAMB_beta_tutorial/MODELS/INLAND/INLAND.nc
#G://CARDAMOM_ILAMB/CARDAMOM_Brazil_1x1_2001_2017_v1.0.nc

#extract NPP flux variable and convert to raster
jules_nppwood <- stack('G://ILAMB_runs_output/CSSP_stippling/JULES/JULES_monthly_1x1_SAmerica.nc',varname="alwood")
class(getZ(jules_nppwood))
plot(jules_nppwood$X2001.01.15)

jules_outputwood <- stack('G://ILAMB_runs_output/CSSP_stippling/JULES/JULES_monthly_1x1_SAmerica.nc',varname="fallw")
class(getZ(jules_outputwood))
plot(jules_outputwood$X2001.01.15)

#extract MTT wood variable and convert to rasterer
jules_cwood <- stack('G://ILAMB_runs_output/CSSP_stippling/JULES/JULES_monthly_1x1_SAmerica.nc',varname="woodbio")
class(getZ(jules_cwood))
plot(jules_cwood$X2001.01.15)

#extract raster images of two periods of interest
jules_nppwood_01_09 <- jules_nppwood[[which(getZ(jules_nppwood) >= as.Date("2001-01-15") & getZ(jules_nppwood) <= as.Date("2009-12-15"))]]
jules_outputwood_01_09 <- jules_outputwood[[which(getZ(jules_outputwood) >= as.Date("2001-01-15") & getZ(jules_outputwood) <= as.Date("2009-12-15"))]]
jules_cwood_01_09 <- jules_cwood[[which(getZ(jules_cwood) >= as.Date("2001-01-15") & getZ(jules_cwood) <= as.Date("2016-12-15"))]]

#calculate mean over time
jules_nppwood_01_09_mean <- stackApply(jules_nppwood_01_09, indices =  rep(1,nlayers(jules_nppwood_01_09)), fun = "mean")
jules_outputwood_01_09_mean <- stackApply(jules_outputwood_01_09, indices =  rep(1,nlayers(jules_outputwood_01_09)), fun = "mean")
jules_cwood_01_09_mean <- stackApply(jules_cwood_01_09, indices =  rep(1,nlayers(jules_cwood_01_09)), fun = "mean")

#crop and mask brzilian amazon
jules_nppwood_01_09_gm2d_amazon_brazil_crop <- crop(jules_nppwood_01_09_mean, extent(amazon_mask_pol))
jules_nppwood_01_09_gm2d_amazon_brazil <- mask(jules_nppwood_01_09_gm2d_amazon_brazil_crop, amazon_mask_pol)

jules_outputwood_01_09_gm2d_amazon_brazil_crop <- crop(jules_outputwood_01_09_mean, extent(amazon_mask_pol))
jules_outputwood_01_09_gm2d_amazon_brazil <- mask(jules_outputwood_01_09_gm2d_amazon_brazil_crop, amazon_mask_pol)

jules_cwood_01_09_gm2d_amazon_brazil_crop <- crop(jules_cwood_01_09_mean, extent(amazon_mask_pol))
jules_cwood_01_09_gm2d_amazon_brazil <- mask(jules_cwood_01_09_gm2d_amazon_brazil_crop, amazon_mask_pol)

par(mfrow=c(2,2))
plot(jules_nppwood_01_09_gm2d_amazon_brazil,main='JULES Brazilian Amazon Wood Productivity 2001-2009 g/m2/day')
plot(jules_outputwood_01_09_gm2d_amazon_brazil,main='JULES Brazilian Amazon Wood Mortality 2001-2009 g/m2/day')
plot(jules_cwood_01_09_gm2d_amazon_brazil,main='JULES Brazilian Amazon Wood carbon 2001-2009 g/m2')

hist(jules_nppwood_01_09_gm2d_amazon_brazil)
hist(jules_outputwood_01_09_gm2d_amazon_brazil)
hist(jules_cwood_01_09_gm2d_amazon_brazil)

cellStats(jules_nppwood_01_09_gm2d_amazon_brazil, 'mean');cellStats(jules_nppwood_01_09_gm2d_amazon_brazil, 'sd')
cellStats(jules_outputwood_01_09_gm2d_amazon_brazil, 'mean');cellStats(jules_outputwood_01_09_gm2d_amazon_brazil, 'sd')
cellStats(jules_cwood_01_09_gm2d_amazon_brazil, 'mean');cellStats(jules_cwood_01_09_gm2d_amazon_brazil, 'sd')

#done

#compare----
leeds_cardamom_npp <- overlay(leeds_woodprod_00_09_gm2d_amazon_brazil, cardamom_nppwood_01_09_gm2d_amazon_brazil, fun=function(x, y) { sqrt((y-x)/x) } )
inland_cwood_012001 <- inland_cwood$X2001.01.15
plot(inland_cwood_012001)
plot(biomass_amazon)

print(amazon_mask)
plot(biomass_amazon_mask_pol,add=T)

#create mask
e <- extent(inland_cwood_012001) #use south america extent
r <- raster(e)
res(r)<- 0.5 #use resolution to mirror other masks in ilamb/ could be 1x1
values(r) <- 1
crs(r)<-crs(amazon_mask)
r
plot(r)

r_amazon <- mask(r, biomass_amazon_mask_pol, updatevalue=0)
r_amazon
summary(r_amazon)
plot(r_amazon)

writeRaster(r_amazon, "ilamb-mask-AMAZONIA.tif", overwrite=TRUE, format="GTiff")
writeRaster(r_amazon, "mask-AMAZONIA-remap.nc", overwrite=TRUE, format="CDF",varname='ids',datatype='LOG1S', xname="lon", yname="lat")

r_southamerica <- mask(r, cardamom_agb_0101_mask_pol, updatevalue=0)
r_southamerica
summary(r_southamerica)
plot(r_southamerica)

#writeRaster(r_amazon, "ilamb-mask-AMAZONIA.tif", overwrite=TRUE, format="GTiff")
writeRaster(r_southamerica, "mask-SOUTHAMERICA-remap.nc", overwrite=TRUE, format="CDF",varname='ids',datatype='LOG1S', xname="lon", yname="lat")

amazonia_mask_ilamb <- mask(amazon_mask, biomass_amazon_mask_pol, updatevalue=0, inverse=T)
summary(amazonia_mask_ilamb)
plot(amazonia_mask_ilamb)

#done
#covert to nc file for ILAMB----
writeRaster(amazon_mask, "r_amazon.nc", overwrite=TRUE, format="CDF", xname="lon", yname="lat")

nc_r_amazon <- nc_open('./r_amazon.nc')
print(nc_r_amazon) #print properties
nc_close(nc_r_amazon) #close nc file

nc_amazonmask <- nc_open('R://ILAMB_tutorial/ilamb-mask-AMAZON.nc')
print(nc_amazonmask) #print properties
ncvar_get(nc_amazonmask)
nc_close(nc_amazonmask) #close nc file

summary(r_amazon)
plot(r_amazon)

ymax(r_amazon)-ymin(r_amazon)

lat_dim= ymax(r_amazon)-ymin(r_amazon)
long_dim= xmax(r_amazon)-xmin(r_amazon)

latitude = seq(ymin(r_amazon),ymax(r_amazon), length.out = lat_dim)
longitude = seq(xmin(r_amazon),xmax(r_amazon), length.out = long_dim)

amazon_mask_nc = array(NA, dim=c(long_dim,lat_dim,nos_quantile))
#check s america files----
nc_card <- nc_open('G://ILAMB_runs_output/CSSP_stippling/DATA/benchmark/CARDAMOM_monthly_1x1_SAmerica.nc')
print(nc_card) #print properties
nc_close(nc_card) #close nc file

cardamom_agb <- stack('G://ILAMB_runs_output/CSSP_stippling/DATA/benchmark/CARDAMOM_monthly_1x1_SAmerica.nc',varname="npptot")
class(getZ(cardamom_agb))
plot(cardamom_agb$X2001.01.15)

cardamom_agb_0101<-cardamom_agb$X2001.01.15
plot(cardamom_agb_0101)
cardamom_agb_0101_mask <- cardamom_agb_0101 > -Inf
cardamom_agb_0101_mask_pol <- rasterToPolygons(cardamom_agb_0101_mask, dissolve=TRUE)
plot(cardamom_agb_0101_mask_pol)
plot(biomass_amazon_mask_pol,add=T)


nc_inland <- nc_open('R://INLAND_monthly_1x1_SAmerica.nc', write=TRUE )
print(nc_inland) #print properties
attributes(nc_inland)$names
attributes(nc_inland$var)$names[17]
#attributes(nc_inland$var)$names[17]
ncatt_get(nc_inland, attributes(nc_inland$var)$names[17])

ncatt_put(nc_inland, "npptot", "units", "g.m-2.d-1")

old_varname <-'npptote'
new_varname <- 'npptot'

nc_inland <- ncvar_rename( nc_inland, old_varname, new_varname )
nc_close(nc_inland) #close nc file

inland_agb <- stack('G://ILAMB_runs_output/CSSP_stippling/INLAND/INLAND_monthly_1x1_SAmerica.nc',varname="woodbio")
class(getZ(inland_agb))
plot(inland_agb$X2001.01.15)

nc_jules <- nc_open('G://ILAMB_runs_output/CSSP_stippling/JULES/JULES_monthly_1x1_SAmerica.nc')
print(nc_jules) #print properties
nc_close(nc_jules) #close nc file

jules_agb <- stack('G://ILAMB_runs_output/CSSP_stippling/JULES/JULES_monthly_1x1_SAmerica.nc',varname="litC")
class(getZ(jules_agb))
plot(jules_agb$X2001.01.15)


#create new file----
xvals <- seq(-84,-32, length.out = 104)
yvals <- seq(-57,14, length.out = 142)
nx <- 104
ny <- 142

cnames   <-c("amazonia")
nstrings <- length(cnames)
## define dimension
xdim <- ncdim_def( "lon", units="degrees_east", xvals)
ydim <- ncdim_def( "lat", units="degrees_north", yvals)

dimnchar   <- ncdim_def("nb",   "", 1:2, create_dimvar=FALSE )
dimregion <- ncdim_def("n", "", 1:nstrings, create_dimvar=FALSE )

## make var
mv <--9.223372e+18 # missing value
var_mask <- ncvar_def( 'ids', "", list(xdim,ydim), mv,prec='integer')
varregion  <- ncvar_def("labels", "", list(dimregion) ,prec="char")

## make output file
output_fname <- 'new_amazon_mask.nc'
mask_new <- nc_create( output_fname, list(var_mask,varregion),force_v4=TRUE)
mask_new

## add data
ncvar_put( mask_new, "labels", cnames, verbose=TRUE )

data_mask <- array(0.,dim=c(72,62))

ncvar_put( mask_new, var_mask, data_mask, start=c(10,73), count=c(72,61))

ncatt_put( mask_new, "ids", "labels", "labels", prec="text" )

attributes(mask_new)$names
attributes(mask_new$var)$names
ncatt_get(mask_new, attributes(mask_new$var)$names[1])

nc_close(mask_new)

amazonia_mask_test <- raster('./new_amazon_mask.nc')
summary(amazonia_mask_test)
plot(amazonia_mask_test)
nc_amazonmask_test <- nc_open('./new_amazon_mask.nc') #gcel
nc_amazonmask_test
nc_close(nc_amazonmask_test)
#varcolors  <- ncvar_def("colours", "", list(dimnchar,dimregion) )


#done