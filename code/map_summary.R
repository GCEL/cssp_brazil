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

woodprod_00_09_gm2d <- calc(woodprod_00_09, thayr_to_gm2day_fun)
woodprod_10_16_gm2d <- calc(woodprod_10_16, thayr_to_gm2day_fun)
biommort_00_09_gm2d <- calc(biommort_00_09, thayr_to_gm2day_fun)
biommort_10_16_gm2d <- calc(biommort_10_16, thayr_to_gm2day_fun)

#summary
summary(woodprod_00_09)
str(woodprod_00_09)

#plots
par(mfrow=c(2,2))
plot(woodprod_00_09,main='Woody Productivity 2000-2009 Mg/Ha/Year')
plot(woodprod_10_16,main='Woody Productivity 2010-2016')
plot(biommort_00_09,main='Biomass Mortality 2000-2009')
plot(biommort_10_16,main='Biomass Mortality 2010-2016')

par(mfrow=c(1,1))
plot(biomass_amazon,main='Leeds Aboveground Biomass')

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


#models----
nc_card <- nc_open('G://CARDAMOM_ILAMB/CARDAMOM_Brazil_1x1_2001_2017_v1.0.nc')
print(nc_card) #print properties
nc_close(nc_card) #close nc file

#G://ILAMB_runs_output/CSSP_stippling/DATA/benchmark/CARDAMOM_monthly_1x1_SAmerica.nc
#G://ILAMB_runs_output/JULES_ORIG/JULES_modeldata_monthly.nc
#G://ILAMB_runs_output/INLAND_ORIG/INLAND_modeldata.nc
#G://CARDAMOM_Brazil/processed/20200924/Brazil_1deg_monthly_nopotAGB_2001_2017.nc
#G://ILAMBbeta/ILAMB_beta_tutorial/MODELS/INLAND/INLAND.nc
#G://CARDAMOM_ILAMB/CARDAMOM_Brazil_1x1_2001_2017_v1.0.nc

#extract NPP flux variable and convert to rasterer
cardamom_nppwood <- stack('G://CARDAMOM_ILAMB/CARDAMOM_Brazil_1x1_2001_2017_v1.0.nc',varname="NPP_wood_flx")
class(getZ(cardamom_nppwood))
plot(cardamom_nppwood$X2001.01.01)

cardamom_nppwood_1101<-cardamom_nppwood$X2001.01.01 #template to get data ignore values (green in plot of previous step)
cardamom_nppwood_1101[cardamom_nppwood_1101 > 500000] <- NA #ignore values above 500000
plot(cardamom_nppwood_1101) #view masked out values

cardamom_nppwood_masked <- mask(cardamom_nppwood, cardamom_nppwood_1101) #mask out data ignore values
plot(cardamom_nppwood_masked$X2001.01.01) #check successful masking

#extract raster images of two periods of interest
cardamom_nppwood_masked_01_10 <- cardamom_nppwood_masked[[which(getZ(cardamom_nppwood_masked) >= as.Date("2001-01-01") & getZ(cardamom_nppwood_masked) <= as.Date("2009-12-01"))]]
cardamom_nppwood_masked_10_16 <- cardamom_nppwood_masked[[which(getZ(cardamom_nppwood_masked) >= as.Date("2010-01-01") & getZ(cardamom_nppwood_masked) <= as.Date("2016-12-01"))]]

#calculate mean over time
cardamom_nppwood_masked_01_10_mean <- stackApply(cardamom_nppwood_masked_01_10, indices =  rep(1,nlayers(cardamom_nppwood_masked_01_10)), fun = "mean")
cardamom_nppwood_masked_10_16_mean <- stackApply(cardamom_nppwood_masked_10_16, indices =  rep(1,nlayers(cardamom_nppwood_masked_10_16)), fun = "mean")

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
amazon_mask <- brick('./data/ilamb-mask-AMAZON.nc')
amazon_mask[amazon_mask < 0] <- NA
amazon_mask <- amazon_mask> -Inf
plot(amazon_mask)
amazon_mask_pol <- rasterToPolygons(amazon_mask, dissolve = TRUE)
plot(amazon_mask_pol,add=T)

# biomass_amazon_mask <- biomass_amazon > -Inf
# biomass_amazon_mask_pol <- rasterToPolygons(biomass_amazon_mask, dissolve=TRUE)
# plot(biomass_amazon_mask_pol)
# 
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