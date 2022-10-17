#Leeds data----
#wood productivity t/ha/year
woodprod_00_09 <- brick('R://brazil_leeds_maps/WoodyProductivity20002009_Mg_perHa_perYear_111km.tif')
woodprod_10_16 <- brick('R://brazil_leeds_maps/WoodyProductivity20102016_Mg_perHa_perYear_111km.tif')
#biomass mortality t/ha/year
biommort_00_09 <- brick('R://brazil_leeds_maps/BiomassMortality_20002009_Mg_perHa_perYear_111km.tif')
biommort_10_16 <- brick('R://brazil_leeds_maps/BiomassMortality_20102016_Mg_perHa_perYear_111km.tif')
#biomass t/ha
biomass_amazon <- brick('R://brazil_leeds_maps/AbovegroundBiomass_Mg_perHa_111km.tif')

#functions for modification
thayr_to_gCm2day_fun <- function(x) {
  x *0.48 * (100/365.25) }
tha_to_gCm2_fun <- function(x) {
  x *0.48 * 100 }
bgb_inclusion <- function(x){
  x+(0.489 * x ** 0.89)
}

#include BGB to wood dynamics
woodprod_00_09 <- calc (woodprod_00_09,bgb_inclusion)
woodprod_10_16 <- calc (woodprod_10_16,bgb_inclusion)
biommort_00_09 <- calc (biommort_00_09,bgb_inclusion)
biommort_10_16 <- calc (biommort_10_16,bgb_inclusion)
biomass_amazon <- calc (biomass_amazon,bgb_inclusion)

#convert to gcm2
woodprod_00_09_gCm2d <- calc(woodprod_00_09, thayr_to_gCm2day_fun)
woodprod_10_16_gCm2d <- calc(woodprod_10_16, thayr_to_gCm2day_fun)
biommort_00_09_gCm2d <- calc(biommort_00_09, thayr_to_gCm2day_fun)
biommort_10_16_gCm2d <- calc(biommort_10_16, thayr_to_gCm2day_fun)
biomass_amazon_gCm2 <- calc(biomass_amazon, tha_to_gCm2_fun)


plot(woodprod_00_09_gCm2d,main='Woody Productivity 2000-2009 gC/m2/day')
plot(woodprod_10_16_gCm2d,main='Woody Productivity 2010-2016 gC/m2/day')
plot(biommort_00_09_gCm2d,main='Woody Mortality 2000-2009 gC/m2/day')
plot(biommort_10_16_gCm2d,main='Woody Mortality 2010-2016 gC/m2/day')
plot(biomass_amazon_gCm2,main='Woody Biomass gC/m2')


# amazon_woodprod_00_09_gm2d <- stack(replicate(10,woodprod_00_09_gm2d))
# amazon_woodprod_10_16_gm2d <- stack(replicate(7,woodprod_10_16_gm2d))
# 
# amazon_woodprod_01_09_monthly <- stack(replicate(108,woodprod_00_09_gm2d))
# amazon_woodprod_10_16_monthly <- stack(replicate(84,woodprod_10_16_gm2d))
# amazon_woodprod_10_monthly <- stack(replicate(12,woodprod_10_16_gm2d))

#Stack annual maps to 
amazon_woodprod_gCm2d_01_09_monthly <- stack(replicate(108,woodprod_00_09_gCm2d))
amazon_woodprod_gCm2d_10_16_monthly <- stack(replicate(84,woodprod_10_16_gCm2d))

amazon_biommort_gCm2d_01_09_monthly <- stack(replicate(108,biommort_00_09_gCm2d))
amazon_biommort_gCm2d_10_16_monthly <- stack(replicate(84,biommort_10_16_gCm2d))

amazon_woodprod_gCm2d_01_16_monthly <- stack(amazon_woodprod_gCm2d_01_09_monthly,amazon_woodprod_gCm2d_10_16_monthly)
amazon_biommort_gCm2d_01_16_monthly <- stack(amazon_biommort_gCm2d_01_09_monthly,amazon_biommort_gCm2d_10_16_monthly)
amazon_biomass_gCm2_01_16_monthly <- stack(replicate(192,biomass_amazon_gCm2))

# amazon_biommort_01_09_monthly <- stack(replicate(108,biommort_00_09_gm2d))
# amazon_biommort_10_16_monthly <- stack(replicate(84,biommort_10_16_gm2d))
# amazon_biommort_10_monthly <- stack(replicate(12,biommort_10_16_gm2d))


# amazon_woodprod_01_10_gm2d <- stack(woodprod_00_09_gm2d,woodprod_00_09_gm2d,woodprod_00_09_gm2d,woodprod_00_09_gm2d,woodprod_00_09_gm2d,
#                                     woodprod_00_09_gm2d,woodprod_00_09_gm2d,woodprod_00_09_gm2d,woodprod_00_09_gm2d,woodprod_10_16_gm2d)
# amazon_biommort_01_10_gm2d <- stack(biommort_00_09_gm2d,biommort_00_09_gm2d,biommort_00_09_gm2d,biommort_00_09_gm2d,biommort_00_09_gm2d,
#                                     biommort_00_09_gm2d,biommort_00_09_gm2d,biommort_00_09_gm2d,biommort_00_09_gm2d,biommort_10_16_gm2d)

# amazon_woodprod_01_16_gm2d <- stack(woodprod_00_09_gm2d,woodprod_00_09_gm2d,woodprod_00_09_gm2d,woodprod_00_09_gm2d,woodprod_00_09_gm2d,
#                                     woodprod_00_09_gm2d,woodprod_00_09_gm2d,woodprod_00_09_gm2d,woodprod_00_09_gm2d,woodprod_10_16_gm2d,
#                                     woodprod_10_16_gm2d,woodprod_10_16_gm2d,woodprod_10_16_gm2d,woodprod_10_16_gm2d,woodprod_10_16_gm2d,
#                                     woodprod_10_16_gm2d)
# amazon_biommort_01_16_gm2d <- stack(biommort_00_09_gm2d,biommort_00_09_gm2d,biommort_00_09_gm2d,biommort_00_09_gm2d,biommort_00_09_gm2d,
#                                     biommort_00_09_gm2d,biommort_00_09_gm2d,biommort_00_09_gm2d,biommort_00_09_gm2d,biommort_10_16_gm2d,
#                                     biommort_10_16_gm2d,biommort_10_16_gm2d,biommort_10_16_gm2d,biommort_10_16_gm2d,biommort_10_16_gm2d,
#                                     biommort_10_16_gm2d)
# 
# amazon_biomass_01_16_monthly <- stack(replicate(192,biomass_amazon_gm2))

# amazon_woodprod_01_16_monthly <- stack(amazon_woodprod_01_09_monthly,amazon_woodprod_10_16_monthly)
# amazon_biommort_01_16_monthly <- stack(amazon_biommort_01_09_monthly,amazon_biommort_10_16_monthly)
# 
# amazon_woodprod_01_10_monthly <- stack(amazon_woodprod_01_09_monthly,amazon_woodprod_10_monthly)
# amazon_biommort_01_10_monthly <- stack(amazon_biommort_01_09_monthly,amazon_biommort_10_monthly)


# rainfor_amazon_01_10 <- stack(amazon_woodprod_01_10_gm2d)


layer_dates <- names(cardamom_nppwood)
layer_dates <- layer_dates[1:192]

# layer_date_format <- getZ(cardamom_nppwood)
# layer_date_format <- layer_date_format[1:192]

#class(layer_date_format)
#layer_dates <- as.Date(layer_dates)

#layer_date_format_01_10 <- getZ(sa_data)

# names(amazon_woodprod_01_16_monthly) <- layer_dates
# names(amazon_biommort_01_16_monthly) <- layer_dates
# names(amazon_biomass_01_16_monthly) <- layer_dates
# 
# names(amazon_woodprod_01_10_monthly) <- layer_dates_01_10
# names(amazon_biommort_01_10_monthly) <- layer_dates_01_10
# 
# amazon_woodprod_01_16_monthly <- setZ(amazon_woodprod_01_16_monthly,layer_date_format,"Date")
# amazon_biommort_01_16_monthly <- setZ(amazon_biommort_01_16_monthly,layer_date_format,"Date")
# amazon_biomass_01_16_monthly <- setZ(amazon_biomass_01_16_monthly,layer_date_format,"Date")
# 
# amazon_woodprod_01_10_monthly <- setZ(amazon_woodprod_01_10_monthly,layer_date_format_01_10,"Date")
# amazon_biommort_01_10_monthly <- setZ(amazon_biommort_01_10_monthly,layer_date_format_01_10,"Date")

names(amazon_woodprod_gCm2d_01_16_monthly) <- layer_dates
names(amazon_biommort_gCm2d_01_16_monthly) <- layer_dates
names(amazon_biomass_gCm2_01_16_monthly) <- layer_dates

amazon_woodprod_gCm2d_01_16_monthly <- setZ(amazon_woodprod_gCm2d_01_16_monthly,layer_date_format,"Date")
amazon_biommort_gCm2d_01_16_monthly <- setZ(amazon_biommort_gCm2d_01_16_monthly,layer_date_format,"Date")
amazon_biomass_gCm2_01_16_monthly <- setZ(amazon_biomass_gCm2_01_16_monthly,layer_date_format,"Date")

# writeRaster(amazon_woodprod_01_10_gm2d, "rainfor_nppwood_01_10.nc", overwrite=TRUE,
#             format="CDF",     varname="nppwood", varunit="g.m-2.d-1",force_v4=TRUE,
#             longname="NPP for wood", xname="lon",   yname="lat", zname='time', zunit='years since 2000-01-01')
# writeRaster(amazon_biommort_01_10_gm2d, "rainfor_outputwood_01_10.nc", overwrite=TRUE,
#             format="CDF",     varname="outputwood", varunit="g.m-2.d-1",force_v4=TRUE,
#             longname="Wood Mortality", xname="lon",   yname="lat", zname='time', zunit='years since 2000-01-01')

# writeRaster(amazon_woodprod_01_16_gm2d, "rainfor_nppwood_01_16.nc", overwrite=TRUE,
#             format="CDF",     varname="nppwood", varunit="g.m-2.d-1",force_v4=TRUE,
#             longname="NPP for wood", xname="lon",   yname="lat", zname='time', zunit='years since 2000-01-01')
# writeRaster(amazon_biommort_01_16_gm2d, "rainfor_outputwood_01_16.nc", overwrite=TRUE,
#             format="CDF",     varname="outputwood", varunit="g.m-2.d-1",force_v4=TRUE,
#             longname="Wood Mortality", xname="lon",   yname="lat", zname='time', zunit='years since 2000-01-01')

# writeRaster(amazon_woodprod_01_16_monthly, "rainfor_nppwood_01_16_monthly.nc", overwrite=TRUE,
#             format="CDF",     varname="nppwood", varunit="g.m-2.d-1",force_v4=TRUE,
#             longname="NPP for wood", xname="lon",   yname="lat", zname='time', zunit='days since 2001-1-15 00:00:00')
# writeRaster(amazon_biommort_01_16_monthly, "rainfor_outputwood_01_16_monthly.nc", overwrite=TRUE,
#             format="CDF",     varname="outputwood", varunit="g.m-2.d-1",force_v4=TRUE,
#             longname="Wood Mortality", xname="lon",   yname="lat", zname='time', zunit='days since 2001-1-15 00:00:00')
# writeRaster(amazon_biomass_01_16_monthly, "rainfor_biomasswood_01_16_monthly.nc", overwrite=TRUE,
#             format="CDF",     varname="WOOD", varunit="g.m-2",force_v4=TRUE,
#             longname="Wood Biomass", xname="lon",   yname="lat", zname='time', zunit='days since 2001-1-15 00:00:00')

# writeRaster(amazon_woodprod_01_10_monthly, "rainfor_nppwood_01_10_monthly.nc", overwrite=TRUE,
#             format="CDF",     varname="NPP_wood_flx", varunit="g.m-2.d-1",force_v4=TRUE,
#             longname="NPP for wood", xname="lon",   yname="lat", zname='time', zunit='days since 2001-1-15 00:00:00')
# writeRaster(amazon_biommort_01_10_monthly, "rainfor_outputwood_01_10_monthly.nc", overwrite=TRUE,
#             format="CDF",     varname="outputwood", varunit="g.m-2.d-1",force_v4=TRUE,
#             longname="Wood Mortality", xname="lon",   yname="lat", zname='time', zunit='days since 2001-1-15 00:00:00')


writeRaster(amazon_woodprod_gCm2d_01_16_monthly, "rainfor_nppwood_gCm2d_01_16_monthly.nc", overwrite=TRUE, 
            format="CDF",     varname="nppwood", varunit="g.m-2.d-1",force_v4=TRUE,
            longname="Wood Productivity", xname="lon",   yname="lat", zname='time', zunit='days since 2001-1-15 00:00:00')
writeRaster(amazon_biommort_gCm2d_01_16_monthly, "rainfor_outputwood_gCm2d_01_16_monthly.nc", overwrite=TRUE, 
            format="CDF",     varname="outputwood", varunit="g.m-2.d-1",force_v4=TRUE,
            longname="Wood Mortality", xname="lon",   yname="lat", zname='time', zunit='days since 2001-1-15 00:00:00')
writeRaster(amazon_biomass_gCm2_01_16_monthly, "rainfor_biomasswood_gCm2_01_16_monthly.nc", overwrite=TRUE, 
            format="CDF",     varname="WOOD", varunit="g.m-2",force_v4=TRUE,
            longname="Wood Biomass", xname="lon",   yname="lat", zname='time', zunit='days since 2001-1-15 00:00:00')

nc_amazon_rainfor_data <- nc_open('./rainfor_outputwood_01_10_monthly.nc') 
nc_amazon_rainfor_data
nc_close(nc_amazon_rainfor_data)
rainfor_outputwood <- stack('./rainfor_outputwood_01_10_monthly.nc',varname="outputwood")
rainfor_outputwood
#summary(rainfor_nppwood$X1)
class(getZ(rainfor_outputwood))

nc_time<-nc_open('./time.nc') 
nc_time
nc_close(nc_time)

nc_wrong_time<-nc_open('./wrongFormat.nc') 
nc_wrong_time
nc_close(nc_wrong_time)

nc_correct_time<-nc_open('./correctFormat.nc') 
nc_correct_time
nc_close(nc_correct_time)
