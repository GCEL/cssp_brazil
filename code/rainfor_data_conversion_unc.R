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
