
###
## Script to process default CARDAMOM output files (RData) into iLAMB digestible netcdf

# load needed libraries
library(ncdf4)

# set working directory
# setwd("M:/CARDAMOM/CARDAMOM_DEV/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/Rainfor_woody_biomass_productivity_and_mortality")
# Load needed functions 
library(zoo) ; library(compiler); library(raster)
source("M:/CARDAMOM/CARDAMOM_DEV/CARDAMOM/R_functions/read_binary_file_format.r")

# # load the CARDAMOM files
# load("./infofile.RData")
# load("./RESULTS_PROCESSED/Rainfor_woody_biomass_productivity_and_mortality_stock_flux.RData")
# load("./RESULTS_PROCESSED/Rainfor_woody_biomass_productivity_and_mortality_parameter_maps.RData")

#getwd()
drivers_270_eso = read_binary_file_format("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/esa_cci_agb_only_subset//DATA/esa_cci_agb_only_subset_00270.bin")
drivers_270_ban = read_binary_file_format("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/Rainfor_woody_biomass_annual_subset/DATA/Rainfor_woody_biomass_annual_subset_00270.bin")
drivers_270_bin = read_binary_file_format("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/Rainfor_woody_biomass_initial_subset/DATA/Rainfor_woody_biomass_initial_subset_00270.bin")
drivers_270_bn = read_binary_file_format("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/Rainfor_woody_biomass_2010_subset/DATA/Rainfor_woody_biomass_2010_subset_00270.bin")


drivers_270_bin$parpriors

drivers_270_eso$obs[,13]
drivers_270_ban$obs[,13]
drivers_270_bn$obs[,13]

esa_cci_agb_17 <- brick('G://AGB/ESA_CCI_BIOMASS/ESA_CCI_AGB_1deg/AGB_map_MgCha_2017.tif')
plot(esa_cci_agb_17)
summary(esa_cci_agb_17)[3]
rainfor_biomass_16 <- brick('G://cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/1deg_converted/wood_biomass/wood_biomass_gCm2_2016.tif')
plot(rainfor_biomass_16)
summary(rainfor_biomass_16)[3]
