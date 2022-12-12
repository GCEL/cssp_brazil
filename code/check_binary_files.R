
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
drivers_270_nwd = read_binary_file_format("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/no_woody_data_subset/DATA/no_woody_data_subset_00270.bin")
drivers_270_eso = read_binary_file_format("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/esa_cci_agb_only_subset//DATA/esa_cci_agb_only_subset_00270.bin")
drivers_270_ban = read_binary_file_format("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/Rainfor_woody_biomass_annual_subset/DATA/Rainfor_woody_biomass_annual_subset_00270.bin")
drivers_270_bin = read_binary_file_format("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/Rainfor_woody_biomass_initial_subset/DATA/Rainfor_woody_biomass_initial_subset_00270.bin")
drivers_270_bn = read_binary_file_format("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/Rainfor_woody_biomass_2010_subset/DATA/Rainfor_woody_biomass_2010_subset_00270.bin")

names(drivers_270_nwd)

drivers_270_nwd$obs

drivers_270_nwd$obs[,13]
drivers_270_eso$obs[,13]
drivers_270_ban$obs[,13]
drivers_270_bn$obs[,13]


drivers_696_nbe = read_binary_file_format("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/nbe_data_alone/DATA/nbe_data_alone_00696.bin")
drivers_696_nbe$obs[,35]

drivers_270_eso$obs[,3]

drivers_696_laic = read_binary_file_format("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/no_woody_data_copernicus/DATA/no_woody_data_copernicus_00696.bin")
drivers_696_laic$obs[,3]
drivers_696_laim = read_binary_file_format("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/no_woody_data_modis/DATA/no_woody_data_modis_00696.bin")
drivers_696_laim$obs[,3]
drivers_270_laim = read_binary_file_format("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/no_woody_data_modis_subset/DATA/no_woody_data_modis_subset_00270.bin")
drivers_270_laim$obs[,3]
