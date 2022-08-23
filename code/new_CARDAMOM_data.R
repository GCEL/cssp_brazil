# load needed libraries
library(ncdf4)

# set working directory
setwd("/home/lsmallma/gcel/CARDAMOM_Brazil/inputs/20200604/")
# Load needed functions 
library(zoo) ; library(compiler)
source("~/WORK/GREENHOUSE/models/CARDAMOM/R_functions/read_binary_file_format.r")

# load the CARDAMOM files
load("R:/CARDAMOM_outputs/DALEC_CDEA_ACM2_BUCKET_MHMCMC/NoRainfor_woody_productivity_mortality/RESULTS_PROCESSED/NoRainfor_woody_productivity_mortality_parameter_maps.RData")

