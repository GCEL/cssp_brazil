
###
## Script to process default CARDAMOM output files (RData) into iLAMB digestible netcdf
### 

# load needed libraries
library(ncdf4)
library(raster)
library(compiler)

# set working directory
#setwd("/home/cnwobi/CARDAMOM/CARDAMOM_DEV/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/NoRainfor_woody_productivity_mortality")
#setwd("/home/cnwobi/CARDAMOM/CARDAMOM_DEV/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/Rainfor_woody_productivity_and_mortality")
# setwd("/home/cnwobi/CARDAMOM/CARDAMOM_DEV/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/Rainfor_woody_biomass_productivity_and_mortality")
#setwd("/home/cnwobi/CARDAMOM/CARDAMOM_DEV/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/Rainfor_woody_productivity")
#setwd("/home/cnwobi/CARDAMOM/CARDAMOM_DEV/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/Rainfor_woody_mortality")
#setwd("/home/cnwobi/CARDAMOM/CARDAMOM_DEV/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/Rainfor_woody_biomass")
#setwd("/home/cnwobi/CARDAMOM/CARDAMOM_DEV/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/Rainfor_woody_biomass_productivity_and_mortality")
#setwd("/home/cnwobi/CARDAMOM/CARDAMOM_DEV/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/NoRainfor_woody_productivity_mortality")
# setwd("/home/cnwobi/CARDAMOM/CARDAMOM_DEV/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/Rainfor_woody_biomass_subset")
# setwd("/home/cnwobi/CARDAMOM/CARDAMOM_DEV/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/Rainfor_woody_productivity_subset")
#setwd("/home/cnwobi/CARDAMOM/CARDAMOM_DEV/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/Rainfor_woody_mortality_subset")
#setwd("/home/cnwobi/CARDAMOM/CARDAMOM_DEV/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/Rainfor_woody_biomass_productivity_and_mortality_subset")
#setwd("/home/cnwobi/CARDAMOM/CARDAMOM_DEV/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/NoRainfor_woody_biomass_productivity_mortality_subset")
#setwd("/home/cnwobi/CARDAMOM/CARDAMOM_DEV/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/Rainfor_woody_productivity_and_mortality_subset")
setwd("M:/CARDAMOM/CARDAMOM_DEV/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/Rainfor_woody_biomass_and_mortality_subset")
#setwd("M:/CARDAMOM/CARDAMOM_DEV/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/Rainfor_woody_biomass_and_productivity_subset")

# Load needed functions 
library(zoo) ; library(compiler) ; library(raster)
source("M:/CARDAMOM/CARDAMOM_DEV/CARDAMOM/R_functions/read_binary_file_format.r")

# load the CARDAMOM files
load("./infofile.RData")
load(paste(PROJECT$results_processedpath,PROJECT$name,"_stock_flux.RData",sep=""))
#load(paste(PROJECT$results_processedpath,PROJECT$name,"_parameter_maps.RData",sep=""))

# set output name
#output_name = paste("/exports/csce/datastore/geos/users/cnwobi/ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/MODELS/",PROJECT$name,"_",PROJECT$start_year,"_",PROJECT$end_year,".nc",sep="")
#output_name = paste("/exports/csce/datastore/geos/users/cnwobi/ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/DATA/",PROJECT$name,"_",PROJECT$start_year,"_",PROJECT$end_year,".nc",sep="")
output_name = paste("/exports/csce/datastore/geos/users/cnwobi/ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/",PROJECT$name,"_",PROJECT$start_year,"_updated_",PROJECT$end_year,".nc",sep="")

# Time information
nos_years = length(c(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)))
steps_per_year = dim(grid_output$lai_m2m2)[3] / nos_years

# create lat / long axes, assumes WGS-84 grid
latitude = seq(PROJECT$latitude[1]+(PROJECT$resolution),PROJECT$latitude[2]-(PROJECT$resolution), length.out = PROJECT$lat_dim)
longitude = seq(PROJECT$longitude[1]+(PROJECT$resolution),PROJECT$longitude[2]-(PROJECT$resolution), length.out = PROJECT$long_dim)

# restricture each of the variable into a complete grid
nos_quantile = length(grid_output$num_quantiles) #  we will output the 2.5 %, 25%, 50 %, 75 %, 97.5 % quantiles

# STATES
LAI = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))
TOT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))
LAB = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))
FOL = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))
ROOT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))
WOOD = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))
LIT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))
SOIL = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))
SOILC = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))

# FLUXES
GPP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))
Ra = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))
Rh = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))
NPP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))
FIR = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))
HAR = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))
RECO = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))
NEE = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))
NBE = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))

# C emergent properties
#CUE = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))

# BIOPHYSICAL
CiCa = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))

# NPP allocation (labile, foliar, fine root, wood, gC/m2/day)
fNPP_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))
laNPP_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))
wNPP_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))
rNPP_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))

# # MRT (foliar, wood, fine root, litter, soil; years)
fMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))
wMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))
rMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))
laMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))
liMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))
sMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))
scMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))

#NAT_OUTPUT (foliar, root, wood; gC/m2/day)
fOUTPUT_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))
laOUTPUT_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))
liOUTPUT_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))
wOUTPUT_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))
rOUTPUT_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))
sOUTPUT_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))
scOUTPUT_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))
OUTPUT_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantile,length(PROJECT$model$timestep_days)))

# Fill the arrays
for (n in seq(1, length(PROJECT$sites))) {

     # Ensure the site has been processed
     if (is.na(grid_output$i_location[n]) == FALSE) {

         # Extract grid position
         i = grid_output$i_location[n]
         j = grid_output$j_location[n]

         # Read in site specific drivers
         drivers = read_binary_file_format(paste(PROJECT$datapath,PROJECT$name,"_",PROJECT$sites[n],".bin",sep=""))
         
         # STATES
         LAI[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$lai_m2m2[n,1:7,] # 2 in file is 25 %, 3 = 50 %, 4 = 75 %
         TOT[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$biomass_gCm2[n,1:7,]
         LAB[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$labile_gCm2[n,1:7,]
         FOL[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$foliage_gCm2[n,1:7,] 
         ROOT[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$roots_gCm2[n,1:7,]
         WOOD[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$wood_gCm2[n,1:7,]
         LIT[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$litter_gCm2[n,1:7,]
         SOIL[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$som_gCm2[n,1:7,]
         SOILC[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$dom_gCm2[n,1:7,]
         
         # FLUXES
         GPP[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$gpp_gCm2day[n,1:7,]
         Ra[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$rauto_gCm2day[n,1:7,]
         Rh[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$rhet_gCm2day[n,1:7,]
         NPP[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$npp_gCm2day[n,1:7,]
         FIR[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$fire_gCm2day[n,1:7,]
         HAR[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$harvest_gCm2day[n,1:7,]
         RECO[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$reco_gCm2day[n,1:7,]
         NEE[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$nee_gCm2day[n,1:7,]
         NBE[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$nbe_gCm2day[n,1:7,]
         
         # BIOPHYSICAL
         CiCa[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$CiCa[n,1:7,]
         
         # NPP (foliar, root, wood; gC/m2/day)
         fNPP_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$combined_alloc_foliage_gCm2day[n,1:7,]
         laNPP_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$alloc_labile_gCm2day[n,1:7,]
         rNPP_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$alloc_roots_gCm2day[n,1:7,]
         wNPP_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$alloc_wood_gCm2day[n,1:7,]
         
         # OUTPUT (foliar, root, wood; gC/m2/day)
         fOUTPUT_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$outflux_foliage_gCm2day[n,1:7,] 
         laOUTPUT_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$outflux_labile_gCm2day[n,1:7,] 
         liOUTPUT_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$outflux_litter_gCm2day[n,1:7,] 
         rOUTPUT_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$outflux_roots_gCm2day[n,1:7,]
         wOUTPUT_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$outflux_wood_gCm2day[n,1:7,]
         sOUTPUT_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$outflux_som_gCm2day[n,1:7,]
         scOUTPUT_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$outflux_dom_gCm2day[n,1:7,]
         OUTPUT_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$outflux_biomass_gCm2day[n,1:7,]
         #,c(1,3,4,5,7),
         # NPP (fraction) and MRT years are requested to have same number of time steps as stocks and fluxes
         # This is awkward as no easy way to repeat specific elements without loop for variables which have no meaningful value at sub-annual timescales
         # (and are therefore calculated as annuals)
         for (q in seq(1, nos_quantile)) {
              # MRT
              fMRT[grid_output$i_location[n],grid_output$j_location[n],q,] = rep(grid_output$MTT_annual_foliage_years[n,q,], each = steps_per_year)
              rMRT[grid_output$i_location[n],grid_output$j_location[n],q,] = rep(grid_output$MTT_annual_roots_years[n,q,], each = steps_per_year)
              wMRT[grid_output$i_location[n],grid_output$j_location[n],q,] = rep(grid_output$MTT_annual_wood_years[n,q,], each = steps_per_year)
              laMRT[grid_output$i_location[n],grid_output$j_location[n],q,] = rep(grid_output$MTT_annual_labile_years[n,q,], each = steps_per_year)
              liMRT[grid_output$i_location[n],grid_output$j_location[n],q,] = rep(grid_output$MTT_annual_litter_years[n,q,], each = steps_per_year)
              sMRT[grid_output$i_location[n],grid_output$j_location[n],q,] = rep(grid_output$MTT_annual_som_years[n,q,], each = steps_per_year)
              scMRT[grid_output$i_location[n],grid_output$j_location[n],q,] = rep(grid_output$MTT_annual_dom_years[n,q,], each = steps_per_year)

              # # For Auguste
              #CUE[grid_output$i_location[n],grid_output$j_location[n],q,] = rep(grid_parameters$parameters[grid_output$i_location[n],grid_output$j_location[n],2,q], each = length(PROJECT$model$timestep_days))
         } # loop quantiles
     } # Does the file exist / has it been processed

} # site loop

#CUE = 1 - CUE 

## define dimension
lat_dimen <- ncdim_def( "lat", units="degree north (-90->90)", latitude )
long_dimen <- ncdim_def( "lon", units="degree east (-180->180)", longitude )
time_dimen <- ncdim_def( "time", units="", 1:length(PROJECT$model$timestep_days))
year_dimen <- ncdim_def( "year", units="", 1:nos_years)

## define output variable
var0 = ncvar_def("Time", units = "d", longname = paste("Monthly time step given in days since 01/01/",PROJECT$start_year,sep=""), dim=list(time_dimen), missval = NA, prec="double", compression = 9)

## STATES
# LAI
var1 = ncvar_def("LAI_2pt5pc", unit="m2.m-2", longname = "Leaf Area Index - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var2 = ncvar_def("LAI_25pc", unit="m2.m-2", longname = "Leaf Area Index - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var3 = ncvar_def("LAI", unit="m2.m-2", longname = "Leaf Area Index - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var4 = ncvar_def("LAI_75pc", unit="m2.m-2", longname = "Leaf Area Index - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var5 = ncvar_def("LAI_97pt5pc", unit="m2.m-2", longname = "Leaf Area Index - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# Labile
var6 = ncvar_def("LAB_2pt5pc", unit="g.m-2", longname = "Labile C - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var7 = ncvar_def("LAB_25pc", unit="g.m-2", longname = "Labile C - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var8 = ncvar_def("LAB", unit="g.m-2", longname = "Labile C - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var9 = ncvar_def("LAB_75pc", unit="g.m-2", longname = "Labile C - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var10 = ncvar_def("LAB_97pt5pc", unit="g.m-2", longname = "Labile C - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# Foliar
var11 = ncvar_def("FOL_2pt5pc", unit="g.m-2", longname = "Foliage C - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var12 = ncvar_def("FOL_25pc", unit="g.m-2", longname = "Foliage C - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var13 = ncvar_def("FOL", unit="g.m-2", longname = "Foliage C - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var14 = ncvar_def("FOL_75pc", unit="g.m-2", longname = "Foliage C - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var15 = ncvar_def("FOL_97pt5pc", unit="g.m-2", longname = "Foliage C - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# Fine root
var16 = ncvar_def("ROOT_2pt5pc", unit="g.m-2", longname = "Fine roots C - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var17 = ncvar_def("ROOT_25pc", unit="g.m-2", longname = "Fine roots C - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var18 = ncvar_def("ROOT", unit="g.m-2", longname = "Fine roots C - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var19 = ncvar_def("ROOT_75pc", unit="g.m-2", longname = "Fine roots C - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var20 = ncvar_def("ROOT_97pt5pc", unit="g.m-2", longname = "Fine roots C - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# Wood
var21 = ncvar_def("WOOD_2pt5pc", unit="g.m-2", longname = "Wood C - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var22 = ncvar_def("WOOD_25pc", unit="g.m-2", longname = "Wood C - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var23 = ncvar_def("WOOD", unit="g.m-2", longname = "Wood C - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var24 = ncvar_def("WOOD_75pc", unit="g.m-2", longname = "Wood C - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var25 = ncvar_def("WOOD_97pt5pc", unit="g.m-2", longname = "Wood C - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# Foliar + fine root litter
var26 = ncvar_def("LIT_2pt5pc", unit="g.m-2", longname = "Fine litter C - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var27 = ncvar_def("LIT_25pc", unit="g.m-2", longname = "Fine litter C - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var28 = ncvar_def("LIT", unit="g.m-2", longname = "Fine litter C - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var29 = ncvar_def("LIT_75pc", unit="g.m-2", longname = "Fine litter C - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var30 = ncvar_def("LIT_97pt5pc", unit="g.m-2", longname = "Fine litter C - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# Soil organic matter
var31 = ncvar_def("SOIL_2pt5pc", unit="g.m-2", longname = "Soil C - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var32 = ncvar_def("SOIL_25pc", unit="g.m-2", longname = "Soil C - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var33 = ncvar_def("SOIL", unit="g.m-2", longname = "Soil C - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var34 = ncvar_def("SOIL_75pc", unit="g.m-2", longname = "Soil C - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var35 = ncvar_def("SOIL_97pt5pc", unit="g.m-2", longname = "Soil C - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# TotalC
var136 = ncvar_def("TOT_2pt5pc", unit="g.m-2", longname = "Total C - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var137 = ncvar_def("TOT_25pc", unit="g.m-2", longname = "Total C - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var138 = ncvar_def("TOT", unit="g.m-2", longname = "Total C - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var139 = ncvar_def("TOT_75pc", unit="g.m-2", longname = "Total C - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var140 = ncvar_def("TOT_97pt5pc", unit="g.m-2", longname = "Total C - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# Dead carbon mass
var141 = ncvar_def("SOIL_C_2pt5pc", unit="g.m-2", longname = "Dead carbon mass - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var142 = ncvar_def("SOIL_C_25pc", unit="g.m-2", longname = "Dead carbon mass - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var143 = ncvar_def("SOIL_C", unit="g.m-2", longname = "Dead carbon mass - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var144 = ncvar_def("SOIL_C_75pc", unit="g.m-2", longname = "Dead carbon mass - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var145 = ncvar_def("SOIL_C_97pt5pc", unit="g.m-2", longname = "Dead carbon mass - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)

## FLUXES
# GPP
var36 = ncvar_def("GPP_2pt5pc", unit="g.m-2.d-1", longname = "Gross Primary Productivity - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var37 = ncvar_def("GPP_25pc", unit="g.m-2.d-1", longname = "Gross Primary Productivity - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var38 = ncvar_def("GPP", unit="g.m-2.d-1", longname = "Gross Primary Productivity - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var39 = ncvar_def("GPP_75pc", unit="g.m-2.d-1", longname = "Gross Primary Productivity - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var40 = ncvar_def("GPP_97pt5pc", unit="g.m-2.d-1", longname = "Gross Primary Productivity - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# Autotrophic respiration
var41 = ncvar_def("Ra_2pt5pc", unit="g.m-2.d-1", longname = "Autotrophic Respiration - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var42 = ncvar_def("Ra_25pc", unit="g.m-2.d-1", longname = "Autotrophic Respiration - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var43 = ncvar_def("Ra", unit="g.m-2.d-1", longname = "Autotrophic Respiration - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var44 = ncvar_def("Ra_75pc", unit="g.m-2.d-1", longname = "Autotrophic Respiration - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var45 = ncvar_def("Ra_97pt5pc", unit="g.m-2.d-1", longname = "Autotrophic Respiration - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# Heterotrophic respiration
var46 = ncvar_def("Rh_2pt5pc", unit="g.m-2.d-1", longname = "Heterotrophic Respiration - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var47 = ncvar_def("Rh_25pc", unit="g.m-2.d-1", longname = "Heterotrophic Respiration - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var48 = ncvar_def("Rh", unit="g.m-2.d-1", longname = "Heterotrophic Respiration - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var49 = ncvar_def("Rh_75pc", unit="g.m-2.d-1", longname = "Heterotrophic Respiration - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var50 = ncvar_def("Rh_97pt5pc", unit="g.m-2.d-1", longname = "Heterotrophic Respiration - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# Net Primary Productivity
var51 = ncvar_def("NPP_2pt5pc", unit="g.m-2.d-1", longname = "Net Primary Productivity - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var52 = ncvar_def("NPP_25pc", unit="g.m-2.d-1", longname = "Net Primary Productivity - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var53 = ncvar_def("NPP", unit="g.m-2.d-1", longname = "Net Primary Productivity - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var54 = ncvar_def("NPP_75pc", unit="g.m-2.d-1", longname = "Net Primary Productivity - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var55 = ncvar_def("NPP_97pt5pc", unit="g.m-2.d-1", longname = "Net Primary Productivity - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# Fire emissions
var56 = ncvar_def("FIR_2pt5pc", unit="g.m-2.d-1", longname = "Fire - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var57 = ncvar_def("FIR_25pc", unit="g.m-2.d-1", longname = "Fire - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var58 = ncvar_def("FIR", unit="g.m-2.d-1", longname = "Fire - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var59 = ncvar_def("FIR_75pc", unit="g.m-2.d-1", longname = "Fire - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var60 = ncvar_def("FIR_97pt5pc", unit="g.m-2.d-1", longname = "Fire - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# Ecosystem respiration
var121 = ncvar_def("RECO_2pt5pc", unit="g.m-2.d-1", longname = "Ecosystem Respiration - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var122 = ncvar_def("RECO_25pc", unit="g.m-2.d-1", longname = "Ecosystem Respiration - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var123 = ncvar_def("RECO", unit="g.m-2.d-1", longname = "Ecosystem Respiration - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var124 = ncvar_def("RECO_75pc", unit="g.m-2.d-1", longname = "Ecosystem Respiration - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var125 = ncvar_def("RECO_97pt5pc", unit="g.m-2.d-1", longname = "Ecosystem Respiration - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# Net Ecosystem Exchange
var126 = ncvar_def("NEE_2pt5pc", unit="g.m-2.d-1", longname = "Net Ecosystem Exchange - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var127 = ncvar_def("NEE_25pc", unit="g.m-2.d-1", longname = "Net Ecosystem Exchange - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var128 = ncvar_def("NEE", unit="g.m-2.d-1", longname = "Net Ecosystem Exchange - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var129 = ncvar_def("NEE_75pc", unit="g.m-2.d-1", longname = "Net Ecosystem Exchange - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var130 = ncvar_def("NEE_97pt5pc", unit="g.m-2.d-1", longname = "Net Ecosystem Exchange - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# Net Biome Exchange
var131 = ncvar_def("NBE_2pt5pc", unit="g.m-2.d-1", longname = "Net Biome Exchange - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var132 = ncvar_def("NBE_25pc", unit="g.m-2.d-1", longname = "Net Biome Exchange - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var133 = ncvar_def("NBE", unit="g.m-2.d-1", longname = "Net Biome Exchange - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var134 = ncvar_def("NBE_75pc", unit="g.m-2.d-1", longname = "Net Biome Exchange - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var135 = ncvar_def("NBE_97pt5pc", unit="g.m-2.d-1", longname = "Net Biome Exchange - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# Harvest emissions
var146 = ncvar_def("HAR_2pt5pc", unit="g.m-2.d-1", longname = "Harvest - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var147 = ncvar_def("HAR_25pc", unit="g.m-2.d-1", longname = "Harvest - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var148 = ncvar_def("HAR", unit="g.m-2.d-1", longname = "Harvest - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var149 = ncvar_def("HAR_75pc", unit="g.m-2.d-1", longname = "Harvest - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var150 = ncvar_def("HAR_97pt5pc", unit="g.m-2.d-1", longname = "Harvest - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)

## Mean Residence Times
# Foliar
var61 = ncvar_def("MTT_fol_2pt5pc", unit="year", longname = "Mean Foliar Transit Time - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var62 = ncvar_def("MTT_fol_25pc", unit="year", longname = "Mean Foliar Transit Time - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var63 = ncvar_def("MTT_fol", unit="year", longname = "Mean Foliar Transit Time - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var64 = ncvar_def("MTT_fol_75pc", unit="year", longname = "Mean Foliar Transit Time - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var65 = ncvar_def("MTT_fol_97pt5pc", unit="year", longname = "Mean Foliar Transit Time - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# Fine root
var66 = ncvar_def("MTT_root_2pt5pc", unit="year", longname = "Mean fine root Transit Time - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var67 = ncvar_def("MTT_root_25pc", unit="year", longname = "Mean fine root Transit Time - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var68 = ncvar_def("MTT_root", unit="year", longname = "Mean fine root Transit Time - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var69 = ncvar_def("MTT_root_75pc", unit="year", longname = "Mean fine root Transit Time - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var70 = ncvar_def("MTT_root_97pt5pc", unit="year", longname = "Mean fine root Transit Time - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# Wood
var71 = ncvar_def("MTT_wood_2pt5pc", unit="year", longname = "Mean wood Transit Time - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var72 = ncvar_def("MTT_wood_25pc", unit="year", longname = "Mean wood Transit Time - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var73 = ncvar_def("MTT_wood", unit="year", longname = "Mean wood Transit Time - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var74 = ncvar_def("MTT_wood_75pc", unit="year", longname = "Mean wood Transit Time - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var75 = ncvar_def("MTT_wood_97pt5pc", unit="year", longname = "Mean wood Transit Time - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# Fine litter (fol + fine root)
var76 = ncvar_def("MTT_lit_2pt5pc", unit="year", longname = "Mean lit+litwood Transit Time - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var77 = ncvar_def("MTT_lit_25pc", unit="year", longname = "Mean lit+litwood Transit Time - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var78 = ncvar_def("MTT_lit", unit="year", longname = "Mean lit+litwood Transit Time - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var79 = ncvar_def("MTT_lit_75pc", unit="year", longname = "Mean lit+litwood Transit Time  - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var80 = ncvar_def("MTT_lit_97pt5pc", unit="year", longname = "Mean lit+litwood Transit Time  - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# Soil
var81 = ncvar_def("MTT_som_2pt5pc", unit="year", longname = "Mean Soil Transit Time - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var82 = ncvar_def("MTT_som_25pc", unit="year", longname = "Mean Soil Transit Time - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var83 = ncvar_def("MTT_som", unit="year", longname = "Mean Soil Transit Time - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var84 = ncvar_def("MTT_som_75pc", unit="year", longname = "Mean Soil Transit Time - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var85 = ncvar_def("MTT_som_97pt5pc", unit="year", longname = "Mean Soil Transit Time - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# Labile
var151 = ncvar_def("MTT_lab_2pt5pc", unit="year", longname = "Mean Labile Transit Time - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var152 = ncvar_def("MTT_lab_25pc", unit="year", longname = "Mean Labile Transit Time - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var153 = ncvar_def("MTT_lab", unit="year", longname = "Mean Labile Transit Time - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var154 = ncvar_def("MTT_lab_75pc", unit="year", longname = "Mean Labile Transit Time  - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var155 = ncvar_def("MTT_lab_97pt5pc", unit="year", longname = "Mean Labile Transit Time  - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# Soil carbon
var176 = ncvar_def("MTT_soilc_2pt5pc", unit="year", longname = "Mean Soil Transit Time - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var177 = ncvar_def("MTT_soilc_25pc", unit="year", longname = "Mean Soil Transit Time - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var178 = ncvar_def("MTT_soilc", unit="year", longname = "Mean Soil Transit Time - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var179 = ncvar_def("MTT_soilc_75pc", unit="year", longname = "Mean Soil Transit Time - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var180 = ncvar_def("MTT_soilc_97pt5pc", unit="year", longname = "Mean Soil Transit Time - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)

## OUTPUT fluxes
# Foliar
var86 = ncvar_def("OUTPUT_fol_flx_2pt5pc", unit="g.m-2.d-1", longname = "Carbon Loss from foliage - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var87 = ncvar_def("OUTPUT_fol_flx_25pc", unit="g.m-2.d-1", longname = "Carbon Loss from foliage - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var88 = ncvar_def("OUTPUT_fol_flx", unit="g.m-2.d-1", longname = "Carbon Loss from foliage - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var89 = ncvar_def("OUTPUT_fol_flx_75pc", unit="g.m-2.d-1", longname = "Carbon Loss from foliage - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var90 = ncvar_def("OUTPUT_fol_flx_97pt5pc", unit="g.m-2.d-1", longname = "Carbon Loss from foliage - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# Fine root
var91 = ncvar_def("OUTPUT_root_flx_2pt5pc", unit="g.m-2.d-1", longname = "Carbon Loss from fine root - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var92 = ncvar_def("OUTPUT_root_flx_25pc", unit="g.m-2.d-1", longname = "Carbon Loss from fine root - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var93 = ncvar_def("OUTPUT_root_flx", unit="g.m-2.d-1", longname = "Carbon Loss from fine root - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var94 = ncvar_def("OUTPUT_root_flx_75pc", unit="g.m-2.d-1", longname = "Carbon Loss from fine root - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var95 = ncvar_def("OUTPUT_root_flx_97pt5pc", unit="g.m-2.d-1", longname = "Carbon Loss from fine root - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# Wood
var96 = ncvar_def("OUTPUT_wood_flx_2pt5pc", unit="g.m-2.d-1", longname = "Carbon Loss from wood - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var97 = ncvar_def("OUTPUT_wood_flx_25pc", unit="g.m-2.d-1", longname = "Carbon Loss from wood - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var98 = ncvar_def("OUTPUT_wood_flx", unit="g.m-2.d-1", longname = "Carbon Loss from wood - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var99 = ncvar_def("OUTPUT_wood_flx_75pc", unit="g.m-2.d-1", longname = "Carbon Loss from wood - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var100 = ncvar_def("OUTPUT_wood_flx_97pt5pc", unit="g.m-2.d-1", longname = "Carbon Loss from wood - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# Labile
var156 = ncvar_def("OUTPUT_lab_flx_2pt5pc", unit="g.m-2.d-1", longname = "Carbon Loss from labile - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var157 = ncvar_def("OUTPUT_lab_flx_25pc", unit="g.m-2.d-1", longname = "Carbon Loss from labile - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var158 = ncvar_def("OUTPUT_lab_flx", unit="g.m-2.d-1", longname = "Carbon Loss from labile - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var159 = ncvar_def("OUTPUT_lab_flx_75pc", unit="g.m-2.d-1", longname = "Carbon Loss from labile - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var160 = ncvar_def("OUTPUT_lab_flx_97pt5pc", unit="g.m-2.d-1", longname = "Carbon Loss from labile - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# Litter
var161 = ncvar_def("OUTPUT_lit_flx_2pt5pc", unit="g.m-2.d-1", longname = "Carbon Loss from litter - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var162 = ncvar_def("OUTPUT_lit_flx_25pc", unit="g.m-2.d-1", longname = "Carbon Loss from litter - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var163 = ncvar_def("OUTPUT_lit_flx", unit="g.m-2.d-1", longname = "Carbon Loss from litter - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var164 = ncvar_def("OUTPUT_lit_flx_75pc", unit="g.m-2.d-1", longname = "Carbon Loss from litter - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var165 = ncvar_def("OUTPUT_lit_flx_97pt5pc", unit="g.m-2.d-1", longname = "Carbon Loss from litter - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# Soil Organic Matter
var166 = ncvar_def("OUTPUT_soil_flx_2pt5pc", unit="g.m-2.d-1", longname = "Carbon Loss from soil organic matter - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var167 = ncvar_def("OUTPUT_soil_flx_25pc", unit="g.m-2.d-1", longname = "Carbon Loss from soil organic matter - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var168 = ncvar_def("OUTPUT_soil_flx", unit="g.m-2.d-1", longname = "Carbon Loss from soil organic matter - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var169 = ncvar_def("OUTPUT_soil_flx_75pc", unit="g.m-2.d-1", longname = "Carbon Loss from soil organic matter - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var170 = ncvar_def("OUTPUT_soil_flx_97pt5pc", unit="g.m-2.d-1", longname = "Carbon Loss from soil organic matter - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# Soil Carbon
var171 = ncvar_def("OUTPUT_soil_c_flx_2pt5pc", unit="g.m-2.d-1", longname = "Soil Carbon Loss- 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var172 = ncvar_def("OUTPUT_soil_c_flx_25pc", unit="g.m-2.d-1", longname = "Soil Carbon Loss - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var173 = ncvar_def("OUTPUT_soil_c_flx", unit="g.m-2.d-1", longname = "Soil Carbon Loss - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var174 = ncvar_def("OUTPUT_soil_c_flx_75pc", unit="g.m-2.d-1", longname = "Soil Carbon Loss - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var175 = ncvar_def("OUTPUT_soil_c_flx_97pt5pc", unit="g.m-2.d-1", longname = "Soil Carbon Loss - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# Total output
var181 = ncvar_def("OUTPUT_flx_2pt5pc", unit="g.m-2.d-1", longname = "Total Carbon Loss- 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var182 = ncvar_def("OUTPUT_flx_25pc", unit="g.m-2.d-1", longname = "Total Carbon Loss - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var183 = ncvar_def("OUTPUT_flx", unit="g.m-2.d-1", longname = "Total Carbon Loss - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var184 = ncvar_def("OUTPUT_flx_75pc", unit="g.m-2.d-1", longname = "Total Carbon Loss - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var185 = ncvar_def("OUTPUT_flx_97pt5pc", unit="g.m-2.d-1", longname = "Total Carbon Loss - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)

## NPP allocation fluxes
# Foliar
var101 = ncvar_def("NPP_fol_flx_2pt5pc", unit="g.m-2.d-1", longname = "Net Primary Productivity to foliage - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var102 = ncvar_def("NPP_fol_flx_25pc", unit="g.m-2.d-1", longname = "Net Primary Productivity to foliage - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var103 = ncvar_def("NPP_fol_flx", unit="g.m-2.d-1", longname = "Net Primary Productivity to foliage - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var104 = ncvar_def("NPP_fol_flx_75pc", unit="g.m-2.d-1", longname = "Net Primary Productivity to foliage - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var105 = ncvar_def("NPP_fol_flx_97pt5pc", unit="g.m-2.d-1", longname = "Net Primary Productivity to foliage - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# Fine root
var106 = ncvar_def("NPP_root_flx_2pt5pc", unit="g.m-2.d-1", longname = "Net Primary Productivity to fine root - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var107 = ncvar_def("NPP_root_flx_25pc", unit="g.m-2.d-1", longname = "Net Primary Productivity to fine root - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var108 = ncvar_def("NPP_root_flx", unit="g.m-2.d-1", longname = "Net Primary Productivity to fine root - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var109 = ncvar_def("NPP_root_flx_75pc", unit="g.m-2.d-1", longname = "Net Primary Productivity to fine root - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var110 = ncvar_def("NPP_root_flx_97pt5pc", unit="g.m-2.d-1", longname = "Net Primary Productivity to fine root - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# Wood
var111 = ncvar_def("NPP_wood_flx_2pt5pc", unit="g.m-2.d-1", longname = "Net Primary Productivity to wood - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var112 = ncvar_def("NPP_wood_flx_25pc", unit="g.m-2.d-1", longname = "Net Primary Productivity to wood - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var113 = ncvar_def("NPP_wood_flx", unit="g.m-2.d-1", longname = "Net Primary Productivity to wood - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var114 = ncvar_def("NPP_wood_flx_75pc", unit="g.m-2.d-1", longname = "Net Primary Productivity to wood - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var115 = ncvar_def("NPP_wood_flx_97pt5pc", unit="g.m-2.d-1", longname = "Net Primary Productivity to wood - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# Labile
var116 = ncvar_def("NPP_lab_flx_2pt5pc", unit="g.m-2.d-1", longname = "Net Primary Productivity to labile - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var117 = ncvar_def("NPP_lab_flx_25pc", unit="g.m-2.d-1", longname = "Net Primary Productivity to labile - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var118 = ncvar_def("NPP_lab_flx", unit="g.m-2.d-1", longname = "Net Primary Productivity to labile - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var119 = ncvar_def("NPP_lab_flx_75pc", unit="g.m-2.d-1", longname = "Net Primary Productivity to labile - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
var120 = ncvar_def("NPP_lab_flx_97pt5pc", unit="g.m-2.d-1", longname = "Net Primary Productivity to labile - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)

# ## C economy
# # CUE
# var116 = ncvar_def("CUE_2pt5pc", unit="1", longname = "Carbon Use Efficiency - 2.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# var117 = ncvar_def("CUE_25pc", unit="1", longname = "Carbon Use Efficiency - 25 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# var118 = ncvar_def("CUE", unit="1", longname = "Carbon Use Efficiency - 50 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# var119 = ncvar_def("CUE_75pc", unit="1", longname = "Carbon Use Efficiency - 75 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
# var120 = ncvar_def("CUE_97pt5pc", unit="1", longname = "Carbon Use Efficiency - 97.5 percentile", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)


# create the empty file
new_file=nc_create(filename=output_name, vars=list(var0,var1,var2,var3,var4,var5,var6,var7,var8,var9,var10,
                                                   var11,var12,var13,var14,var15,var16,var17,var18,var19,var20,
                                                   var21,var22,var23,var24,var25,var26,var27,var28,var29,var30,
                                                   var31,var32,var33,var34,var35,var36,var37,var38,var39,var40,
                                                   var41,var42,var43,var44,var45,var46,var47,var48,var49,var50,
                                                   var51,var52,var53,var54,var55,var56,var57,var58,var59,var60,
                                                   var61,var62,var63,var64,var65,var66,var67,var68,var69,var70,
                                                   var71,var72,var73,var74,var75,var76,var77,var78,var79,var80,
                                                   var81,var82,var83,var84,var85,var86,var87,var88,var89,var90,
                                                   var91,var92,var93,var94,var95,var96,var97,var98,var99,var100,
                                                   var101,var102,var103,var104,var105,var106,var107,var108,var109,var110,
                                                   var111,var112,var113,var114,var115,var116,var117,var118,var119,var120,
                                                   var121,var122,var123,var124,var125,var126,var127,var128,var129,var130,
                                                   var131,var132,var133,var134,var135,var136,var137,var138,var139,var140,
                                                   var141,var142,var143,var144,var145,var146,var147,var148,var149,var150,
                                                   var151,var152,var153,var154,var155,var156,var157,var158,var159,var160,
                                                   var161,var162,var163,var164,var165,var166,var167,var168,var169,var170,
                                                   var171,var172,var173,var174,var175,var176,var177,var178,var179,var180,
                                                   var181,var182,var183,var184,var185), force_v4 = TRUE)
                                                   
###
# Load data into output variable
###

## TIMING
ncvar_put(new_file, var0, (cumsum(PROJECT$model$timestep_days)-PROJECT$model$timestep_days[1]+1))

## STATES
# LAI
ncvar_put(new_file, var1, LAI[,,1,])
ncvar_put(new_file, var2, LAI[,,3,])
ncvar_put(new_file, var3, LAI[,,4,])
ncvar_put(new_file, var4, LAI[,,5,])
ncvar_put(new_file, var5, LAI[,,7,])
# LAB
ncvar_put(new_file, var6, LAB[,,1,])
ncvar_put(new_file, var7, LAB[,,3,])
ncvar_put(new_file, var8, LAB[,,4,])
ncvar_put(new_file, var9, LAB[,,5,])
ncvar_put(new_file, var10, LAB[,,7,])
# FOL
ncvar_put(new_file, var11, FOL[,,1,])
ncvar_put(new_file, var12, FOL[,,3,])
ncvar_put(new_file, var13, FOL[,,4,])
ncvar_put(new_file, var14, FOL[,,5,])
ncvar_put(new_file, var15, FOL[,,7,])
# ROOT
ncvar_put(new_file, var16, ROOT[,,1,])
ncvar_put(new_file, var17, ROOT[,,3,])
ncvar_put(new_file, var18, ROOT[,,4,])
ncvar_put(new_file, var19, ROOT[,,5,])
ncvar_put(new_file, var20, ROOT[,,7,])
# WOOD
ncvar_put(new_file, var21, WOOD[,,1,])
ncvar_put(new_file, var22, WOOD[,,3,])
ncvar_put(new_file, var23, WOOD[,,4,])
ncvar_put(new_file, var24, WOOD[,,5,])
ncvar_put(new_file, var25, WOOD[,,7,])
# LIT
ncvar_put(new_file, var26, LIT[,,1,])
ncvar_put(new_file, var27, LIT[,,3,])
ncvar_put(new_file, var28, LIT[,,4,])
ncvar_put(new_file, var29, LIT[,,5,])
ncvar_put(new_file, var30, LIT[,,7,])
# SOIL
ncvar_put(new_file, var31, SOIL[,,1,])
ncvar_put(new_file, var32, SOIL[,,3,])
ncvar_put(new_file, var33, SOIL[,,4,])
ncvar_put(new_file, var34, SOIL[,,5,])
ncvar_put(new_file, var35, SOIL[,,7,])
# TOT
ncvar_put(new_file, var136, TOT[,,1,])
ncvar_put(new_file, var137, TOT[,,3,])
ncvar_put(new_file, var138, TOT[,,4,])
ncvar_put(new_file, var139, TOT[,,5,])
ncvar_put(new_file, var140, TOT[,,7,])
# SOIL_C
ncvar_put(new_file, var141, SOILC[,,1,])
ncvar_put(new_file, var142, SOILC[,,3,])
ncvar_put(new_file, var143, SOILC[,,4,])
ncvar_put(new_file, var144, SOILC[,,5,])
ncvar_put(new_file, var145, SOILC[,,7,])

## FLUXES
# GPP
ncvar_put(new_file, var36, GPP[,,1,])
ncvar_put(new_file, var37, GPP[,,3,])
ncvar_put(new_file, var38, GPP[,,4,])
ncvar_put(new_file, var39, GPP[,,5,])
ncvar_put(new_file, var40, GPP[,,7,])
# Ra
ncvar_put(new_file, var41, Ra[,,1,])
ncvar_put(new_file, var42, Ra[,,3,])
ncvar_put(new_file, var43, Ra[,,4,])
ncvar_put(new_file, var44, Ra[,,5,])
ncvar_put(new_file, var45, Ra[,,7,])
# Rh
ncvar_put(new_file, var46, Rh[,,1,])
ncvar_put(new_file, var47, Rh[,,3,])
ncvar_put(new_file, var48, Rh[,,4,])
ncvar_put(new_file, var49, Rh[,,5,])
ncvar_put(new_file, var50, Rh[,,7,])
# NPP
ncvar_put(new_file, var51, NPP[,,1,])
ncvar_put(new_file, var52, NPP[,,3,])
ncvar_put(new_file, var53, NPP[,,4,])
ncvar_put(new_file, var54, NPP[,,5,])
ncvar_put(new_file, var55, NPP[,,7,])
# FIR
ncvar_put(new_file, var56, FIR[,,1,])
ncvar_put(new_file, var57, FIR[,,3,])
ncvar_put(new_file, var58, FIR[,,4,])
ncvar_put(new_file, var59, FIR[,,5,])
ncvar_put(new_file, var60, FIR[,,7,])
# RECO
ncvar_put(new_file, var121, RECO[,,1,])
ncvar_put(new_file, var122, RECO[,,3,])
ncvar_put(new_file, var123, RECO[,,4,])
ncvar_put(new_file, var124, RECO[,,5,])
ncvar_put(new_file, var125, RECO[,,7,])
# NEE
ncvar_put(new_file, var126, NEE[,,1,])
ncvar_put(new_file, var127, NEE[,,3,])
ncvar_put(new_file, var128, NEE[,,4,])
ncvar_put(new_file, var129, NEE[,,5,])
ncvar_put(new_file, var130, NEE[,,7,])
# NBE
ncvar_put(new_file, var131, NBE[,,1,])
ncvar_put(new_file, var132, NBE[,,3,])
ncvar_put(new_file, var133, NBE[,,4,])
ncvar_put(new_file, var134, NBE[,,5,])
ncvar_put(new_file, var135, NBE[,,7,])
# HAR
ncvar_put(new_file, var146, HAR[,,1,])
ncvar_put(new_file, var147, HAR[,,3,])
ncvar_put(new_file, var148, HAR[,,4,])
ncvar_put(new_file, var149, HAR[,,5,])
ncvar_put(new_file, var150, HAR[,,7,])

## MTT - time series
# FOL
ncvar_put(new_file, var61, fMRT[,,1,])
ncvar_put(new_file, var62, fMRT[,,3,])
ncvar_put(new_file, var63, fMRT[,,4,])
ncvar_put(new_file, var64, fMRT[,,5,])
ncvar_put(new_file, var65, fMRT[,,7,])
# ROOT
ncvar_put(new_file, var66, rMRT[,,1,])
ncvar_put(new_file, var67, rMRT[,,3,])
ncvar_put(new_file, var68, rMRT[,,4,])
ncvar_put(new_file, var69, rMRT[,,5,])
ncvar_put(new_file, var70, rMRT[,,7,])
# WOOD
ncvar_put(new_file, var71, wMRT[,,1,])
ncvar_put(new_file, var72, wMRT[,,3,])
ncvar_put(new_file, var73, wMRT[,,4,])
ncvar_put(new_file, var74, wMRT[,,5,])
ncvar_put(new_file, var75, wMRT[,,7,])
# LIT
ncvar_put(new_file, var76, liMRT[,,1,])
ncvar_put(new_file, var77, liMRT[,,3,])
ncvar_put(new_file, var78, liMRT[,,4,])
ncvar_put(new_file, var79, liMRT[,,5,])
ncvar_put(new_file, var80, liMRT[,,7,])
# SOIL
ncvar_put(new_file, var81, sMRT[,,1,])
ncvar_put(new_file, var82, sMRT[,,3,])
ncvar_put(new_file, var83, sMRT[,,4,])
ncvar_put(new_file, var84, sMRT[,,5,])
ncvar_put(new_file, var85, sMRT[,,7,])
# LABILE
ncvar_put(new_file, var151, laMRT[,,1,])
ncvar_put(new_file, var152, laMRT[,,3,])
ncvar_put(new_file, var153, laMRT[,,4,])
ncvar_put(new_file, var154, laMRT[,,5,])
ncvar_put(new_file, var155, laMRT[,,7,])
# SOIL
ncvar_put(new_file, var176, scMRT[,,1,])
ncvar_put(new_file, var177, scMRT[,,3,])
ncvar_put(new_file, var178, scMRT[,,4,])
ncvar_put(new_file, var179, scMRT[,,5,])
ncvar_put(new_file, var180, scMRT[,,7,])

## OUTPUT fluxes
# FOL
ncvar_put(new_file, var86, fOUTPUT_FLX[,,1,])
ncvar_put(new_file, var87, fOUTPUT_FLX[,,3,])
ncvar_put(new_file, var88, fOUTPUT_FLX[,,4,])
ncvar_put(new_file, var89, fOUTPUT_FLX[,,5,])
ncvar_put(new_file, var90, fOUTPUT_FLX[,,7,])
# ROOT
ncvar_put(new_file, var91, rOUTPUT_FLX[,,1,])
ncvar_put(new_file, var92, rOUTPUT_FLX[,,3,])
ncvar_put(new_file, var93, rOUTPUT_FLX[,,4,])
ncvar_put(new_file, var94, rOUTPUT_FLX[,,5,])
ncvar_put(new_file, var95, rOUTPUT_FLX[,,7,])
# WOOD
ncvar_put(new_file, var96, wOUTPUT_FLX[,,1,])
ncvar_put(new_file, var97, wOUTPUT_FLX[,,3,])
ncvar_put(new_file, var98, wOUTPUT_FLX[,,4,])
ncvar_put(new_file, var99, wOUTPUT_FLX[,,5,])
ncvar_put(new_file, var100, wOUTPUT_FLX[,,7,])
# LABILE
ncvar_put(new_file, var156, laOUTPUT_FLX[,,1,])
ncvar_put(new_file, var157, laOUTPUT_FLX[,,3,])
ncvar_put(new_file, var158, laOUTPUT_FLX[,,4,])
ncvar_put(new_file, var159, laOUTPUT_FLX[,,5,])
ncvar_put(new_file, var160, laOUTPUT_FLX[,,7,])
# LITTER
ncvar_put(new_file, var161, liOUTPUT_FLX[,,1,])
ncvar_put(new_file, var162, liOUTPUT_FLX[,,3,])
ncvar_put(new_file, var163, liOUTPUT_FLX[,,4,])
ncvar_put(new_file, var164, liOUTPUT_FLX[,,5,])
ncvar_put(new_file, var165, liOUTPUT_FLX[,,7,])
# SOIL ORGANIC MATTER
ncvar_put(new_file, var166, sOUTPUT_FLX[,,1,])
ncvar_put(new_file, var167, sOUTPUT_FLX[,,3,])
ncvar_put(new_file, var168, sOUTPUT_FLX[,,4,])
ncvar_put(new_file, var169, sOUTPUT_FLX[,,5,])
ncvar_put(new_file, var170, sOUTPUT_FLX[,,7,])
# SOIL CARBON
ncvar_put(new_file, var171, scOUTPUT_FLX[,,1,])
ncvar_put(new_file, var172, scOUTPUT_FLX[,,3,])
ncvar_put(new_file, var173, scOUTPUT_FLX[,,4,])
ncvar_put(new_file, var174, scOUTPUT_FLX[,,5,])
ncvar_put(new_file, var175, scOUTPUT_FLX[,,7,])
# Total
ncvar_put(new_file, var181, OUTPUT_FLX[,,1,])
ncvar_put(new_file, var182, OUTPUT_FLX[,,3,])
ncvar_put(new_file, var183, OUTPUT_FLX[,,4,])
ncvar_put(new_file, var184, OUTPUT_FLX[,,5,])
ncvar_put(new_file, var185, OUTPUT_FLX[,,7,])

## NPP fluxes
# FOL
ncvar_put(new_file, var101, fNPP_FLX[,,1,])
ncvar_put(new_file, var102, fNPP_FLX[,,3,])
ncvar_put(new_file, var103, fNPP_FLX[,,4,])
ncvar_put(new_file, var104, fNPP_FLX[,,5,])
ncvar_put(new_file, var105, fNPP_FLX[,,7,])
# ROOT
ncvar_put(new_file, var106, rNPP_FLX[,,1,])
ncvar_put(new_file, var107, rNPP_FLX[,,3,])
ncvar_put(new_file, var108, rNPP_FLX[,,4,])
ncvar_put(new_file, var109, rNPP_FLX[,,5,])
ncvar_put(new_file, var110, rNPP_FLX[,,7,])
# WOOD
ncvar_put(new_file, var111, wNPP_FLX[,,1,])
ncvar_put(new_file, var112, wNPP_FLX[,,3,])
ncvar_put(new_file, var113, wNPP_FLX[,,4,])
ncvar_put(new_file, var114, wNPP_FLX[,,5,])
ncvar_put(new_file, var115, wNPP_FLX[,,7,])
# LABILE
ncvar_put(new_file, var116, laNPP_FLX[,,1,])
ncvar_put(new_file, var117, laNPP_FLX[,,3,])
ncvar_put(new_file, var118, laNPP_FLX[,,4,])
ncvar_put(new_file, var119, laNPP_FLX[,,5,])
ncvar_put(new_file, var120, laNPP_FLX[,,7,])

# ## Carbon economy
# # CUE
# ncvar_put(new_file, var116, CUE[,,1,])
# ncvar_put(new_file, var117, CUE[,,3,])
# ncvar_put(new_file, var118, CUE[,,4,])
# ncvar_put(new_file, var119, CUE[,,5,])
# ncvar_put(new_file, var120, CUE[,,7,])

## close the file to write to disk
nc_close(new_file)
setwd("/home/cnwobi/CARDAMOM/CARDAMOM_DEV/CARDAMOM/cssp_rainfor_brazil")
