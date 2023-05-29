
###
## Script for the creation of a latex script to generate C-budgets for DALEC
## Author: T L Smallman (t.l.smallman@ed.ac.uk)
## Created: 12/07/2022
## Last updated: 29/05/2023
## Updated by: C J Nwobi (nwobicj@gmail.com) 
## NOTES: 
## 1) This version presents a combined labile and foliage pool. 
## 2) Assumes a CDEA style phenology with both direct GPP to foliage and via labile pathways
## 3) This includes CUE and Mean residence time (MRT) rather than Meam Transient time (MTT)
## 4) Only use if Fire and Biomass removals have been excluded from CARDAMOM-DALEC runs
## 5) Creates different latex files based on number of sub-regions
###

###
## Prepare the work space

# Set working directory
setwd("/home/cnwobi/CARDAMOM/CARDAMOM/")

# Load R libraries
library(rgdal);library(raster)

# Load user defined functions
source("./R_functions/load_all_cardamom_functions.r")

###
## Load files of entire region over which CARDAMOM ran and sub regions interested in (found in shape_files_for_analysis)  
amazon_nw_poly <- shapefile("/exports/csce/datastore/geos/groups/gcel/cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/amazon_nw.shp")
amazon_sw_poly <- shapefile("/exports/csce/datastore/geos/groups/gcel/cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/amazon_sw.shp")
amazon_ec_poly <- shapefile("/exports/csce/datastore/geos/groups/gcel/cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/amazon_ec.shp")
amazon_bs_poly <- shapefile("/exports/csce/datastore/geos/groups/gcel/cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/amazon_bs.shp")
amazon_gs_poly <- shapefile("/exports/csce/datastore/geos/groups/gcel/cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/amazon_gs.shp")

#load raster files of extrapolated RAINFOR data (contact university of leeds for data David Galbraith- D.R.Galbraith@leeds.ac.uk)
biomass_amazon <- brick('/exports/csce/datastore/geos/groups/gcel/cssp_rainfor_amazon_brazil/rainfor_leeds_data/AbovegroundBiomass_Mg_perHa_111km.tif')
biomass_amazon_C_ifl<-brick("/exports/csce/datastore/geos/groups/gcel/cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/biomass_ifl_subset.tif")

#create shape file from RAINFOR data
biomass_amazon_mask <- biomass_amazon > -Inf
biomass_amazon_mask_pol <- rasterToPolygons(biomass_amazon_mask, dissolve=TRUE)

#re-create raster of amazonian region and assign value
amazon_region_new <- extent(biomass_amazon_C_ifl)
r_amazon_region_new <- raster(amazon_region_new,res=1)
values(r_amazon_region_new) <- 0
r_amazon_region_final <- rasterize(biomass_amazon_mask_pol,r_amazon_region_new)

#function to subset amazonia region based on ecoregions shape files
#a<-amazonian extent, s<-amazonia ecoregion
country_to_amazon_crop_fun <- function(a,s) { 
  cropped <- crop(a, extent(s))
  masked1 <- mask(cropped, s)
  masked2 <- masked1 > -Inf
  poly <- rasterToPolygons(masked2, dissolve=TRUE)
}

amazon_nw_pol_mask <- country_to_amazon_crop_fun(biomass_amazon_C_ifl,amazon_nw_poly)
amazon_sw_pol_mask <- country_to_amazon_crop_fun(biomass_amazon_C_ifl,amazon_sw_poly)
amazon_ec_pol_mask <- country_to_amazon_crop_fun(biomass_amazon_C_ifl,amazon_ec_poly)
amazon_bs_pol_mask <- country_to_amazon_crop_fun(biomass_amazon_C_ifl,amazon_bs_poly)
amazon_gs_pol_mask <- country_to_amazon_crop_fun(biomass_amazon_C_ifl,amazon_gs_poly)

#create raster of amazonia and ecoregions with values assigned to each region
r_amazonia_nw <- mask(r_amazon_region_final, amazon_nw_pol_mask, updatevalue=0)
r_amazonia_sw <- mask(r_amazonia_nw, amazon_sw_pol_mask, inverse=T, updatevalue=2)
r_amazonia_bs <- mask(r_amazonia_sw, amazon_bs_pol_mask, inverse=T, updatevalue=3)
r_amazonia_ec <- mask(r_amazonia_bs, amazon_ec_pol_mask, inverse=T, updatevalue=4)
r_amazonia_ecoregions <- mask(r_amazonia_ec, amazon_gs_pol_mask, inverse=T, updatevalue=5)

#create matix to match matrix of CARDAMOM outputs and quantiles wanted
amazonia_ecoregions_matrix <- t(as.matrix(flip(r_amazonia_ecoregions,2)))

amazonia_ecoregions_matrices <- array(numeric(),c(37,32,3)) 
amazonia_ecoregions_matrices[,,1]<-amazonia_ecoregions_matrix
amazonia_ecoregions_matrices[,,2]<-amazonia_ecoregions_matrix
amazonia_ecoregions_matrices[,,3]<-amazonia_ecoregions_matrix

ecoregion_names <- c('northwest','southwest','brazilshield','eastcentral','guyanashield') #ecoregion names

###loop run for each ecoregion
for (e in ecoregion_names) {
  if (e =='northwest') {
    region_mask <- amazonia_ecoregions_matrices==1
  }
  else if (e =='southwest') {
    region_mask <- amazonia_ecoregions_matrices==2
  }
  else if (e =='brazilshield') {
    region_mask <- amazonia_ecoregions_matrices==3
  }
  else if (e =='eastcentral') {
    region_mask <- amazonia_ecoregions_matrices==4
  }
  else {
    region_mask <- amazonia_ecoregions_matrices==5
  }
  
###
## Load files from which C-budget is extracted, prepare C-budget values
  # Load information file found in amazonia_ifl_cardamom_runs zip folder
# load("/home/cnwobi/CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/amazonia_ifl_esa_cci_agb_nomngt/infofile.RData")
# load("/home/cnwobi/CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/amazonia_ifl_rainfor_biomass_productivity_2005_nomngt/infofile.RData")
# load("/home/cnwobi/CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/amazonia_ifl_rainfor_biomass_annual_productivity_nomngt/infofile.RData")
load("/home/cnwobi/CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/amazonia_ifl_rainfor_biomass_annual_nomngt/infofile.RData")
  
# Define the output file name for the created script
outfilename = paste(e,"_tha_wcue_wMRT_",PROJECT$name,"_latex_C_budget.tex",sep="")

#
# If this is a site analysis 
# which site are we using?
site_nos = 1

###
## Define figure labels

# The caption must be written in correct latex syntex
# NOTE: if the latex language requires use of "\" ensure it is a "\\".
figure_caption = paste("Steady state C-budget for assimilation of ",e," ",chartr("_"," ",PROJECT$name),". Numbers show median estimate of fluxes (alongside arrows) and of stocks (in boxes). Units are MgC ha^\\({-1}\\) for stocks and MgC ha\\(^{-1}\\) yr^\\({-1}\\) for fluxes. 95\\% confidence intervals are shown in a fractional form with 2.5 and 97.5 percentiles as numerator and denominator. Black fluxes are biogenic, including net primary production (\\(NPP\\)), mortality (\\(Mort\\)), autotrophic respiration (\\(Ra\\)) and heterotrophic respiration (\\(Rh\\)). \\(NEE = Ra + Rh - GPP\\). \\(NBE = NEE + E_{total}\\). However, \\(NEE = NBE\\) as there is no fire emission.",sep="")
# The label will be used for referencing the figure in the latex document
figure_label = "SIFig:fcp_budget"
# Desired precision, i.e. decimal places
dp = 1

if (PROJECT$spatial_type == "grid") {
    # A gridded analysis
    
    # Specifiy quantiles we want to extract, should only be 3 (low/median/high).
    # From gridded analysis this must be from thos available in the file. 
    # Check grid_output$num_quantiles for available
    quantiles_wanted = c(1,4,7)
    
    # Load the processed site file
    load(paste(PROJECT$results_processedpath,PROJECT$name,"_stock_flux.RData",sep=""))

    #Masking phase
    # NATURAL FLUXES
    grid_output$mean_gpp_gCm2day[,,quantiles_wanted][!region_mask]<-NA
    grid_output$mean_cue[,,quantiles_wanted][!region_mask]<-NA
    grid_output$mean_rauto_gCm2day[,,quantiles_wanted][!region_mask]<-NA
    grid_output$mean_rhet_litter_gCm2day[,,quantiles_wanted][!region_mask]<-NA
    grid_output$mean_rhet_som_gCm2day[,,quantiles_wanted][!region_mask]<-NA
    grid_output$mean_rhet_gCm2day[,,quantiles_wanted][!region_mask]<-NA
    grid_output$mean_npp_gCm2day[,,quantiles_wanted][!region_mask]<-NA
    grid_output$mean_combined_alloc_foliage_gCm2day[,,quantiles_wanted][!region_mask]<-NA
    grid_output$mean_alloc_roots_gCm2day[,,quantiles_wanted][!region_mask]<-NA
    grid_output$mean_alloc_wood_gCm2day[,,quantiles_wanted][!region_mask]<-NA
    grid_output$mean_foliage_to_litter_gCm2day[,,quantiles_wanted][!region_mask]<-NA
    grid_output$mean_roots_to_litter_gCm2day[,,quantiles_wanted][!region_mask]<-NA
    grid_output$mean_wood_to_litter_gCm2day[,,quantiles_wanted][!region_mask]<-NA
    grid_output$mean_litter_to_som_gCm2day[,,quantiles_wanted][!region_mask]<-NA
    # FIRE FLUXES
    grid_output$mean_FIRElitter_labile_gCm2day[,,quantiles_wanted][!region_mask]<-NA
    grid_output$mean_FIRElitter_foliage_gCm2day[,,quantiles_wanted][!region_mask]<-NA
    grid_output$mean_FIRElitter_roots_gCm2day[,,quantiles_wanted][!region_mask]<-NA
    grid_output$mean_FIRElitter_wood_gCm2day[,,quantiles_wanted][!region_mask]<-NA
    grid_output$mean_FIRElitter_litter_gCm2day[,,quantiles_wanted][!region_mask]<-NA
    grid_output$mean_FIREemiss_labile_gCm2day[,,quantiles_wanted] [!region_mask]<-NA
    grid_output$mean_FIREemiss_foliage_gCm2day[,,quantiles_wanted][!region_mask]<-NA
    grid_output$mean_FIREemiss_roots_gCm2day[,,quantiles_wanted][!region_mask]<-NA
    grid_output$mean_FIREemiss_wood_gCm2day[,,quantiles_wanted][!region_mask]<-NA
    grid_output$mean_FIREemiss_litter_gCm2day[,,quantiles_wanted][!region_mask]<-NA
    grid_output$mean_FIREemiss_som_gCm2day[,,quantiles_wanted][!region_mask]<-NA
    grid_output$mean_fire_gCm2day[,,quantiles_wanted][!region_mask]<-NA
    # Net fluxes
    grid_output$mean_nbe_gCm2day[,,quantiles_wanted][!region_mask]<-NA
    grid_output$mean_nee_gCm2day[,,quantiles_wanted][!region_mask]<-NA
    # STOCKS
    grid_output$mean_labile_gCm2[,,quantiles_wanted] [!region_mask]<-NA
    grid_output$mean_foliage_gCm2[,,quantiles_wanted][!region_mask]<-NA
    grid_output$mean_roots_gCm2[,,quantiles_wanted][!region_mask]<-NA
    grid_output$mean_wood_gCm2[,,quantiles_wanted][!region_mask]<-NA
    grid_output$mean_litter_gCm2[,,quantiles_wanted][!region_mask]<-NA
    grid_output$mean_som_gCm2[,,quantiles_wanted][!region_mask]<-NA
    # MEAN TRANSIENT TIMES
    grid_output$MTT_wood_years[,,quantiles_wanted][!region_mask]<-NA
    grid_output$MTT_foliage_years[,,quantiles_wanted][!region_mask]<-NA
    grid_output$MTT_labile_years[,,quantiles_wanted][!region_mask]<-NA
    grid_output$MTT_roots_years[,,quantiles_wanted][!region_mask]<-NA
    grid_output$MTT_litter_years[,,quantiles_wanted][!region_mask]<-NA
    grid_output$MTT_som_years[,,quantiles_wanted][!region_mask]<-NA

    
    # Extract or calculate required derived values
    # NOTE: unit conversion from gC/m2/day to gC/m2/yr

    # NATURAL FLUXES
    gpp_gCm2yr = format(round(apply(grid_output$mean_gpp_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
    cue = format(round(apply(grid_output$mean_cue[,,quantiles_wanted],3,mean, na.rm=TRUE), digits = 2), nsmall = dp)
    rauto_gCm2yr = format(round(apply(grid_output$mean_rauto_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
    rhet_litter_gCm2yr = format(round(apply(grid_output$mean_rhet_litter_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
    rhet_som_gCm2yr = format(round(apply(grid_output$mean_rhet_som_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
    rhet_gCm2yr = format(round(apply(grid_output$mean_rhet_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
    npp_gCm2yr = format(round(apply(grid_output$mean_npp_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
    npp_labilefoliage_gCm2yr = format(round(apply(grid_output$mean_combined_alloc_foliage_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
    npp_roots_gCm2yr = format(round(apply(grid_output$mean_alloc_roots_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
    npp_wood_gCm2yr = format(round(apply(grid_output$mean_alloc_wood_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
    foliage_to_litter_gCm2yr = format(round(apply(grid_output$mean_foliage_to_litter_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
    roots_to_litter_gCm2yr = format(round(apply(grid_output$mean_roots_to_litter_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
    wood_to_litter_gCm2yr = format(round(apply(grid_output$mean_wood_to_litter_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
    litter_to_som_gCm2yr = format(round(apply(grid_output$mean_litter_to_som_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
    # FIRE FLUXES
    FIRElitter_labile_foliage_gCm2yr = format(round(apply(grid_output$mean_FIRElitter_labile_gCm2day[,,quantiles_wanted]+grid_output$mean_FIRElitter_foliage_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
    FIRElitter_roots_gCm2yr = format(round(apply(grid_output$mean_FIRElitter_roots_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
    FIRElitter_wood_gCm2yr = format(round(apply(grid_output$mean_FIRElitter_wood_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
    FIRElitter_litter_gCm2yr = format(round(apply(grid_output$mean_FIRElitter_litter_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
    FIREemiss_labile_foliage_gCm2yr = format(round(apply(grid_output$mean_FIREemiss_labile_gCm2day[,,quantiles_wanted]+grid_output$mean_FIREemiss_foliage_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
    FIREemiss_roots_gCm2yr = format(round(apply(grid_output$mean_FIREemiss_roots_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
    FIREemiss_wood_gCm2yr = format(round(apply(grid_output$mean_FIREemiss_wood_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
    FIREemiss_litter_gCm2yr = format(round(apply(grid_output$mean_FIREemiss_litter_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
    FIREemiss_som_gCm2yr = format(round(apply(grid_output$mean_FIREemiss_som_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
    fire_gCm2yr = format(round(apply(grid_output$mean_fire_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
    # Net fluxes
    nbe_gCm2yr = format(round(apply(grid_output$mean_nbe_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
    nee_gCm2yr = format(round(apply(grid_output$mean_nee_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
    # STOCKS
    labile_foliage_gCm2 = format(round(apply(grid_output$mean_labile_gCm2[,,quantiles_wanted]+grid_output$mean_foliage_gCm2[,,quantiles_wanted],3,mean, na.rm=TRUE) /100, digits = dp), nsmall = dp)
    roots_gCm2 = format(round(apply(grid_output$mean_roots_gCm2[,,quantiles_wanted],3,mean, na.rm=TRUE) /100, digits = dp), nsmall = dp)
    wood_gCm2 = format(round(apply(grid_output$mean_wood_gCm2[,,quantiles_wanted],3,mean, na.rm=TRUE) /100, digits = dp), nsmall = dp)
    litter_gCm2 = format(round(apply(grid_output$mean_litter_gCm2[,,quantiles_wanted],3,mean, na.rm=TRUE) /100, digits = dp), nsmall = dp)
    som_gCm2 = format(round(apply(grid_output$mean_som_gCm2[,,quantiles_wanted],3,mean, na.rm=TRUE) /100, digits = dp), nsmall = dp)
    # MEAN TRANSIENT TIMES
    mtt_wood_years = format(round(apply(grid_output$MTT_wood_years[,,quantiles_wanted],3,mean, na.rm=TRUE), digits = dp), nsmall = dp)
    mtt_fol_years = format(round(apply(grid_output$MTT_foliage_years[,,quantiles_wanted],3,mean, na.rm=TRUE), digits = dp), nsmall = dp)
    mtt_lab_years = format(round(apply(grid_output$MTT_labile_years[,,quantiles_wanted],3,mean, na.rm=TRUE), digits = dp), nsmall = dp)
    mtt_root_years = format(round(apply(grid_output$MTT_roots_years[,,quantiles_wanted],3,mean, na.rm=TRUE), digits = dp), nsmall = dp)
    mtt_litter_years = format(round(apply(grid_output$MTT_litter_years[,,quantiles_wanted],3,mean, na.rm=TRUE), digits = dp), nsmall = dp)
    mtt_som_years = format(round(apply(grid_output$MTT_som_years[,,quantiles_wanted],3,mean, na.rm=TRUE), digits = dp), nsmall = dp)
    
    
}
else if (PROJECT$spatial_type == "site") {
  # A site analysis
  
  # Specifiy quantiles we want to extract, should only be 3 (low/median/high)
  quantiles_wanted = c(0.025,0.50,0.975)    
  
  # Load the processed site file
  load(paste(PROJECT$results_processedpath,PROJECT$sites[site_nos],".RData",sep=""))
  
  # Extract or calculate required derived values
  # NOTE: unit conversion from gC/m2/day to gC/m2/yr
  
  # NATURAL FLUXES
  gpp_gCm2yr = format(round(quantile(apply(states_all$gpp_gCm2day,1,mean), prob=quantiles_wanted) * (365.25/100), digits = dp), nsmall = dp)
  rauto_gCm2yr = format(round(quantile(apply(states_all$rauto_gCm2day,1,mean), prob=quantiles_wanted) * (365.25/100), digits = dp), nsmall = dp)
  rhet_litter_gCm2yr = format(round(quantile(apply(states_all$rhet_litter_gCm2day,1,mean), prob=quantiles_wanted) * (365.25/100), digits = dp), nsmall = dp)
  rhet_som_gCm2yr = format(round(quantile(apply(states_all$rhet_som_gCm2day,1,mean), prob=quantiles_wanted) * (365.25/100), digits = dp), nsmall = dp)
  rhet_gCm2yr = format(round(quantile(apply(states_all$rhet_litter_gCm2day+states_all$rhet_som_gCm2day,1,mean), prob=quantiles_wanted) * (365.25/100), digits = dp), nsmall = dp)
  npp_gCm2yr = format(round(quantile(apply(states_all$gpp_gCm2day-states_all$rauto_gCm2day,1,mean), prob=quantiles_wanted) * (365.25/100), digits = dp), nsmall = dp)
  npp_labilefoliage_gCm2yr = format(round(quantile(apply(states_all$alloc_foliage_gCm2day+states_all$labile_to_foliage_gCm2day,1,mean), prob=quantiles_wanted) * (365.25/100), digits = dp), nsmall = dp)
  npp_roots_gCm2yr = format(round(quantile(apply(states_all$alloc_roots_gCm2day,1,mean), prob=quantiles_wanted) * (365.25/100), digits = dp), nsmall = dp)
  npp_wood_gCm2yr = format(round(quantile(apply(states_all$alloc_wood_gCm2day,1,mean), prob=quantiles_wanted) * (365.25/100), digits = dp), nsmall = dp)
  foliage_to_litter_gCm2yr = format(round(quantile(apply(states_all$foliage_to_litter_gCm2day,1,mean), prob=quantiles_wanted) * (365.25/100), digits = dp), nsmall = dp)
  roots_to_litter_gCm2yr = format(round(quantile(apply(states_all$roots_to_litter_gCm2day,1,mean), prob=quantiles_wanted) * (365.25/100), digits = dp), nsmall = dp)
  wood_to_litter_gCm2yr = format(round(quantile(apply(states_all$wood_to_litter_gCm2day,1,mean), prob=quantiles_wanted) * (365.25/100), digits = dp), nsmall = dp)
  litter_to_som_gCm2yr = format(round(quantile(apply(states_all$litter_to_som_gCm2day,1,mean), prob=quantiles_wanted) * (365.25/100), digits = dp), nsmall = dp)
  # FIRE FLUXES
  FIRElitter_labile_foliage_gCm2yr = format(round(quantile(apply(states_all$FIRElitter_labile_gCm2day + states_all$FIRElitter_foliage_gCm2day,1,mean), prob=quantiles_wanted) * (365.25/100), digits = dp), nsmall = dp)
  FIRElitter_roots_gCm2yr = format(round(quantile(apply(states_all$FIRElitter_roots_gCm2day,1,mean), prob=quantiles_wanted) * (365.25/100), digits = dp), nsmall = dp)
  FIRElitter_wood_gCm2yr = format(round(quantile(apply(states_all$FIRElitter_wood_gCm2day,1,mean), prob=quantiles_wanted) * (365.25/100), digits = dp), nsmall = dp)
  FIRElitter_litter_gCm2yr = format(round(quantile(apply(states_all$FIRElitter_litter_gCm2day,1,mean), prob=quantiles_wanted) * (365.25/100), digits = dp), nsmall = dp)
  FIREemiss_labile_foliage_gCm2yr = format(round(quantile(apply(states_all$FIREemiss_labile_gCm2day + states_all$FIREemiss_foliage_gCm2day,1,mean), prob=quantiles_wanted) * (365.25/100), digits = dp), nsmall = dp)
  FIREemiss_roots_gCm2yr = format(round(quantile(apply(states_all$FIREemiss_roots_gCm2day,1,mean), prob=quantiles_wanted) * (365.25/100), digits = dp), nsmall = dp)
  FIREemiss_wood_gCm2yr = format(round(quantile(apply(states_all$FIREemiss_wood_gCm2day,1,mean), prob=quantiles_wanted) * (365.25/100), digits = dp), nsmall = dp)
  FIREemiss_litter_gCm2yr = format(round(quantile(apply(states_all$FIREemiss_litter_gCm2day,1,mean), prob=quantiles_wanted) * (365.25/100), digits = dp), nsmall = dp)
  FIREemiss_som_gCm2yr = format(round(quantile(apply(states_all$FIREemiss_som_gCm2day,1,mean), prob=quantiles_wanted) * (365.25/100), digits = dp), nsmall = dp)
  fire_gCm2yr = format(round(quantile(apply(states_all$fire_gCm2day,1,mean), prob=quantiles_wanted) * (365.25/100), digits = dp), nsmall = dp)
  # Net fluxes
  nbe_gCm2yr = format(round(quantile(apply((states_all$fire_gCm2day+states_all$rhet_litter_gCm2day+states_all$rhet_som_gCm2day+states_all$rauto_gCm2day)-states_all$gpp_gCm2day,1,mean), prob=quantiles_wanted) * (365.25/100), digits = dp), nsmall = dp)
  nee_gCm2yr = format(round(quantile(apply((states_all$rhet_litter_gCm2day+states_all$rhet_som_gCm2day+states_all$rauto_gCm2day)-states_all$gpp_gCm2day,1,mean), prob=quantiles_wanted) * (365.25/100), digits = dp), nsmall = dp)
  # STOCKS
  labile_foliage_gCm2 = format(round(quantile(apply(states_all$labile_gCm2 + states_all$foliage_gCm2,1,mean), prob=quantiles_wanted) /100, digits = dp), nsmall = dp)
  roots_gCm2 = format(round(quantile(apply(states_all$roots_gCm2,1,mean), prob=quantiles_wanted) /100, digits = dp), nsmall = dp)
  wood_gCm2 = format(round(quantile(apply(states_all$wood_gCm2,1,mean), prob=quantiles_wanted) /100, digits = dp), nsmall = dp)
  litter_gCm2 = format(round(quantile(apply(states_all$litter_gCm2,1,mean), prob=quantiles_wanted) /100, digits = dp), nsmall = dp)
  som_gCm2 = format(round(quantile(apply(states_all$som_gCm2,1,mean), prob=quantiles_wanted) /100, digits = dp), nsmall = dp)
  
}
else {
  # We have a compatibility problem
  stop("PROJECT$spatial_type does not have a compatible value, i.e. grid or site")
} # end if is_grid

###
## Begin writting the latex code

col_sep = ""
nos_cols = 20

write(    c("\\documentclass{article}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = FALSE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\title{C budget template figure}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\usepackage[T1]{fontenc}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\usepackage{verbatim}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\usepackage{color}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\usepackage{hyperref}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\usepackage{cleveref}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\usepackage{fixmath}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\usepackage{ulem}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\usepackage{lscape}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\usepackage{subfigure}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\usepackage{array,multirow,graphicx}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\usepackage{chngcntr}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\usepackage[final]{changes}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\usepackage[A4,landscape, margin=2cm]{geometry}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\usepackage{helvet}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\renewcommand{\\familydefault}{\\sfdefault}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\begin{document}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("% Biogenic flux and emissions figure with budgets"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\begin{figure}[]"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("% NOTE: \\put(x coord,y coord){ ... } where to put ..."), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("%       \\vector(x slope,y slope){length} used for arrows"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("%       \\line(x slope,y slope){length} used for lines"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("%       \\framebox(x,y){...} puts ... at box centre"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("   \\centering"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("   \\fbox{"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\setlength{\\unitlength}{0.60cm}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\begin{picture}(33,12)"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % GPP"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(0.25,4.50){\\vector(1,0){3.25}} % arrow for GPP"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(1.0,4.7){$GPP$}               % Label"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(0.25,3.9){\\small ",gpp_gCm2yr[2],"}                 % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(1.90,3.925){\\scriptsize $\\frac{",gpp_gCm2yr[1],"}{",gpp_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Ra"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(26.0,0.78){$Ra$}              % Label"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(25.7,0.1){\\small ",rauto_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(27.4,0.2){\\scriptsize $\\frac{",rauto_gCm2yr[1],"}{",rauto_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % CUE and associated arrows         "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(3.67,4.7){$CUE$} % Label"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(3.47,3.9){\\small ",cue[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(4.465,3.925){\\scriptsize $\\frac{",cue[1],"}{",cue[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(3.50,3.70){\\dashbox{0.2}(1.75,1.65)} % Box around CUE"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(4.375,3.70){\\line(0,-1){3.05}} % Vertical line down from box"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(4.375,0.65){\\vector(1,0){21.625}} % horizontal line to Ra"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)         
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Add labels and box for the internal C-cycle dynamics"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(3.30,7.70){\\textit{Internal carbon rates}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(3.25,0){\\dashbox{0.2}(22.25,8.2)}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Add labels for input / output"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(0.25,11.5){\\textit{Input}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(0.25,11.0){\\textit{carbon}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(0.25,10.5){\\textit{rates}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(28.5,11.5){\\textit{Output}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(28.5,11.0){\\textit{carbon}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(28.5,10.5){\\textit{rates}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Total NPP"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(5.25,4.5){\\vector(1,0){3.0}} % arrow for NPP"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(6.0,4.7){$NPP$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(5.35,3.9){\\small ",npp_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(6.9,3.925){\\scriptsize $\\frac{",npp_gCm2yr[1],"}{",npp_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Add partitioning point in graphic"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(8.70,4.5){\\small \\rotatebox[origin=c]{90}{NPP allocation}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(8.25,1.5){\\dashbox{0.2}(1.2,6)} % Add box around label"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % NPP to labile + foliage - change by model?"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(9.45,7){\\vector(1,0){3.3}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(9.65,7.2){\\small $NPP_{fol+lab}$}    % Label"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(9.88,6.4){\\small ",npp_labilefoliage_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(11.02,6.45){\\scriptsize $\\frac{",npp_labilefoliage_gCm2yr[1],"}{",npp_labilefoliage_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % NPP to fine root"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(9.45,4.5){\\vector(1,0){3.3}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(9.75,4.7){\\small $NPP_{root}$} % Label"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(9.85,3.9){\\small ",npp_roots_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(11.07,3.95){\\scriptsize $\\frac{",npp_roots_gCm2yr[1],"}{",npp_roots_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % NPP to wood"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(9.45,2){\\vector(1,0){3.3}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(9.75,2.2){\\small $NPP_{wood}$} % Label"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(9.85,1.45){\\small ",npp_wood_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(11.07,1.50){\\scriptsize $\\frac{",npp_wood_gCm2yr[1],"}{",npp_wood_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Foliage + labile C pools"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(12.70,6.0){\\framebox(2.5,2.05)}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(12.95,7.45){$C_{fol+lab}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)              
write(paste("         \\put(13.30,6.80){\\small ",labile_foliage_gCm2[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(13.40,6.20){\\scriptsize $\\frac{",labile_foliage_gCm2[1],"}{",labile_foliage_gCm2[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Fine root C pools"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(12.70,3.5){\\framebox(2.5,2.05)}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(13.20,4.95){$C_{root}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)              
write(paste("         \\put(13.30,4.35){\\small ",roots_gCm2[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(13.40,3.75){\\scriptsize $\\frac{",roots_gCm2[1],"}{",roots_gCm2[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Wood C pools"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(12.70,1){\\framebox(2.5,2.05)}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(13.1,2.45){$C_{wood}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)              
write(paste("         \\put(13.10,1.85){\\small ",wood_gCm2[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(13.20,1.25){\\scriptsize $\\frac{",wood_gCm2[1],"}{",wood_gCm2[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(    c("         % litter/mortality, natural and fire driven (red)"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Labile + foliage"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(15.7,7.2){\\small $Mort_{fol+lab}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(15.65,6.4){\\small ",foliage_to_litter_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(16.85,6.45){\\scriptsize $\\frac{",foliage_to_litter_gCm2yr[1],"}{",foliage_to_litter_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(paste("         \\put(18.0,6.4){\\color{red}{\\small ",FIRElitter_labile_foliage_gCm2yr[2],"}}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(paste("         \\put(18.95,6.45){\\color{red}{\\scriptsize $\\frac{",FIRElitter_labile_foliage_gCm2yr[1],"}{",FIRElitter_labile_foliage_gCm2yr[3],"}$}} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Assign arrow to Clitter"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(15.25,7.0){\\line(1,0){4.45}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(22.20,7.0){\\vector(1,-1){0.6}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Fine root"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(15.7,4.65){\\small $Mort_{root}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(15.65,3.9){\\small ",roots_to_litter_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(16.85,3.95){\\scriptsize $\\frac{",roots_to_litter_gCm2yr[1],"}{",roots_to_litter_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(paste("         \\put(18.0,3.9){\\color{red}{\\small ",FIRElitter_roots_gCm2yr[2],"}}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(paste("         \\put(18.95,3.95){\\color{red}{\\scriptsize $\\frac{",FIRElitter_roots_gCm2yr[1],"}{",FIRElitter_roots_gCm2yr[3],"}$}} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Assign arrow to Clitter"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(15.25,4.5){\\line(1,0){4.45}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(22.20,4.5){\\vector(1,2){0.6}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Wood"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(15.70,2.2){\\small $Mort_{wood}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(15.65,1.45){\\small ",wood_to_litter_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(16.85,1.50){\\scriptsize $\\frac{",wood_to_litter_gCm2yr[1],"}{",wood_to_litter_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(paste("         \\put(18.0,1.45){\\color{red}{\\small ",FIRElitter_wood_gCm2yr[2],"}}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(paste("         \\put(18.95,1.50){\\color{red}{\\scriptsize $\\frac{",FIRElitter_wood_gCm2yr[1],"}{",FIRElitter_wood_gCm2yr[3],"}$}} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Assign arrow to Csom"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(15.25,2.0){\\line(1,0){4.45}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(22.20,2.0){\\vector(1,1){0.6}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Mean Transient Times Live Pools"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % MRT foliage"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(20.95,7.0){\\color{black}\\circle{2.5}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(20,7.3){\\small $MRT_{fol}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)              
write(paste("         \\put(20.7,6.7){\\small ",mtt_fol_years[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(20.7,6.2){\\scriptsize $\\frac{",mtt_fol_years[1],"}{",mtt_fol_years[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % MRT root"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(20.95,4.5){\\color{black}\\circle{2.5}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(20,4.7){\\small $MRT_{root}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)              
write(paste("         \\put(20.7,4.1){\\small ",mtt_root_years[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(20.7,3.6){\\scriptsize $\\frac{",mtt_root_years[1],"}{",mtt_root_years[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % MRT wood"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(20.95,2.0){\\color{black}\\circle{2.5}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(19.8,2.2){\\small $MRT_{wood}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)              
write(paste("         \\put(20.5,1.6){\\small ",mtt_wood_years[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(20.5,1.1){\\scriptsize $\\frac{",mtt_wood_years[1],"}{",mtt_wood_years[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Fire emission fluxes"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Labile + foliage"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(    c("         \\put(13.37,9.4){\\small \\rotatebox[origin=c]{90}{$E_{fol+lab}$}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(    c("         \\put(14.0,8.1){\\color{red}{\\line(0,1){2.9}}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(paste("         \\put(14.1,9.7){\\color{red}{\\small ",FIREemiss_labile_foliage_gCm2yr[2],"}}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(paste("         \\put(14.1,9.1){\\color{red}{\\scriptsize $\\frac{",FIREemiss_labile_foliage_gCm2yr[1],"}{",FIREemiss_labile_foliage_gCm2yr[3],"}$}} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Fine roots"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(    c("         \\put(11.85,9.4){\\small \\rotatebox[origin=c]{90}{$E_{root}$}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(    c("         \\put(12.45,5.0){\\color{red}{\\line(1,0){0.25}}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(    c("         \\put(12.45,5.0){\\color{red}{\\line(0,1){6.0}}} "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(paste("         \\put(12.50,9.7){\\color{red}{\\small ",FIREemiss_roots_gCm2yr[2],"}}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(paste("         \\put(12.45,9.1){\\color{red}{\\scriptsize $\\frac{",FIREemiss_roots_gCm2yr[1],"}{",FIREemiss_roots_gCm2yr[3],"}$}} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Wood"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(    c("         \\put(15.05,9.4){\\small \\rotatebox[origin=c]{90}{$E_{wood}$}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(    c("         \\put(15.2,2.7){\\color{red}{\\line(1,0){0.35}}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(    c("         \\put(15.55,2.7){\\color{red}{\\line(0,1){8.3}}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(paste("         \\put(15.65,9.7){\\color{red}{\\small ",FIREemiss_wood_gCm2yr[2],"}}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(paste("         \\put(15.70,9.1){\\color{red}{\\scriptsize $\\frac{",FIREemiss_wood_gCm2yr[1],"}{",FIREemiss_wood_gCm2yr[3],"}$}} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Litter"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(    c("         \\put(22.75,9.4){\\small \\rotatebox[origin=c]{90}{$E_{litter}$}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(    c("         \\put(23.25,7.2){\\color{red}{\\line(0,1){3.8}}} "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(paste("         \\put(23.35,9.7){\\color{red}{\\small ",FIREemiss_litter_gCm2yr[2],"}}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(paste("         \\put(23.35,9.2){\\color{red}{\\scriptsize $\\frac{",FIREemiss_litter_gCm2yr[1],"}{",FIREemiss_litter_gCm2yr[3],"}$}} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % SOM"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(    c("         \\put(25,9.4){\\small \\rotatebox[origin=c]{90}{$E_{som}$}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(    c("         \\put(25.25,3.4){\\color{red}{\\line(1,0){0.3}}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(    c("         \\put(25.55,3.4){\\color{red}{\\line(0,1){7.6}}} "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(paste("         \\put(25.6,9.7){\\color{red}{\\small ",FIREemiss_som_gCm2yr[2],"}}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(paste("         \\put(25.6,9.2){\\color{red}{\\scriptsize $\\frac{",FIREemiss_som_gCm2yr[1],"}{",FIREemiss_som_gCm2yr[3],"}$}} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Fire emissions total"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(    c("         \\put(25.9,11.2){$E_{total}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(    c("         \\put(12.45,11.0){\\color{red}{\\vector(1,0){13.45}}} "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(paste("         \\put(25.9,10.5){\\color{red}{\\small ",fire_gCm2yr[2],"}}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(paste("         \\put(27.10,10.5){\\color{red}{\\scriptsize $\\frac{",fire_gCm2yr[1],"}{",fire_gCm2yr[3],"}$}} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Litter C pool"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(22.75,5.2){\\framebox(2.5,2)}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(23.15,6.55){$C_{litter}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(23.4,5.95){\\small ",litter_gCm2[2],"}               % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(23.4,5.45){\\scriptsize $\\frac{",litter_gCm2[1],"}{",litter_gCm2[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Som C pool"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(22.75,1.0){\\framebox(2.5,2.75)}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(23.15,2.9){$C_{som}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(23,2.1){\\small ",som_gCm2[2],"}                 % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(23.2,1.45){\\scriptsize $\\frac{",som_gCm2[1],"}{",som_gCm2[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Decomposition - natural and fire"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(24.0,5.2){\\vector(0,-1){1.40}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Decomposition - natural"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(24.2,4.65){\\small ",litter_to_som_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(24.15,4.1){\\scriptsize $\\frac{",litter_to_som_gCm2yr[1],"}{",litter_to_som_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Decomposition - combusted litter to som"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(paste("         \\put(22.9,4.65){\\color{red}{\\small ",FIRElitter_litter_gCm2yr[2],"}}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
# write(paste("         \\put(23,4.1){\\color{red}{\\scriptsize $\\frac{",FIRElitter_litter_gCm2yr[1],"}{",FIRElitter_litter_gCm2yr[3],"}$}} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Heterotrophic respiration of litter"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(25.25,6.0){\\line(1,0){2.25}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(25.55,6.3){\\small $Rh_{litter}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(25.55,5.45){\\small ",rhet_litter_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(25.6,4.95){\\scriptsize $\\frac{",rhet_litter_gCm2yr[1],"}{",rhet_litter_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(30.0,6.0){\\vector(1,-2){0.70}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Heterotrophic respiration of som"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(25.25,2.75){\\line(1,0){2.25}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(25.55,3.0){\\small $Rh_{som}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(25.55,2.2){\\small ",rhet_som_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(26.15,1.75){\\scriptsize $\\frac{",rhet_som_gCm2yr[1],"}{",rhet_som_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(30.0,2.75){\\vector(1,1){0.85}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Mean Transient Times Dead Pools"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % MRT litter"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(28.75,6.0){\\color{black}\\circle{2.5}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(27.55,6.0){\\small $MRT_{litter}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)              
write(paste("         \\put(28.4,5.5){\\small ",mtt_litter_years[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(28.4,5.0){\\scriptsize $\\frac{",mtt_litter_years[1],"}{",mtt_litter_years[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % MRT som"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(28.75,2.75){\\color{black}\\circle{2.5}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(27.7,3){\\small $MRT_{som}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)              
write(paste("         \\put(28.3,2.4){\\small ",mtt_som_years[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(28.3,1.9){\\scriptsize $\\frac{",mtt_som_years[1],"}{",mtt_som_years[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Total heterotrophic respiration"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(31.0,4.6){$Rh$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(30.2,4.0){\\small ",rhet_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(31.9,4.0){\\scriptsize $\\frac{",rhet_gCm2yr[1],"}{",rhet_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Net Biome Exchange "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(3.25,8.5){\\dashbox{0.2}(2.2,2.2)}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(3.75,10.05){NBE}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(3.75,9.4){\\small ",nbe_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(3.75,8.8){\\scriptsize $\\frac{",nbe_gCm2yr[1],"}{",nbe_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Net Ecosystem Exchange "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(6.25,8.5){\\dashbox{0.2}(2.2,2.2)}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(6.75,10.05){NEE}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(6.75,9.4){\\small ",nee_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(6.75,8.8){\\scriptsize $\\frac{",nee_gCm2yr[1],"}{",nee_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\end{picture}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("        }"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(paste("   \\caption{",figure_caption,"}", sep="")), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(paste("   \\label{",figure_label,"}"), sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\end{figure}   "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\end{document}   "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
print(paste(e," ",chartr("_"," ",PROJECT$name)," DONE",sep=""))
}

setwd("/home/cnwobi/CARDAMOM/cssp_rainfor_brazil")
# DONE!
