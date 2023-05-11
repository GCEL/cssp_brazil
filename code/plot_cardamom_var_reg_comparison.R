

reg_compare_var_names <- c('GPP','NPP','CUE','NPP_wood_flx','WOOD','OUTPUT_wood_flx','MTT_wood')
amazonia_ifl_layer <- shapefile("R:/brazil_leeds_maps/ifl_2000_amazonia.shp")
amazonia_extent <- shapefile("./data/amazonia_extent")
prefix <- 'R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/amazonia_ifl_'
suffix <- '_nomngt_2000_updated_2010.nc'

reg_compare_mod_var <- c('rainfor_biomass_annual_productivity')

extract_var_median <- function (cardamom_var,model_variant){
  for (i in cardamom_var) {
    data <- brick(paste(prefix,model_variant,suffix,sep=""),varname=cardamom_var)
    data_med <- calc(data, median)
    if (cardamom_var == 'MTT_wood' | cardamom_var == 'CUE') {
      data_med <- data_med
    }
    else if (cardamom_var == 'WOOD'){
      data_med <- data_med/100
    }
    else {
      data_med <- data_med * (365.25/100)
    }
  }
  plot(data_med);plot(sa_extent_poly,add=T);plot(amazon_extent_nw_poly,add=T);plot(amazon_extent_sw_poly,add=T);plot(amazon_extent_ec_poly,add=T);plot(amazon_extent_bs_poly,add=T);plot(amazon_extent_gs_poly,add=T)
}
extract_var_median(reg_compare_var_names[1],reg_compare_mod_var)
extract_var_median(reg_compare_var_names[2],reg_compare_mod_var)
extract_var_median(reg_compare_var_names[3],reg_compare_mod_var)
extract_var_median(reg_compare_var_names[4],reg_compare_mod_var)
extract_var_median(reg_compare_var_names[5],reg_compare_mod_var)
extract_var_median(reg_compare_var_names[6],reg_compare_mod_var)
extract_var_median(reg_compare_var_names[7],reg_compare_mod_var)

load("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/amazonia_ifl_rainfor_biomass_annual_productivity_nomngt//infofile.RData")
load(paste("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/amazonia_ifl_rainfor_biomass_annual_productivity_nomngt/RESULTS_PROCESSED/","amazonia_ifl_rainfor_biomass_annual_productivity_nomngt","_stock_flux.RData",sep=""))

# NATURAL FLUXES
gpp_MgChayr = format(round(apply(grid_output$mean_gpp_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
cue = format(round(apply(grid_output$mean_cue[,,quantiles_wanted],3,mean, na.rm=TRUE), digits = 2), nsmall = dp)
npp_MgChayr = format(round(apply(grid_output$mean_npp_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
npp_wood_MgChayr = format(round(apply(grid_output$mean_alloc_wood_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
wood_to_litter_gCm2yr = format(round(apply(grid_output$mean_wood_to_litter_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25, digits = dp), nsmall = dp)
# STOCKS
wood_MgCha = format(round(apply(grid_output$mean_wood_gCm2[,,quantiles_wanted],3,mean, na.rm=TRUE) / 100, digits = dp), nsmall = dp)
# MEAN TRANSIENT TIMES
mrt_wood_years = format(round(apply(grid_output$MTT_wood_years[,,quantiles_wanted],3,mean, na.rm=TRUE), digits = dp), nsmall = dp)

gpp_schem <- flip(raster(t(grid_output$mean_gpp_gCm2day[,,4]* (365.25/100))),2)
extent(gpp_schem) <- extent(biomass_amazon_C_ifl)
crs(gpp_schem) <- crs(biomass_amazon_C_ifl)
plot(gpp_schem);plot(sa_extent_poly,add=T);plot(amazon_extent_nw_poly,add=T);plot(amazon_extent_sw_poly,add=T);plot(amazon_extent_ec_poly,add=T);plot(amazon_extent_bs_poly,add=T);plot(amazon_extent_gs_poly,add=T)

extract_var_grid_mean <- function (cardamom_var,model_ref) {
  for (v in cardamom_var) {
    if (cardamom_var== 'GPP') {
      r_schem <- flip(raster(t(grid_output$mean_gpp_gCm2day[,,4]* (365.25/100))),2)
    }
    else if (cardamom_var== 'NPP') {
      r_schem <- flip(raster(t(grid_output$mean_npp_gCm2day[,,4]* (365.25/100))),2)
    }
    else if (cardamom_var== 'CUE') {
      r_schem <- flip(raster(t(grid_output$mean_cue[,,4])),2)
    }
    else if (cardamom_var== 'NPP_wood_flx') {
      r_schem <- flip(raster(t(grid_output$mean_alloc_wood_gCm2day[,,4]* (365.25/100))),2)
    }
    else if (cardamom_var== 'WOOD') {
      r_schem <- flip(raster(t(grid_output$mean_wood_gCm2[,,4]/100)),2)
    }
    else if (cardamom_var== 'OUTPUT_wood_flx') {
      r_schem <- flip(raster(t(grid_output$mean_wood_to_litter_gCm2day[,,4]* (365.25/100))),2)
    }
    else {
      r_schem <- flip(raster(t(grid_output$MTT_wood_years[,,4])),2)
    }
    extent(r_schem) <- extent(model_ref)
    crs(r_schem) <- crs(model_ref)
    plot(r_schem);plot(sa_extent_poly,add=T);plot(amazon_extent_nw_poly,add=T);plot(amazon_extent_sw_poly,add=T);plot(amazon_extent_ec_poly,add=T);plot(amazon_extent_bs_poly,add=T);plot(amazon_extent_gs_poly,add=T)
  }
}
extract_var_grid_mean(reg_compare_var_names[1],biomass_amazon_C_ifl)
extract_var_grid_mean(reg_compare_var_names[2],biomass_amazon_C_ifl)
extract_var_grid_mean(reg_compare_var_names[3],biomass_amazon_C_ifl)
extract_var_grid_mean(reg_compare_var_names[4],biomass_amazon_C_ifl)
extract_var_grid_mean(reg_compare_var_names[5],biomass_amazon_C_ifl)
extract_var_grid_mean(reg_compare_var_names[6],biomass_amazon_C_ifl)
extract_var_grid_mean(reg_compare_var_names[7],biomass_amazon_C_ifl)
