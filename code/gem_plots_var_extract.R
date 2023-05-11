##########Masks for Creating Ecoregions##################
amazon_region_gem <- extent(biomass_C_amazon)
r_amazon_region_gem <- raster(amazon_region_gem,res=1)
values(r_amazon_region_gem) <- 0
r_amazon_region_gem <- rasterize(biomass_amazon_mask_pol,r_amazon_region_gem)
plot(r_amazon_region_gem)

allpahuayo_poly <- shapefile("./data/Allpahuayo.shp")
caxiuana_poly <- shapefile("./data/Caxiuana.shp")
kenia_poly <- shapefile("./data/Kenia.shp")
tambopata_poly <- shapefile("./data/Tambopata.shp")
tanguro_poly <- shapefile("./data/Tanguro.shp")

plot(biomass_C_amazon);plot(amazon_nw_poly,add=T);plot(amazon_sw_poly,add=T);plot(amazon_ec_poly,add=T);plot(amazon_bs_poly,add=T);plot(amazon_gs_poly,add=T)
plot(allpahuayo_poly,add=T)
plot(caxiuana_poly,add=T)
plot(kenia_poly,add=T)
plot(tambopata_poly,add=T)
plot(tanguro_poly,add=T)

allpahuayo_pol_mask <- country_to_amazon_crop_fun(biomass_C_amazon,allpahuayo_poly)
caxiuana_pol_mask <- country_to_amazon_crop_fun(biomass_C_amazon,caxiuana_poly)
kenia_pol_mask <- country_to_amazon_crop_fun(biomass_C_amazon,kenia_poly)
tambopata_pol_mask <- country_to_amazon_crop_fun(biomass_C_amazon,tambopata_poly)
tanguro_pol_mask <- country_to_amazon_crop_fun(biomass_C_amazon,tanguro_poly)

r_allpahuayo <- mask(r_amazon_region_gem, allpahuayo_pol_mask, updatevalue=0)
r_caxiuana <- mask(r_allpahuayo, caxiuana_pol_mask, inverse=T, updatevalue=2)
r_kenia <- mask(r_caxiuana, tambopata_pol_mask, inverse=T, updatevalue=3)
r_tambopata <- mask(r_kenia, kenia_pol_mask, inverse=T, updatevalue=4)
r_amazonia_gem <- mask(r_tambopata, tanguro_pol_mask, inverse=T, updatevalue=5)
# plot(r_amazonia_gem)
# plot(flip(r_amazonia_gem,2))
# plot(t(r_amazonia_gem))
# 
# plot(allpahuayo_poly,add=T);plot(caxiuana_poly,add=T);plot(kenia_poly,add=T);plot(tambopata_poly,add=T);plot(tanguro_poly,add=T)
# plot(amazonia_extent,add=T)

amazonia_gem_matrix <- t(as.matrix(flip(r_amazonia_gem,2)))

amazonia_gem_matrices <- array(numeric(),c(37,32,3)) 
amazonia_gem_matrices[,,1]<-amazonia_gem_matrix
amazonia_gem_matrices[,,2]<-amazonia_gem_matrix
amazonia_gem_matrices[,,3]<-amazonia_gem_matrix

# test_gem <- amazonia_gem_matrices==1
# dim(test_gem)
gem_plot_names <- c('allpahuayo','caxiuana','kenia','tambopata','tanguro')
gem_var_names <- c('gpp_MgChayr','cue','npp_MgChayr','npp_wood_MgChayr','wood_MgCha','mrt_wood_years')
gem_prefix <- "M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/amazonia_"
gem_suffix_1 <- "_nomngt/infofile.RData"
gem_midfix <- "_nomngt/RESULTS_PROCESSED/amazonia_"
gem_suffix_2 <- "_nomngt_stock_flux.RData"

for (i in mod_var) {
  gem_plots_df <- data.frame(matrix(nrow=6,ncol=length(gem_plot_names)))
  colnames(gem_plots_df) <- gem_plot_names
  rownames(gem_plots_df) <- gem_var_names

for (e in gem_plot_names) {
  # Load the processed site file
  load(paste(gem_prefix,i,gem_suffix_1,sep=""))
  load(paste(gem_prefix,i,gem_midfix,i,gem_suffix_2,sep=""))
  
  if (e =='allpahuayo') {
    gem_plot_mask <- amazonia_gem_matrices==1
  }
  else if (e =='caxiuana') {
    gem_plot_mask <- amazonia_gem_matrices==2
  }
  else if (e =='kenia') {
    gem_plot_mask <- amazonia_gem_matrices==3
  }
  else if (e =='tambopata') {
    gem_plot_mask <- amazonia_gem_matrices==4
  }
  else {
    gem_plot_mask <- amazonia_gem_matrices==5
  }
dp = 2
quantiles_wanted = c(1,4,7)

#Masking phase
# NATURAL FLUXES
grid_output$mean_gpp_gCm2day[,,quantiles_wanted][!gem_plot_mask]<-NA
grid_output$mean_cue[,,quantiles_wanted][!gem_plot_mask]<-NA
grid_output$mean_npp_gCm2day[,,quantiles_wanted][!gem_plot_mask]<-NA
grid_output$mean_alloc_wood_gCm2day[,,quantiles_wanted][!gem_plot_mask]<-NA
# STOCKS
grid_output$mean_wood_gCm2[,,quantiles_wanted][!gem_plot_mask]<-NA
# MEAN TRANSIENT TIMES
grid_output$MTT_wood_years[,,quantiles_wanted][!gem_plot_mask]<-NA

# Extract or calculate required derived values
# NOTE: unit conversion from gC/m2/day to MgC/ha/yr
# NATURAL FLUXES
gpp_MgChayr = format(round(apply(grid_output$mean_gpp_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
cue = format(round(apply(grid_output$mean_cue[,,quantiles_wanted],3,mean, na.rm=TRUE), digits = 2), nsmall = dp)
npp_MgChayr = format(round(apply(grid_output$mean_npp_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
npp_wood_MgChayr = format(round(apply(grid_output$mean_alloc_wood_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
# STOCKS
wood_MgCha = format(round(apply(grid_output$mean_wood_gCm2[,,quantiles_wanted],3,mean, na.rm=TRUE) / 100, digits = dp), nsmall = dp)
# MEAN TRANSIENT TIMES
mtt_wood_years = format(round(apply(grid_output$MTT_wood_years[,,quantiles_wanted],3,mean, na.rm=TRUE), digits = dp), nsmall = dp)

#collate all variables
all_var_means <- as.numeric(c(gpp_MgChayr[2],cue[2],npp_MgChayr[2],npp_wood_MgChayr[2],wood_MgCha[2],mtt_wood_years[2]))
#add to dataframe
gem_plots_df[e]<-all_var_means
}
  assign(paste("gem_plots_df_",i,sep=""),gem_plots_df)
}
write.csv(gem_plots_df_esa_cci_agb,"gem_plots_df_esa_cci_agb.csv")
write.csv(gem_plots_df_rainfor_biomass_productivity_2005,"gem_plots_df_rainfor_biomass_productivity_2005.csv")
write.csv(gem_plots_df_rainfor_biomass_annual_productivity,"gem_plots_df_rainfor_biomass_annual_productivity.csv")



