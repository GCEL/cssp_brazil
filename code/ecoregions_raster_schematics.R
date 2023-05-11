##########Masks for Creating Ecoregions##################
amazon_region_new <- extent(biomass_amazon_C_ifl)
r_amazon_region_new <- raster(amazon_region_new,res=1)
values(r_amazon_region_new) <- 0
r_amazon_region_final <- rasterize(biomass_amazon_mask_pol,r_amazon_region_new)
plot(r_amazon_region_final)

amazon_nw_pol_mask <- country_to_amazon_crop_fun(biomass_amazon_C_ifl,amazon_nw_poly)
amazon_sw_pol_mask <- country_to_amazon_crop_fun(biomass_amazon_C_ifl,amazon_sw_poly)
amazon_ec_pol_mask <- country_to_amazon_crop_fun(biomass_amazon_C_ifl,amazon_ec_poly)
amazon_bs_pol_mask <- country_to_amazon_crop_fun(biomass_amazon_C_ifl,amazon_bs_poly)
amazon_gs_pol_mask <- country_to_amazon_crop_fun(biomass_amazon_C_ifl,amazon_gs_poly)

r_amazonia_nw <- mask(r_amazon_region_final, amazon_nw_pol_mask, updatevalue=0)
r_amazonia_sw <- mask(r_amazonia_nw, amazon_sw_pol_mask, inverse=T, updatevalue=2)
r_amazonia_bs <- mask(r_amazonia_sw, amazon_bs_pol_mask, inverse=T, updatevalue=3)
r_amazonia_ec <- mask(r_amazonia_bs, amazon_ec_pol_mask, inverse=T, updatevalue=4)
r_amazonia_ecoregions <- mask(r_amazonia_ec, amazon_gs_pol_mask, inverse=T, updatevalue=5)
plot(r_amazonia_ecoregions)
plot(flip(r_amazonia_ecoregions,2))
plot(t(r_amazonia_ecoregions))

plot(amazon_nw_poly,add=T);plot(amazon_sw_poly,add=T);plot(amazon_ec_poly,add=T);plot(amazon_bs_poly,add=T);plot(amazon_gs_poly,add=T)
plot(amazonia_extent,add=T)

###convert to array###
# r_amazonia_ecoregions
# dim(matrix(r_amazonia_ecoregions))
# dim(array(r_amazonia_ecoregions))
# dim(r_amazonia_ecoregions)
# dim(grid_output$landmask)
# array(r_amazonia_ecoregions,dim=c(37,32))

amazonia_ecoregions_matrix <- t(as.matrix(flip(r_amazonia_ecoregions,2)))
test_single <- amazonia_ecoregions_matrix==1
  
amazonia_ecoregions_matrices <- array(numeric(),c(37,32,3)) 
amazonia_ecoregions_matrices[,,1]<-amazonia_ecoregions_matrix
amazonia_ecoregions_matrices[,,2]<-amazonia_ecoregions_matrix
amazonia_ecoregions_matrices[,,3]<-amazonia_ecoregions_matrix

test_region <- amazonia_ecoregions_matrices==1
dim(test_region)

load("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/amazonia_rainfor_biomass_annual_productivity_nomngt//infofile.RData")
load(paste("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/amazonia_rainfor_biomass_annual_productivity_nomngt/RESULTS_PROCESSED/","amazonia_rainfor_biomass_annual_productivity_nomngt","_stock_flux.RData",sep=""))

grid_output$landmask
dim(grid_output$landmask)

nw_landmask <- grid_output$landmask
nw_landmask[!test_single]<- NA

new_extent <- extent(biomass_ifl_subset)

r_landmask <- raster(t(grid_output$landmask))
r_landmask_nw <- raster(t(nw_landmask))

extent(r_landmask) <- new_extent
extent(r_landmask_nw) <- new_extent
# c(1,4,7)

plot(r_landmask)
plot(flip(r_landmask,2))
# plot(t(r_landmask))
plot(r_landmask_nw)
plot(flip(r_landmask_nw,2))

dim(grid_output$mean_gpp_gCm2day[,,c(1,4,7)])
grid_output$mean_gpp_gCm2day[,,c(1,4,7)][!test_region]<- NA

apply(grid_output$mean_gpp_gCm2day[,,c(1,4,7)][!test_region]<- NA,3,mean, na.rm=TRUE)
