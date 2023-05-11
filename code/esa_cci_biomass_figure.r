
esa_cci_data = brick("G://AGB/ESA_CCI_BIOMASS/ESA_CCI_AGB_1deg/AGB_map_MgCha_2010.tif")
target = brick("G://cssp_rainfor_amazon_brazil/rainfor_leeds_data/AbovegroundBiomass_Mg_perHa_111km.tif")
esa_cci_data_resamp = resample(esa_cci_data, target, method = "bilinear", na.rm=TRUE)
plot(target);plot(esa_cci_data_resamp)
plot(esa_cci_data)

esa_cci_data_resamp_gCm2 <- calc(esa_cci_data_resamp, tha_to_gCm2_fun)

esa_cci_biomass_C_2010_amazonia_extent <- extract_subset(amazonia_extent,esa_cci_data_resamp)
# esa_cci_biomass_C_2010_amazonia_extent <- esa_cci_biomass_2010_amazonia_extent*0.48
plot(esa_cci_biomass_C_2010_amazonia_extent)

esa_cci_biomass_C_2010_amazonia_ifl <- extract_subset(amazonia_ifl_layer,esa_cci_biomass_2010_amazonia_extent)
# esa_cci_biomass_C_2010_amazonia_ifl <- esa_cci_biomass_2010_amazonia_ifl*0.48
plot(esa_cci_biomass_C_2010_amazonia_ifl)

amazon_extent_nw_poly <- shapefile("C://Users/nwobi/OneDrive - University of Edinburgh/brazil_leeds_maps/amazonia_extent_nw.shp")
amazon_extent_sw_poly <- shapefile("C://Users/nwobi/OneDrive - University of Edinburgh/brazil_leeds_maps/amazonia_extent_sw.shp")
amazon_extent_ec_poly <- shapefile("C://Users/nwobi/OneDrive - University of Edinburgh/brazil_leeds_maps/amazonia_extent_ec.shp")
amazon_extent_bs_poly <- shapefile("C://Users/nwobi/OneDrive - University of Edinburgh/brazil_leeds_maps/amazonia_extent_bs.shp")
amazon_extent_gs_poly <- shapefile("C://Users/nwobi/OneDrive - University of Edinburgh/brazil_leeds_maps/amazonia_extent_gs.shp")
sa_extent_poly <- shapefile("C://Users/nwobi/OneDrive - University of Edinburgh/brazil_leeds_maps/sa_extent_dis.shp")

plot(biomass_amazon_C_ifl);plot(sa_extent_poly,add=T)
plot(biomass_C_amazon);plot(sa_extent_poly,add=T)
plot(amazon_extent_nw_poly,add=T);plot(amazon_extent_sw_poly,add=T);plot(amazon_extent_ec_poly,add=T);plot(amazon_extent_bs_poly,add=T);plot(amazon_extent_gs_poly,add=T)
plot(allpahuayo_poly,add=T);plot(caxiuana_poly,add=T);plot(kenia_poly,add=T);plot(tambopata_poly,add=T);plot(tanguro_poly,add=T)
