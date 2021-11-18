#############Prepare Mask Template for ILAMB#####################################

#create mask
SA_extent <- extent(c(-84,-32,-57,14)) #use south america extent from INLAND SA
r <- raster(SA_extent)
res(r)<- 0.5 #use resolution to mirror other masks in ilamb/ could be 1x1
values(r) <- 1
crs(r)<-crs(amazon_mask)
r
plot(r)

biomass_amazon <- brick('G://amazon_rainfor_leeds_maps/AbovegroundBiomass_Mg_perHa_111km.tif')
biomass_amazon_mask <- biomass_amazon > -Inf
biomass_amazon_mask_pol <- rasterToPolygons(biomass_amazon_mask, dissolve=TRUE)
plot(biomass_amazon_mask_pol)

r_amazon <- mask(r, biomass_amazon_mask_pol, updatevalue=0) #create mask from polygon of RAINFOR amazon region
r_amazon
summary(r_amazon)
plot(r_amazon)

writeRaster(r_amazon, "mask-AMAZONIA-remap.nc", overwrite=TRUE, format="CDF",varname='ids',datatype='LOG1S', xname="lon", yname="lat")

cardamom_agb <- stack('G://ILAMB_runs_output/CSSP_stippling/DATA/benchmark/CARDAMOM_monthly_1x1_SAmerica.nc',varname="npptot")
class(getZ(cardamom_agb))
plot(cardamom_agb$X2001.01.15)

cardamom_agb_0101<-cardamom_agb$X2001.01.15
plot(cardamom_agb_0101)
cardamom_agb_0101_mask <- cardamom_agb_0101 > -Inf
cardamom_agb_0101_mask_pol <- rasterToPolygons(cardamom_agb_0101_mask, dissolve=TRUE)
plot(cardamom_agb_0101_mask_pol)

r_southamerica <- mask(r, cardamom_agb_0101_mask_pol, updatevalue=0)
r_southamerica
summary(r_southamerica)
plot(r_southamerica)

#writeRaster(r_amazon, "ilamb-mask-AMAZONIA.tif", overwrite=TRUE, format="GTiff")
writeRaster(r_southamerica, "mask-SOUTHAMERICA-remap.nc", overwrite=TRUE, format="CDF",varname='ids',datatype='LOG1S', xname="lon", yname="lat")
