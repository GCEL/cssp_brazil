######Read in shape files ###########
amazon_nw_poly <- shapefile("./amazon_nw.shp");amazon_sw_poly <- shapefile("./amazon_sw.shp");amazon_ec_poly <- shapefile("./amazon_ec.shp")
amazon_bs_poly <- shapefile("./amazon_bs.shp");amazon_gs_poly <- shapefile("./amazon_gs.shp")
amazonia_poly <-rasterToPolygons(biommort_10_16)
# amazon_nw_raster<-crop(biommort_10_16,amazon_nw_poly)
# amazon_sw_raster<-crop(biommort_10_16,amazon_sw_poly)
# amazon_ec_raster<-crop(biommort_10_16,amazon_ec_poly)
# amazon_bs_raster<-crop(biommort_10_16,amazon_bs_poly)
# amazon_gs_raster<-crop(biommort_10_16,amazon_gs_poly)

# function to extract country amazon regions
country_to_amazon_crop_fun <- function(a,s) { 
  cropped <- crop(a, extent(s))
  masked1 <- mask(cropped, s)
  masked2 <- masked1 > -Inf
  poly <- rasterToPolygons(masked2, dissolve=TRUE)
}

amazon_nw_pol <- country_to_amazon_crop_fun(biommort_10_16,amazon_nw_poly)
amazon_sw_pol <- country_to_amazon_crop_fun(biommort_10_16,amazon_sw_poly)
amazon_ec_pol <- country_to_amazon_crop_fun(biommort_10_16,amazon_ec_poly)
amazon_bs_pol <- country_to_amazon_crop_fun(biommort_10_16,amazon_bs_poly)
amazon_gs_pol <- country_to_amazon_crop_fun(biommort_10_16,amazon_gs_poly)

plot(amazon_nw_pol,add=T, col='pink');plot(amazon_sw_pol,add=T, col='white');plot(amazon_ec_pol,add=T, col='green')
plot(amazon_bs_pol,add=T, col='red');plot(amazon_gs_pol,add=T, col='orange')

#create mask
mask_remap_fun <- function(p) {
amazon_extent <- extent(c(-79.7704,-42.87695,-20.93973,10.96843)) #use south america extent from INLAND SA
r <- raster(amazon_extent)
res(r)<- 1 #use resolution to mirror other masks in ilamb/ could be 1x1
values(r) <- 1
crs(r) <- "+proj=longlat +datum=WGS84 +no_defs"
r_amazon <- mask(r, p, updatevalue=0)#create mask from polygon of RAINFOR amazon region
}

r_amazon_nw_pol <- mask_remap_fun(amazon_nw_pol)
writeRaster(r_amazon_nw_pol, "mask-NORTHWEST-remap.nc", overwrite=TRUE, format="CDF",varname='ids',datatype='LOG1S', xname="lon", yname="lat")
r_amazon_sw_pol <- mask_remap_fun(amazon_sw_pol)
writeRaster(r_amazon_sw_pol, "mask-SOUTHWEST-remap.nc", overwrite=TRUE, format="CDF",varname='ids',datatype='LOG1S', xname="lon", yname="lat")
r_amazon_ec_pol <- mask_remap_fun(amazon_ec_pol)
writeRaster(r_amazon_ec_pol, "mask-EASTCENTRAL-remap.nc", overwrite=TRUE, format="CDF",varname='ids',datatype='LOG1S', xname="lon", yname="lat")
r_amazon_bs_pol <- mask_remap_fun(amazon_bs_pol)
writeRaster(r_amazon_bs_pol, "mask-BRAZILSHIELD-remap.nc", overwrite=TRUE, format="CDF",varname='ids',datatype='LOG1S', xname="lon", yname="lat")
r_amazon_gs_pol <- mask_remap_fun(amazon_gs_pol)
writeRaster(r_amazon_gs_pol, "mask-GUIANASHIELD-remap.nc", overwrite=TRUE, format="CDF",varname='ids',datatype='LOG1S', xname="lon", yname="lat")
r_amazonia_pol <- mask_remap_fun(amazonia_poly)
writeRaster(r_amazonia_pol, "mask-AMAZONIA-remap.nc", overwrite=TRUE, format="CDF",varname='ids',datatype='LOG1S', xname="lon", yname="lat")

nc_masks <- nc_open('./mask-GUIANASHIELD-remap.nc') #gcel
nc_masks
nc_close(nc_masks)
