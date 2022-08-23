
library(ncdf4) ; library(raster)
biomass_amazon <- brick('G://amazon_rainfor_leeds_maps/AbovegroundBiomass_Mg_perHa_111km.tif')
biomass_amazon #leeds amazon biomass use extent for amazon extent
plot(biomass_amazon)

tha_to_gm2_fun <- function(x) { x * 100 }
biomass_amazon_gm2 <- calc(biomass_amazon, tha_to_gm2_fun)
plot(biomass_amazon_gm2)

#convert to a polygon
biomass_amazon_mask <- biomass_amazon > -Inf
biomass_amazon_mask_pol <- rasterToPolygons(biomass_amazon_mask, dissolve=TRUE)
plot(biomass_amazon_mask_pol)

amazon_mask <- brick('G://ILAMB_beta_devel/ILAMB_beta_tutorial/ilamb-mask-AMAZON.nc')
amazon_mask #aim to mirror this brazil amazon mask in properties
summary(amazon_mask)

#create mask
SA_extent <- extent(c(-84,-32,-57,17)) #use south america extent from INLAND SA
SA_extent_r <- raster(SA_extent)
res(SA_extent_r)<- 0.5 #use resolution to mirror other masks in ilamb/ could be 1x1
values(SA_extent_r) <- -9.223372e+18 #set value to non-mask area
crs(SA_extent_r)<-"+proj=longlat +datum=WGS84 +no_defs" 
plot(SA_extent_r)

amazon_region_mask <- mask(SA_extent_r, biomass_amazon_mask_pol, updatevalue=0, inverse=T)
amazon_region_mask
summary(amazon_region_mask)
plot(amazon_region_mask)
summary(amazon_region_mask)

writeRaster(amazonia_mask, "amazonia_mask.nc", overwrite=TRUE, format="CDF",varname='ids',xname="lon",yname="lat")
nc_amazonia_mask <- nc_open('./amazonia_mask.nc',write=TRUE)
nc_close(nc_amazonia_mask)

#to nc
writeRaster(amazonia_mask_ilamb, "ilamb-mask-AMAZONIA.nc", overwrite=TRUE, format="CDF",varname='ids', xname="lon", yname="lat")

nc_leeds_amazonmask <- nc_open('./ilamb-mask-AMAZONIA.nc',write=TRUE)
print(nc_leeds_amazonmask) #print properties
attributes(nc_leeds_amazonmask$dim)$names

# xdim2 <- nc_leeds_amazonmask$dim[['lon']]
# ydim2 <- nc_leeds_amazonmask$dim[['lat']]
# 
# mv2 <- "labels"
# var_q <- ncvar_def( 'ids', '', list(xdim2,ydim2), mv2 , prec="integer")
# var_cf <- ncvar_def( 'labels', '', list() )

cnames   <-c("amazonia")
nstrings <- length(cnames)
dimnchar   <- ncdim_def("nb",   "", 1:12, create_dimvar=FALSE )
dimregion <- ncdim_def("n", "", 1:nstrings, create_dimvar=FALSE )
#varcolors  <- ncvar_def("colours", "", list(dimnchar,dimregion) )
varregion  <- ncvar_def("labels", "", list(dimnchar,dimregion) ,prec="char")

#varcolors  <- ncvar_def("colors", "", list(dimnchar, dimcolorno),prec="char" )

#c_leeds_amazonmask <- ncvar_add( nc_leeds_amazonmask, varcolors)	# NOTE this returns a modified netcdf file handle

nc_leeds_amazonmask <- ncvar_add(nc_leeds_amazonmask, varamazon)

ncvar_put(varregion, "labels", cnames, verbose=TRUE)

ncid <- nc_create( "regionnames.nc", list(varregion))
nc_amazonia_names <- nc_open('./regionnames.nc',write=TRUE)
nc_close(nc_amazonia_names)

ncvar_put( ncid, "labels", cnames, verbose=TRUE )

attributes(nc_leeds_amazonmask)$names
attributes(nc_leeds_amazonmask$var)$names
ncatt_get(nc_leeds_amazonmask, attributes(nc_leeds_amazonmask$var)$names[3])

nc_close(nc_leeds_amazonmask) #close nc file

amazonia_mask <- brick('./amazon_leeds_mask.nc')
amazonia_mask
summary(amazonia_mask)

writeRaster(biomass_amazon_gm2, "leeds_biomass_amazon.nc", overwrite=TRUE, format="CDF",varname="WOOD", varunit="g.m-2", longname="Wood C", xname="lon", yname="lat")
nc_leeds_agb_amazon <- nc_open('./leeds_biomass_amazon.nc')
print(nc_leeds_agb_amazon) #print properties
nc_close(nc_leeds_agb_amazon) #close nc file

#default brazilian amazon mask in ilamb devel
nc_amazonmask <- nc_open('G:/ILAMB_runs_output/CSSP_stippling/CARDAMOM_monthly_1x1_SAmerica.nc') #gcel
nc_amazonmask
# summary(nc_amazonmask)
# print(nc_amazonmask) #print properties
# attributes(nc_amazonmask)$names
# attributes(nc_amazonmask$var)$names
# ncatt_get(nc_amazonmask, attributes(nc_amazonmask$var)$names[1])
nc_close(nc_amazonmask) #close nc file


amazon_mask <- raster('./mask-AMAZONIA-remap.nc')
amazon_mask
plot(amazon_mask)
summary(amazon_mask)
