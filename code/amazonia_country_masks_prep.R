######Read in shape files ###########
bolivia_poly <- shapefile("./BOL_adm0.shp");brazil_poly <- shapefile("./BRA_adm0.shp");colombia_poly <- shapefile("./COL_adm0.shp")
ecuador_poly <- shapefile("./ECU_adm0.shp");freguy_poly <- shapefile("./GUF_adm0.shp");guyana_poly <- shapefile("./GUY_adm0.shp")
peru_poly <- shapefile("./PER_adm0.shp");suriname_poly <- shapefile("./SUR_adm0.shp");venezuela_poly <- shapefile("./VEN_adm0.shp")
amazon_north_poly <- shapefile("./amazon_north_dis.shp");amazon_south_poly <- shapefile("./amazon_south_dis.shp");amazon_east_poly <- shapefile("./amazon_east_dis.shp")
amazon_west_poly <- shapefile("./amazon_west_dis.shp");amazon_central_poly <- shapefile("./amazon_central_dis.shp")

plot(sa_data$GDAL.Band.Number.1)
plot(biomass_amazon_mask_pol,add=T, col='red')
plot(bolivia_poly,add=T);plot(brazil_poly,add=T);plot(colombia_poly,add=T)
plot(ecuador_poly,add=T);plot(freguy_poly,add=T);plot(guyana_poly,add=T)
plot(peru_poly,add=T);plot(suriname_poly,add=T);plot(venezuela_poly,add=T)
plot(amazon_north_poly,add=T);plot(amazon_south_poly,add=T);plot(amazon_east_poly,add=T);plot(amazon_west_poly,add=T);plot(amazon_central_poly,add=T)
biomass_amazon <- disaggregate(biomass_amazon, fact=2)

# function to extract country amazon regions
country_to_amazon_crop_fun <- function(a,s) { 
  cropped <- crop(a, extent(s))
  masked1 <- mask(cropped, s)
  masked2 <- masked1 > -Inf
  poly <- rasterToPolygons(masked2, dissolve=TRUE)
}
bolivia_amazon_pol <- country_to_amazon_crop_fun(biomass_amazon,bolivia_poly)
brazil_amazon_pol <- country_to_amazon_crop_fun(biomass_amazon,brazil_poly)
colombia_amazon_pol <- country_to_amazon_crop_fun(biomass_amazon,colombia_poly)
ecuador_amazon_pol <- country_to_amazon_crop_fun(biomass_amazon,ecuador_poly)
freguy_amazon_pol <- country_to_amazon_crop_fun(biomass_amazon,freguy_poly)
guyana_amazon_pol <- country_to_amazon_crop_fun(biomass_amazon,guyana_poly)
peru_amazon_pol <- country_to_amazon_crop_fun(biomass_amazon,peru_poly)
suriname_amazon_pol <- country_to_amazon_crop_fun(biomass_amazon,suriname_poly)
venezuela_amazon_pol <- country_to_amazon_crop_fun(biomass_amazon,venezuela_poly)
amazon_north <- country_to_amazon_crop_fun(biomass_amazon,amazon_north_poly)
amazon_south <- country_to_amazon_crop_fun(biomass_amazon,amazon_south_poly)
amazon_east <- country_to_amazon_crop_fun(biomass_amazon,amazon_east_poly)
amazon_west <- country_to_amazon_crop_fun(biomass_amazon,amazon_west_poly)
amazon_central <- country_to_amazon_crop_fun(biomass_amazon,amazon_central_poly)

bolivia_amazon_pol2 <- country_to_amazon_crop_fun(r_southamerica,bolivia_poly)

trial1<-crop(biomass_amazon,bolivia_poly)
plot(trial1)
trial2<-mask(trial1,bolivia_poly)
plot(trial2)
trial2 <- disaggregate(trial2, fact=2)
res(trial2)
trial2
extent(trial2) <- c(-84.0,-32.0,-57.0,13.5)
extent(trial2)
writeRaster(trial2, "trial3.nc", overwrite=TRUE, format="CDF",varname='ids', xname="lon", yname="lat")
nc_trial <- nc_open('./mask-BOLIVIAMAZON-remap.nc') #gcel
nc_trial
ids_array <- ncvar_get(nc_trial,'ids')
# dlname <- ncatt_get(ncin,dname,"long_name")
# dunits <- ncatt_get(ncin,dname,"units")
ids_fillvalue <- ncatt_get(nc_trial,'ids',"_FillValue")
dim(ids_array)
# replace netCDF fill values with NA's
ids_array[ids_array==ids_fillvalue$value] <- NA
nc_close(nc_trial)

plot(sa_data$X2001.01.15)
#plot(biomass_amazon_mask_pol,add=T, col='red')
plot(bolivia_amazon_pol,add=T, col='red');plot(brazil_amazon_pol,add=T, col='orange');plot(colombia_amazon_pol,add=T, col='black')
plot(ecuador_amazon_pol,add=T, col='purple');plot(freguy_amazon_pol,add=T, col='blue');plot(guyana_amazon_pol,add=T, col='yellow')
plot(peru_amazon_pol,add=T, col='pink');plot(suriname_amazon_pol,add=T, col='white');plot(venezuela_amazon_pol,add=T, col='green')

plot(amazon_north,add=T, col='pink');plot(amazon_south,add=T, col='white');plot(amazon_east,add=T, col='green');plot(amazon_west,add=T, col='red');plot(amazon_central,add=T, col='orange')
#create mask
mask_remap_fun <- function(p) {
SA_extent <- extent(c(-83.75,-31.75,-56.75,13.75)) #use south america extent from INLAND SA
r <- raster(SA_extent)
res(r)<- 0.5 #use resolution to mirror other masks in ilamb/ could be 1x1
values(r) <- 1
crs(r) <- "+proj=longlat +datum=WGS84 +no_defs"
r_amazon <- mask(r, p, updatevalue=0)#create mask from polygon of RAINFOR amazon region
}

r_bolivia_amazon_pol <- mask_remap_fun(bolivia_poly)
writeRaster(r_bolivia_amazon_pol, "mask-BOLIVIAMAZON-remap.nc", overwrite=TRUE, format="CDF",varname='ids',datatype='LOG1S', xname="lon", yname="lat")
r_brazil_amazon_pol <- mask_remap_fun(brazil_amazon_pol)
writeRaster(r_brazil_amazon_pol, "mask-BRAZILAMAZON-remap.nc", overwrite=TRUE, format="CDF",varname='ids',datatype='LOG1S', xname="lon", yname="lat")
r_colombia_amazon_pol <- mask_remap_fun(colombia_amazon_pol)
writeRaster(r_colombia_amazon_pol, "mask-COLOMBIAMAZON-remap.nc", overwrite=TRUE, format="CDF",varname='ids',datatype='LOG1S', xname="lon", yname="lat")
r_ecuador_amazon_pol <- mask_remap_fun(ecuador_amazon_pol)
writeRaster(r_ecuador_amazon_pol, "mask-ECUADORAMAZON-remap.nc", overwrite=TRUE, format="CDF",varname='ids',datatype='LOG1S', xname="lon", yname="lat")
r_freguy_amazon_pol <- mask_remap_fun(freguy_amazon_pol)
writeRaster(r_freguy_amazon_pol, "mask-FGUYAMAZON-remap.nc", overwrite=TRUE, format="CDF",varname='ids',datatype='LOG1S', xname="lon", yname="lat")
r_guyana_amazon_pol <- mask_remap_fun(guyana_amazon_pol)
writeRaster(r_guyana_amazon_pol, "mask-GUYAMAZON-remap.nc", overwrite=TRUE, format="CDF",varname='ids',datatype='LOG1S', xname="lon", yname="lat")
r_peru_amazon_pol <- mask_remap_fun(peru_amazon_pol)
writeRaster(r_peru_amazon_pol, "mask-PERUAMAZON-remap.nc", overwrite=TRUE, format="CDF",varname='ids',datatype='LOG1S', xname="lon", yname="lat")
r_suriname_amazon_pol <- mask_remap_fun(suriname_amazon_pol)
writeRaster(r_suriname_amazon_pol, "mask-SURINAMAZON-remap.nc", overwrite=TRUE, format="CDF",varname='ids',datatype='LOG1S', xname="lon", yname="lat")
r_venezuela_amazon_pol <- mask_remap_fun(venezuela_amazon_pol)
writeRaster(r_venezuela_amazon_pol, "mask-VENEZUELAMAZON-remap.nc", overwrite=TRUE, format="CDF",varname='ids',datatype='LOG1S', xname="lon", yname="lat")

r_amazon_north_pol <- mask_remap_fun(amazon_north)
writeRaster(r_amazon_north_pol, "mask-NORTHAMAZON-remap.nc", overwrite=TRUE, format="CDF",varname='ids',datatype='LOG1S', xname="lon", yname="lat")
r_amazon_east_pol <- mask_remap_fun(amazon_east)
writeRaster(r_amazon_east_pol, "mask-EASTAMAZON-remap.nc", overwrite=TRUE, format="CDF",varname='ids',datatype='LOG1S', xname="lon", yname="lat")
r_amazon_south_pol <- mask_remap_fun(amazon_south)
writeRaster(r_amazon_south_pol, "mask-SOUTHAMAZON-remap.nc", overwrite=TRUE, format="CDF",varname='ids',datatype='LOG1S', xname="lon", yname="lat")
r_amazon_central_pol <- mask_remap_fun(amazon_central)
writeRaster(r_amazon_central_pol, "mask-CENTRALAMAZON-remap.nc", overwrite=TRUE, format="CDF",varname='ids',datatype='LOG1S', xname="lon", yname="lat")
r_amazon_west_pol <- mask_remap_fun(amazon_west)
writeRaster(r_amazon_west_pol, "mask-WESTAMAZON-remap.nc", overwrite=TRUE, format="CDF",varname='ids',datatype='LOG1S', xname="lon", yname="lat")

plot(r_amazon_north_pol)

nc_trial <- nc_open('R:/ILAMB_beta_devel/RAINFOR_leeds_run/masks.notworking/ilamb-mask-BOLIVIAMAZON.nc') #gcel
nc_trial
ids_array <- ncvar_get(nc_trial,'ids')
# dlname <- ncatt_get(ncin,dname,"long_name")
# dunits <- ncatt_get(ncin,dname,"units")
ids_fillvalue <- ncatt_get(nc_trial,'ids',"_FillValue")
dim(ids_array)
# replace netCDF fill values with NA's
ids_array[ids_array==0] <- NA

# get longitude and latitude
mask_lon <- ncvar_get(nc_trial,"lon")
nmask_lon <- dim(mask_lon)
head(mask_lon)

mask_lat <- ncvar_get(nc_trial,"lat")
nmask_lat <- dim(mask_lat)
head(mask_lat)

print(c(nmask_lon,nmask_lat))

nc_close(nc_trial)

bol_mask_data <- raster('R:/ILAMB_beta_devel/RAINFOR_leeds_run/masks.notworking/ilamb-mask-AMAZONIA.nc')
bol_mask_data
plot(bol_mask_data)

nc_southamazon <- nc_open('./ilamb-mask-SOUTHAMAZON.nc') #gcel
nc_southamazon
nc_close(nc_southamazon) #close nc file

southamazon_data <- raster('./ilamb-mask-SOUTHAMAZON.nc')
southamazon_data
plot(southamazon_data)
summary(southamazon_data)

nc_amazon <- nc_open('./mask-AMAZONIA-remap.nc') #gcel
nc_amazon
nc_close(nc_amazon) #close nc file

amazon_data <- raster('./AMAZONIA.nc')
amazon_data
plot(amazon_data)
summary(amazon_data)

nc_brazil <- nc_open('./ilamb-mask-BRAZIL.nc') #gcel
nc_brazil
nc_close(nc_brazil) #close nc file
brazil_data <- raster('./ilamb-mask-BRAZIL.nc')
brazil_data
plot(brazil_data)
summary(brazil_data)

# get longitude and latitude
mask_lon <- ncvar_get(nc_sa,"longitude")
nmask_lon <- dim(mask_lon)
head(mask_lon)

mask_lat <- ncvar_get(nc_sa,"latitude")
nmask_lat <- dim(mask_lat)
head(mask_lat)

print(c(nmask_lon,nmask_lat))

# get temperature
ids_mask_array <- ncvar_get(nc_sa,'Band1')
dlname <- ncatt_get(nc_sa,dname,"long_name")
dunits <- ncatt_get(nc_sa,dname,"units")
fillvalue <- ncatt_get(nc_sa,dname,"_FillValue")
dim(tmp_array)

r_southamerica
r_sa_bol <- crop(trial2, extent(r_southamerica))
