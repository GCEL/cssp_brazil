# Load libraries
library(fields) ; library(ncdf4) ; library(raster)

# Open the file
#cardamom_data = nc_open("G:/ECMWF/ERA5/0.125deg_global/precipitation_daily_mean_200410.nc")
cardamom_data = nc_open("G://CARDAMOM_ILAMB/CARDAMOM_Brazil_1x1_2001_2017_v1.0.nc")
# Read in the precipitation information (kg/m2/s)
npp_flux = ncvar_get(cardamom_data, "NPP_wood_flx")
# Average across the month
npp_flux = apply(npp_flux,c(1,2), mean)
# Read in the lat / long information
latitude = ncvar_get(cardamom_data, "lat")
longitude = ncvar_get(cardamom_data, "lon")
# Tidy up
nc_close(cardamom_data)

# To help visualise the next step plot them up
par(mfrow=c(2,3)) ; plot(latitude) ; plot(longitude) ; image.plot(npp_flux)
# Note the direction of the latitude change, i.e. does it move from north to south or vice versa
# How does this compare with the npp_fluxitation map? They should be in matching orientation, which makes the map upside down but consistent with its lat / long information.

# Flip the latitude axis of both to make them appear the correct way around
latitude = latitude[length(latitude):1]
npp_flux = npp_flux[,dim(npp_flux)[2]:1]
# Now when you plot them they should be the correct way round
image.plot(npp_flux)

# Create the full lat / long grid
latitude = array(rep(latitude, each = length(longitude)), dim=dim(npp_flux))
longitude = array(longitude, dim=dim(npp_flux))
# Check they match the npp_flux
image.plot(latitude) ; image.plot(longitude)

# Create the target raster structure
npp_flux = rasterFromXYZ(data.frame(x = as.vector(longitude),y = as.vector(latitude),z = as.vector(npp_flux)), crs = CRS("+init=epsg:4326"))
par(mfrow=c(1,1));plot(npp_flux)
# You can then use extract or mask function based on other raster or shape file to extract the area of interest.

