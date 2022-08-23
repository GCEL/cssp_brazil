
###create array#####
# generate lons, lats and set time
lon3 <- as.array(seq(-83.75,-31.75,0.5))
nlon3 <- 105
lat3 <- as.array(seq(-56.75,13.75,0.5))
nlat3 <- 142

# create arrays
# nlon * nlat * nt array
fillvalue <- 1
Band1_array3 <- array(fillvalue, dim=c(nlon3,nlat3))


# path and file name, set dname
ncpath <- "./trial_nc_files/"
ncname <- "SA_bolivia_map"  
ncfname <- paste(ncpath, ncname, ".nc", sep="")
dname <- "Band1"  # note: tmp means temperature (not temporary)

# create and write the netCDF file -- ncdf4 version
# define dimensions
londim <- ncdim_def("lon","degrees_east",as.double(lon3)) 
latdim <- ncdim_def("lat","degrees_north",as.double(lat3)) 
#timedim <- ncdim_def("time",tunits3,as.double(time3))

# define variables
fillvalue <- 1e20
dlname <- "GDAL Band Number 1"
Band1_def <- ncvar_def("ids","",list(londim,latdim),fillvalue,dlname,prec="double")
print(c(nlon3,nlat3))
# create netCDF file and put arrays
ncout <- nc_create(ncfname,list(Band1_def),force_v4=TRUE)
#ncout=NULL
# put variables
ncvar_put(ncout,Band1_def,ids_array,start=c(2,2))
#ncvar_put(ncout,Band1_def,1,start=c(30,))

# put additional attributes into dimension and data variables
ncatt_put(ncout,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout,"lat","axis","Y")
#ncatt_put(ncout,"time","axis","T")

# Get a summary of the created file:
ncout
nc_close(ncout)
nc_sa <- nc_open('./trial_nc_files/SA_bolivia_map.nc') #gcel
nc_sa
nc_close(nc_sa)
sa_data <- raster('./trial_nc_files/SA_bolivia_map.nc')
sa_data
plot(sa_data)
summary(sa_data)