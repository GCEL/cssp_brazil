
###
## Load Trendy v11 models and combine in a friendly object matching time periods (2000-2021).
## These data should all be in monthly time step
## 12 Models:
## CARDAMOM, CLASSIC, CML5.0, ISBA-CTRIP, JSBACH, JULES, LPJ-GUESS, LPX-Bern, OCRHIDEE, SDGVM, VISIT-WIES, VISIT
## The following models are regridded from their native resolutions to 1x1 for consistency across the analysis
###

# Load needed libraries
library(ncdf4)
library(fields)
library(raster)
library(compiler)
library(marmap)
library(zoo)
# Load any CARDAMOM libraries needed
source("M://CARDAMOM/CARDAMOM/R_functions/generate_wgs_grid.r")
source("M://CARDAMOM/CARDAMOM/R_functions/determine_lat_long_needed.r")

# Define local function
regrid_func<-function(var1_in, lat_in, long_in, cardamom_ext) {

   # Set flags
   lat_done = FALSE

   # Loop through each timestep in the year
   for (t in seq(1, dim(var1_in)[3])) {
        # Convert to a raster, assuming standad WGS84 grid
        var1 = data.frame(x = as.vector(long_in), y = as.vector(lat_in), z = as.vector(var1_in[,,t]))
        if (length(unique(diff(input_long[,1]))) > 1 | length(unique(diff(input_lat[1,]))) > 1) {
            var1 = griddify(var1, dim(var1_in)[1], dim(var1_in)[2])
        } else {
            var1 = rasterFromXYZ(var1, crs = ("+init=epsg:4326"))
        }

        # Create raster with the target crs (technically this bit is not required)
        target = raster(crs = ("+init=epsg:4326"), ext = extent(var1), resolution = res(var1))
        # Check whether the target and actual analyses have the same CRS
        if (compareCRS(var1,target) == FALSE) {
            # Resample to correct grid
            var1 = resample(var1, target, method="ngb") ; gc() ; removeTmpFiles()
        }
        # Extend the extent of the overall grid to the analysis domain
        var1 = extend(var1,cardamom_ext) 
        # Trim the extent of the overall grid to the analysis domain
        var1 = crop(var1,cardamom_ext)

        # If this is a gridded analysis and the desired CARDAMOM resolution is coarser than the currently provided then aggregate here.
        # Despite creation of a cardamom_ext for a site run do not allow aggragation here as tis will damage the fine resolution datasets
        if (res(var1)[1] != res(cardamom_ext)[1] | res(var1)[2] != res(cardamom_ext)[2]) {

            # Create raster with the target resolution
            target = raster(crs = crs(cardamom_ext), ext = extent(cardamom_ext), resolution = res(cardamom_ext))
            # Resample to correct grid
            var1 = resample(var1, target, method="bilinear") ; gc() ; removeTmpFiles()

        } # Aggrgeate to resolution

        if (lat_done == FALSE) {
            # extract dimension information for the grid, note the axis switching between raster and actual array
            xdim = dim(var1)[2] ; ydim = dim(var1)[1]
            # extract the lat / long information needed
            long = coordinates(var1)[,1] ; lat = coordinates(var1)[,2]
            # restructure into correct orientation
            long = array(long, dim=c(xdim,ydim))
            lat = array(lat, dim=c(xdim,ydim))
        }
        # break out from the rasters into arrays which we can manipulate
        var1 = array(as.vector(unlist(var1)), dim=c(xdim,ydim))

        # vectorise at this time
        if (lat_done == FALSE) {
            var_out = as.vector(var1)
        } else {
            var_out = append(var_out,as.vector(var1))
        }

        # update flag for lat / long load
        if (lat_done == FALSE) {lat_done = TRUE}

   } # Within variable time loop

   # restructure
   var_out = array(var_out, dim=c(xdim,ydim,dim(var1_in)[3]))

   # Return to user
   return(list(var = var_out, lat = lat, long = long))
   
} # end function regrid_func

# Set working directory
setwd("G://Trendy_v11/models/")

# Assuming all grids are the same, therefore have the same area.
# Therefore, we can use CARDAMOM's equivalent for this
load("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/nbe_data_alone/infofile.RData")
year = as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)
year = seq(as.numeric(PROJECT$start_year),as.numeric(PROJECT$end_year)+1, length.out = 1+(length(year)*12))
year = year[-length(year)]
area_m2 = array(PROJECT$area_m2, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
land_fraction = array(PROJECT$landsea, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
vegetated_cover = array(0, dim=c(PROJECT$long_dim,PROJECT$lat_dim))

# Use CARDAMOM's allowed pixels to determine a common set of pixels which have green stuff
vegetated_cover[as.numeric(PROJECT$sites)] = 1

# Flip on the latitude axis, for consistency with how the processed model maps will be put out
area_m2 = area_m2[,PROJECT$lat_dim:1]
land_fraction = land_fraction[,PROJECT$lat_dim:1]
vegetated_cover = vegetated_cover[,PROJECT$lat_dim:1]

# Determine the cardamom_ext needed to allow for regridding 
output = determine_lat_long_needed(PROJECT$latitude,PROJECT$longitude,PROJECT$resolution,PROJECT$grid_type,PROJECT$waterpixels)
cardamom_ext = output$cardamom_ext

# Set number of models
nos_models = 12
model_list = c("CARDAMOM","CLASSIC","CLM5.0","ISBA-CTRIP","JSBACH","JULES","LPJ-GUESS","LPX-Bern","OCRHIDEE","SDGVM","VISIT","VISIT-NIES")
if (length(model_list) != nos_models) {stop("Error specifying nos_models and model_list")}

# Define model list variable
trendy = list(model_list = model_list, year = year,
              area_m2 = area_m2, land_fraction = land_fraction, vegetated_cover = vegetated_cover,
              latitude = output$obs_lat_grid, longitude = output$obs_long_grid,
              lai_m2m2     = array(NA, dim = c(PROJECT$long_dim,PROJECT$lat_dim,22*12)),
              gpp_gCm2day  = array(NA, dim = c(PROJECT$long_dim,PROJECT$lat_dim,22*12)),
              ra_gCm2day   = array(NA, dim = c(PROJECT$long_dim,PROJECT$lat_dim,22*12)),
              rh_gCm2day   = array(NA, dim = c(PROJECT$long_dim,PROJECT$lat_dim,22*12)),
              fire_gCm2day = array(NA, dim = c(PROJECT$long_dim,PROJECT$lat_dim,22*12)),
              nbp_gCm2day  = array(NA, dim = c(PROJECT$long_dim,PROJECT$lat_dim,22*12)),          
              nbe_gCm2day  = array(NA, dim = c(PROJECT$long_dim,PROJECT$lat_dim,22*12)),          
              soil_gCm2    = array(NA, dim = c(PROJECT$long_dim,PROJECT$lat_dim,22*12)),              
              wood_gCm2    = array(NA, dim = c(PROJECT$long_dim,PROJECT$lat_dim,22*12)),              
              veg_gCm2     = array(NA, dim = c(PROJECT$long_dim,PROJECT$lat_dim,22*12)))

# tidy
rm(PROJECT,output)
###
## Output to combined netcdf file

output_prefix = "Trendyv11_ensemble_" # follow with "_"
output_suffix = "" # begin with "_"
output_dir = "./individual/"


for (i in model_list){
  if (i == "CARDAMOM") {
    ###
    ## CARDAMOM
    
    # model_counter = 1
    
    # LAI (m2m2)
    input = nc_open("./CARDAMOM/CARDAMOM_S3_lai.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "lai")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180) 
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$lai_m2m2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # GPP (kgC/m2/s -> gC/m2/d)
    input = nc_open("./CARDAMOM/CARDAMOM_S3_gpp.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "gpp")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180) 
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$gpp_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Ra (kgC/m2/s -> gC/m2/d)
    input = nc_open("./CARDAMOM/CARDAMOM_S3_ra.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "ra")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180) 
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$ra_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Rh (kgC/m2/s -> gC/m2/d)
    input = nc_open("./CARDAMOM/CARDAMOM_S3_rh.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "rh")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180) 
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$rh_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Fire (kgC/m2/s -> gC/m2/d)
    input = nc_open("./CARDAMOM/CARDAMOM_S3_fFire.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "fFire")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180) 
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$fire_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # NBP (kgC/m2/s -> gC/m2/d)
    input = nc_open("./CARDAMOM/CARDAMOM_S3_nbp.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "nbp")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180) 
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$nbp_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # NBE (kgC/m2/s -> gC/m2/d)
    input = nc_open("./CARDAMOM/CARDAMOM_S3_nbe.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "nbe")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180) 
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$nbe_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Veg (kgC/m2 -> gC/m2)
    input = nc_open("./CARDAMOM/CARDAMOM_S3_cVeg.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "cVeg")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180) 
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$veg_gCm2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Soil (kgC/m2 -> gC/m2)
    input = nc_open("./CARDAMOM/CARDAMOM_S3_cSoil.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "cSoil")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180) 
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$soil_gCm2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Wood (kgC/m2 -> gC/m2)
    input = nc_open("./CARDAMOM/CARDAMOM_S3_cWoodTotal.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "cWoodTotal")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180) 
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$wood_gCm2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    ## define dimension
    lat_dimen <- ncdim_def( "lat", units="degree north (-90->90)", trendy$latitude[1,] )
    long_dimen <- ncdim_def( "lon", units="degree east (-180->180)", trendy$longitude[,1] )
    time_dimen <- ncdim_def( "time", units="", trendy$year)
    # model_dimen <- ncdim_def( "model", units="", c(1:length(trendy$model_list)))
    dimnchar <- ncdim_def("nchar", "", 1:18, create_dimvar=FALSE ) # Maximum number of characters used in string vector
    
    ## define output variable
    # var0 = ncvar_def("model_list", units = "-", longname = "Model name", 
    #                  dim=list(dimnchar,model_dimen), prec="char", compression = 9)
    var0 = ncvar_def("area", units = "m2", longname = paste("Pixel area",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    var1 = ncvar_def("land_fraction", units = "1", longname = paste("Fraction of pixel which is land",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    var2 = ncvar_def("vegetated_cover", units = "1", longname = paste("Pixels determined to have vegetation in CARDAMOM grid",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    
    # Define the output file name
    output_name = paste(output_dir,output_prefix,i,"_trendy_output",output_suffix,".nc", sep="")
    # Delete if the file currently exists
    if (file.exists(output_name)) {file.remove(output_name)}
    # Create the empty file space
    new_file=nc_create(filename=output_name, vars=list(var0,var1,var2), force_v4 = TRUE)
    # Load first variable into the file
    # # Model list
    # ncvar_put(new_file, var0, trendy$model_list)
    # Grid area 
    ncvar_put(new_file, var0, trendy$area_m2)
    # Land fraction
    ncvar_put(new_file, var1, trendy$land_fraction)
    # Land fraction
    ncvar_put(new_file, var2, trendy$vegetated_cover)
    # Close the existing file to ensure its written to file
    nc_close(new_file)
    
    ## Re-open the file so that we can add to it a variable at a time
    new_file <- nc_open( output_name, write=TRUE )
    
    # Leaf area index
    var_new  = ncvar_def("LAI", unit="m2.m-2", longname = "Leaf Area Index", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$lai_m2m2)
    
    # Gross Primary Productivity
    var_new  = ncvar_def("GPP", unit="gC.m-2.d-1", longname = "Gross Primary Productivity", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$gpp_gCm2day)
    
    # Autotrophic respiration
    var_new  = ncvar_def("Ra", unit="gC.m-2.d-1", longname = "Autotrophic respiration", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$ra_gCm2day)
    
    # Heterotrophic respiration
    var_new  = ncvar_def("Rh", unit="gC.m-2.d-1", longname = "Heterotrophic respiration", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$rh_gCm2day)
    
    # Fire
    var_new  = ncvar_def("FIRE", unit="gC.m-2.d-1", longname = "Fire C emissions", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$fire_gCm2day)
    
    # NBP
    var_new  = ncvar_def("NBP", unit="gC.m-2.d-1", longname = "Net Biome Productivity", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$nbp_gCm2day)
    
    # NBE
    var_new  = ncvar_def("NBE", unit="gC.m-2.d-1", longname = "Net Biome Exchange", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$nbe_gCm2day)
    
    # Soil
    var_new  = ncvar_def("SOIL", unit="gC.m-2", longname = "Soil C stocks", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$soil_gCm2)
    
    # Vegetation
    var_new  = ncvar_def("TOT", unit="gC.m-2", longname = "Vegetation C stocks", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$veg_gCm2)
    
    # Wood
    var_new  = ncvar_def("WOOD", unit="gC.m-2", longname = "Wood C stocks", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$wood_gCm2)
    
    # Close file
    nc_close(new_file)
  }
  if (i == "CLASSIC") {
    ###
    ## CLASSIC
    
    # model_counter = 2
    
    # LAI (m2m2)
    input = nc_open("./CLASSIC/CLASSIC_S3_lai.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "lai")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$lai_m2m2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # GPP (kgC/m2/s -> gC/m2/d)
    input = nc_open("./CLASSIC/CLASSIC_S3_gpp.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "gpp")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$gpp_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Ra (kgC/m2/s -> gC/m2/d)
    input = nc_open("./CLASSIC/CLASSIC_S3_ra.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "ra")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$ra_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Rh (kgC/m2/s -> gC/m2/d)
    input = nc_open("./CLASSIC/CLASSIC_S3_rh.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "rh")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$rh_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Fire (kgC/m2/s -> gC/m2/d)
    input = nc_open("./CLASSIC/CLASSIC_S3_fFire.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "fFire")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$fire_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # NBP (kgC/m2/s -> gC/m2/d)
    input = nc_open("./CLASSIC/CLASSIC_S3_nbp.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "nbp")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data_nbp = input_data * 1e3 * 86400 #save for NBE calculation
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$nbp_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # fLuc (kgC/m2/s -> gC/m2/d)
    input = nc_open("./CLASSIC/CLASSIC_S3_fLuc.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "fLuc")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180) 
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data_fluc = input_data * 1e3 * 86400
    # # Tidy
    # nc_close(input) ; rm(input_data) ; gc()
    
    # NBE (kgC/m2/s -> gC/m2/d)
    # Unit convertion
    input_data_nbe = input_data_nbp + input_data_fluc
    input_data_nbe = input_data_nbe * -1
    # Begin regridding
    input_data_nbe = regrid_func(input_data_nbe,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$nbe_gCm2day[,,] = input_data_nbe$var
    # Tidy
    nc_close(input_data_nbp) ; rm(input_data_fluc); rm(input_data); rm(input_data_nbe) ; gc()
    
    # Veg (kgC/m2 -> gC/m2)
    input = nc_open("./CLASSIC/CLASSIC_S3_cVeg.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "cVeg")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$veg_gCm2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Soil (kgC/m2 -> gC/m2)
    input = nc_open("./CLASSIC/CLASSIC_S3_cSoil.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "cSoil")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$soil_gCm2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Wood (kgC/m2 -> gC/m2)
    input = nc_open("./CLASSIC/CLASSIC_S3_cWood.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "cWood")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180) 
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$wood_gCm2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    ## define dimension
    lat_dimen <- ncdim_def( "lat", units="degree north (-90->90)", trendy$latitude[1,] )
    long_dimen <- ncdim_def( "lon", units="degree east (-180->180)", trendy$longitude[,1] )
    time_dimen <- ncdim_def( "time", units="", trendy$year)
    # model_dimen <- ncdim_def( "model", units="", c(1:length(trendy$model_list)))
    dimnchar <- ncdim_def("nchar", "", 1:18, create_dimvar=FALSE ) # Maximum number of characters used in string vector
    
    ## define output variable
    # var0 = ncvar_def("model_list", units = "-", longname = "Model name", 
    #                  dim=list(dimnchar,model_dimen), prec="char", compression = 9)
    var0 = ncvar_def("area", units = "m2", longname = paste("Pixel area",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    var1 = ncvar_def("land_fraction", units = "1", longname = paste("Fraction of pixel which is land",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    var2 = ncvar_def("vegetated_cover", units = "1", longname = paste("Pixels determined to have vegetation in CARDAMOM grid",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    
    # Define the output file name
    output_name = paste(output_dir,output_prefix,i,"_trendy_output",output_suffix,".nc", sep="")
    # Delete if the file currently exists
    if (file.exists(output_name)) {file.remove(output_name)}
    # Create the empty file space
    new_file=nc_create(filename=output_name, vars=list(var0,var1,var2), force_v4 = TRUE)
    # Load first variable into the file
    # # Model list
    # ncvar_put(new_file, var0, trendy$model_list)
    # Grid area 
    ncvar_put(new_file, var0, trendy$area_m2)
    # Land fraction
    ncvar_put(new_file, var1, trendy$land_fraction)
    # Land fraction
    ncvar_put(new_file, var2, trendy$vegetated_cover)
    # Close the existing file to ensure its written to file
    nc_close(new_file)
    
    ## Re-open the file so that we can add to it a variable at a time
    new_file <- nc_open( output_name, write=TRUE )
    
    # Leaf area index
    var_new  = ncvar_def("LAI", unit="m2.m-2", longname = "Leaf Area Index", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$lai_m2m2)
    
    # Gross Primary Productivity
    var_new  = ncvar_def("GPP", unit="gC.m-2.d-1", longname = "Gross Primary Productivity", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$gpp_gCm2day)
    
    # Autotrophic respiration
    var_new  = ncvar_def("Ra", unit="gC.m-2.d-1", longname = "Autotrophic respiration", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$ra_gCm2day)
    
    # Heterotrophic respiration
    var_new  = ncvar_def("Rh", unit="gC.m-2.d-1", longname = "Heterotrophic respiration", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$rh_gCm2day)
    
    # Fire
    var_new  = ncvar_def("FIRE", unit="gC.m-2.d-1", longname = "Fire C emissions", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$fire_gCm2day)
    
    # NBP
    var_new  = ncvar_def("NBP", unit="gC.m-2.d-1", longname = "Net Biome Productivity", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$nbp_gCm2day)
    
    # NBE
    var_new  = ncvar_def("NBE", unit="gC.m-2.d-1", longname = "Net Biome Exchange", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$nbe_gCm2day)
    
    # Soil
    var_new  = ncvar_def("SOIL", unit="gC.m-2", longname = "Soil C stocks", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$soil_gCm2)
    
    # Vegetation
    var_new  = ncvar_def("TOT", unit="gC.m-2", longname = "Vegetation C stocks", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$veg_gCm2)
    
    # Wood
    var_new  = ncvar_def("WOOD", unit="gC.m-2", longname = "Wood C stocks", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$wood_gCm2)
    
    # Close file
    nc_close(new_file)
    
  }
  if (i == "CLM5.0") {
    ###
    ## CLM5.0
    
    # model_counter = 3
    
    # LAI (m2m2)
    input = nc_open("./CLM5.0/CLM5.0_S3_lai.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "lai")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$lai_m2m2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # GPP (kgC/m2/s -> gC/m2/d)
    input = nc_open("./CLM5.0/CLM5.0_S3_gpp.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "gpp")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$gpp_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Ra (kgC/m2/s -> gC/m2/d)
    input = nc_open("./CLM5.0/CLM5.0_S3_ra.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "ra")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$ra_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Rh (kgC/m2/s -> gC/m2/d)
    input = nc_open("./CLM5.0/CLM5.0_S3_rh.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "rh")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$rh_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Fire (kgC/m2/s -> gC/m2/d)
    input = nc_open("./CLM5.0/CLM5.0_S3_fFire.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "fFire")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$fire_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # NBP (kgC/m2/s -> gC/m2/d)
    input = nc_open("./CLM5.0/CLM5.0_S3_nbp.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "nbp")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data_nbp = input_data * 1e3 * 86400 #save for NBE calculation
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$nbp_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # fLuc (kgC/m2/s -> gC/m2/d)
    input = nc_open("./CLM5.0/CLM5.0_S3_fLuc.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "fLuc")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180) 
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data_fluc = input_data * 1e3 * 86400
    # Tidy
    rm(input_data) ; gc()
    
    # NBE (kgC/m2/s -> gC/m2/d)
    # Unit convertion
    input_data = (input_data_nbp+input_data_fluc)*-1
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$nbe_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input_data_nbp) ; rm(input); rm(input_data_fluc); rm(input_data) ; gc()
    
    # Veg (kgC/m2 -> gC/m2)
    input = nc_open("./CLM5.0/CLM5.0_S3_cVeg.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "cVeg")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$veg_gCm2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Soil (kgC/m2 -> gC/m2)
    input = nc_open("./CLM5.0/CLM5.0_S3_cSoil.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "cSoil")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$soil_gCm2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Wood (kgC/m2 -> gC/m2)
    input = nc_open("./CLM5.0/CLM5.0_S3_cWood.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "cWood")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180) 
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$wood_gCm2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    ## define dimension
    lat_dimen <- ncdim_def( "lat", units="degree north (-90->90)", trendy$latitude[1,] )
    long_dimen <- ncdim_def( "lon", units="degree east (-180->180)", trendy$longitude[,1] )
    time_dimen <- ncdim_def( "time", units="", trendy$year)
    # model_dimen <- ncdim_def( "model", units="", c(1:length(trendy$model_list)))
    dimnchar <- ncdim_def("nchar", "", 1:18, create_dimvar=FALSE ) # Maximum number of characters used in string vector
    
    ## define output variable
    # var0 = ncvar_def("model_list", units = "-", longname = "Model name", 
    #                  dim=list(dimnchar,model_dimen), prec="char", compression = 9)
    var0 = ncvar_def("area", units = "m2", longname = paste("Pixel area",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    var1 = ncvar_def("land_fraction", units = "1", longname = paste("Fraction of pixel which is land",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    var2 = ncvar_def("vegetated_cover", units = "1", longname = paste("Pixels determined to have vegetation in CARDAMOM grid",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    
    # Define the output file name
    output_name = paste(output_dir,output_prefix,i,"_trendy_output",output_suffix,".nc", sep="")
    # Delete if the file currently exists
    if (file.exists(output_name)) {file.remove(output_name)}
    # Create the empty file space
    new_file=nc_create(filename=output_name, vars=list(var0,var1,var2), force_v4 = TRUE)
    # Load first variable into the file
    # # Model list
    # ncvar_put(new_file, var0, trendy$model_list)
    # Grid area 
    ncvar_put(new_file, var0, trendy$area_m2)
    # Land fraction
    ncvar_put(new_file, var1, trendy$land_fraction)
    # Land fraction
    ncvar_put(new_file, var2, trendy$vegetated_cover)
    # Close the existing file to ensure its written to file
    nc_close(new_file)
    
    ## Re-open the file so that we can add to it a variable at a time
    new_file <- nc_open( output_name, write=TRUE )
    
    # Leaf area index
    var_new  = ncvar_def("LAI", unit="m2.m-2", longname = "Leaf Area Index", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$lai_m2m2)
    
    # Gross Primary Productivity
    var_new  = ncvar_def("GPP", unit="gC.m-2.d-1", longname = "Gross Primary Productivity", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$gpp_gCm2day)
    
    # Autotrophic respiration
    var_new  = ncvar_def("Ra", unit="gC.m-2.d-1", longname = "Autotrophic respiration", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$ra_gCm2day)
    
    # Heterotrophic respiration
    var_new  = ncvar_def("Rh", unit="gC.m-2.d-1", longname = "Heterotrophic respiration", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$rh_gCm2day)
    
    # Fire
    var_new  = ncvar_def("FIRE", unit="gC.m-2.d-1", longname = "Fire C emissions", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$fire_gCm2day)
    
    # NBP
    var_new  = ncvar_def("NBP", unit="gC.m-2.d-1", longname = "Net Biome Productivity", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$nbp_gCm2day)
    
    # NBE
    var_new  = ncvar_def("NBE", unit="gC.m-2.d-1", longname = "Net Biome Exchange", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$nbe_gCm2day)
    
    # Soil
    var_new  = ncvar_def("SOIL", unit="gC.m-2", longname = "Soil C stocks", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$soil_gCm2)
    
    # Vegetation
    var_new  = ncvar_def("TOT", unit="gC.m-2", longname = "Vegetation C stocks", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$veg_gCm2)
    
    # Wood
    var_new  = ncvar_def("WOOD", unit="gC.m-2", longname = "Wood C stocks", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$wood_gCm2)
    
    # Close file
    nc_close(new_file)
    
  }
  if (i == "ISBA-CTRIP") {
    ###
    ## ISBA-CTRIP
    
    # model_counter = 4
    
    # LAI (m2m2)
    input = nc_open("./ISBA-CTRIP/ISBA-CTRIP_S3_lai.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "lai")
    input_lat = ncvar_get(input,"lat_FULL") ; input_long = ncvar_get(input, "lon_FULL")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$lai_m2m2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # GPP (kgC/m2/s -> gC/m2/d)
    input = nc_open("./ISBA-CTRIP/ISBA-CTRIP_S3_gpp.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "gpp")
    input_lat = ncvar_get(input,"lat_FULL") ; input_long = ncvar_get(input, "lon_FULL")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$gpp_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Ra (kgC/m2/s -> gC/m2/d)
    input = nc_open("./ISBA-CTRIP/ISBA-CTRIP_S3_ra.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "ra")
    input_lat = ncvar_get(input,"lat_FULL") ; input_long = ncvar_get(input, "lon_FULL")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$ra_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Rh (kgC/m2/s -> gC/m2/d)
    input = nc_open("./ISBA-CTRIP/ISBA-CTRIP_S3_rh.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "rh")
    input_lat = ncvar_get(input,"lat_FULL") ; input_long = ncvar_get(input, "lon_FULL")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$rh_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Fire (kgC/m2/s -> gC/m2/d)
    input = nc_open("./ISBA-CTRIP/ISBA-CTRIP_S3_fFire.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "fFire")
    input_lat = ncvar_get(input,"lat_FULL") ; input_long = ncvar_get(input, "lon_FULL")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$fire_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # NBP (kgC/m2/s -> gC/m2/d)
    input = nc_open("./ISBA-CTRIP/ISBA-CTRIP_S3_nbp.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "nbp")
    input_lat = ncvar_get(input,"lat_FULL") ; input_long = ncvar_get(input, "lon_FULL")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data_nbp = input_data * 1e3 * 86400 #save for NBE calculation
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$nbp_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # fLuc (kgC/m2/s -> gC/m2/d)
    input = nc_open("./ISBA-CTRIP/ISBA-CTRIP_S3_fLuc.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "fLuc")
    input_lat = ncvar_get(input,"lat_FULL") ; input_long = ncvar_get(input, "lon_FULL")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180) 
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data_fluc = input_data * 1e3 * 86400
    # Tidy
    rm(input_data) ; gc()
    
    # NBE (kgC/m2/s -> gC/m2/d)
    # Unit convertion
    input_data = (input_data_nbp+input_data_fluc)*-1
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$nbe_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Veg (kgC/m2 -> gC/m2)
    input = nc_open("./ISBA-CTRIP/ISBA-CTRIP_S3_cVeg.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "cVeg")
    input_lat = ncvar_get(input,"lat_FULL") ; input_long = ncvar_get(input, "lon_FULL")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$veg_gCm2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Soil (kgC/m2 -> gC/m2)
    input = nc_open("./ISBA-CTRIP/ISBA-CTRIP_S3_cSoil.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "cSoil")
    input_lat = ncvar_get(input,"lat_FULL") ; input_long = ncvar_get(input, "lon_FULL")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$soil_gCm2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Wood (kgC/m2 -> gC/m2)
    input = nc_open("./ISBA-CTRIP/ISBA-CTRIP_S3_cWood.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "cWood")
    input_lat = ncvar_get(input,"lat_FULL") ; input_long = ncvar_get(input, "lon_FULL")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180) 
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$wood_gCm2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    ## define dimension
    lat_dimen <- ncdim_def( "lat", units="degree north (-90->90)", trendy$latitude[1,] )
    long_dimen <- ncdim_def( "lon", units="degree east (-180->180)", trendy$longitude[,1] )
    time_dimen <- ncdim_def( "time", units="", trendy$year)
    # model_dimen <- ncdim_def( "model", units="", c(1:length(trendy$model_list)))
    dimnchar <- ncdim_def("nchar", "", 1:18, create_dimvar=FALSE ) # Maximum number of characters used in string vector
    
    ## define output variable
    # var0 = ncvar_def("model_list", units = "-", longname = "Model name", 
    #                  dim=list(dimnchar,model_dimen), prec="char", compression = 9)
    var0 = ncvar_def("area", units = "m2", longname = paste("Pixel area",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    var1 = ncvar_def("land_fraction", units = "1", longname = paste("Fraction of pixel which is land",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    var2 = ncvar_def("vegetated_cover", units = "1", longname = paste("Pixels determined to have vegetation in CARDAMOM grid",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    
    # Define the output file name
    output_name = paste(output_dir,output_prefix,i,"_trendy_output",output_suffix,".nc", sep="")
    # Delete if the file currently exists
    if (file.exists(output_name)) {file.remove(output_name)}
    # Create the empty file space
    new_file=nc_create(filename=output_name, vars=list(var0,var1,var2), force_v4 = TRUE)
    # Load first variable into the file
    # # Model list
    # ncvar_put(new_file, var0, trendy$model_list)
    # Grid area 
    ncvar_put(new_file, var0, trendy$area_m2)
    # Land fraction
    ncvar_put(new_file, var1, trendy$land_fraction)
    # Land fraction
    ncvar_put(new_file, var2, trendy$vegetated_cover)
    # Close the existing file to ensure its written to file
    nc_close(new_file)
    
    ## Re-open the file so that we can add to it a variable at a time
    new_file <- nc_open( output_name, write=TRUE )
    
    # Leaf area index
    var_new  = ncvar_def("LAI", unit="m2.m-2", longname = "Leaf Area Index", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$lai_m2m2)
    
    # Gross Primary Productivity
    var_new  = ncvar_def("GPP", unit="gC.m-2.d-1", longname = "Gross Primary Productivity", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$gpp_gCm2day)
    
    # Autotrophic respiration
    var_new  = ncvar_def("Ra", unit="gC.m-2.d-1", longname = "Autotrophic respiration", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$ra_gCm2day)
    
    # Heterotrophic respiration
    var_new  = ncvar_def("Rh", unit="gC.m-2.d-1", longname = "Heterotrophic respiration", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$rh_gCm2day)
    
    # Fire
    var_new  = ncvar_def("FIRE", unit="gC.m-2.d-1", longname = "Fire C emissions", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$fire_gCm2day)
    
    # NBP
    var_new  = ncvar_def("NBP", unit="gC.m-2.d-1", longname = "Net Biome Productivity", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$nbp_gCm2day)
    
    # NBE
    var_new  = ncvar_def("NBE", unit="gC.m-2.d-1", longname = "Net Biome Exchange", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$nbe_gCm2day)
    
    # Soil
    var_new  = ncvar_def("SOIL", unit="gC.m-2", longname = "Soil C stocks", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$soil_gCm2)
    
    # Vegetation
    var_new  = ncvar_def("TOT", unit="gC.m-2", longname = "Vegetation C stocks", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$veg_gCm2)
    
    # Wood
    var_new  = ncvar_def("WOOD", unit="gC.m-2", longname = "Wood C stocks", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$wood_gCm2)
    
    # Close file
    nc_close(new_file)
    
  }
  if (i == "JSBACH") {
    ###
    ## JSBACH
    
    # model_counter = 5
    
    # LAI (m2m2)
    input = nc_open("./JSBACH/JSBACH_S3_lai.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "lai")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$lai_m2m2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # GPP (kgC/m2/s -> gC/m2/d)
    input = nc_open("./JSBACH/JSBACH_S3_gpp.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "gpp")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$gpp_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Ra (kgC/m2/s -> gC/m2/d)
    input = nc_open("./JSBACH/JSBACH_S3_ra.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "ra")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$ra_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Rh (kgC/m2/s -> gC/m2/d)
    input = nc_open("./JSBACH/JSBACH_S3_rh.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "rh")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$rh_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Fire (kgC/m2/s -> gC/m2/d)
    input = nc_open("./JSBACH/JSBACH_S3_fFire.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "fFire")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$fire_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # NBP (kgC/m2/s -> gC/m2/d)
    input = nc_open("./JSBACH/JSBACH_S3_nbp.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "nbp")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data_nbp = input_data * 1e3 * 86400 #save for NBE calculation
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$nbp_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # fLuc (kgC/m2/s -> gC/m2/d)
    input = nc_open("./JSBACH/JSBACH_S3_fLuc.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "fLuc")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180) 
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data_fluc = input_data * 1e3 * 86400
    # Tidy
    rm(input_data) ; gc()
    
    # NBE (kgC/m2/s -> gC/m2/d)
    # Unit convertion
    input_data = (input_data_nbp+input_data_fluc)*-1
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$nbe_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Veg (kgC/m2 -> gC/m2)
    input = nc_open("./JSBACH/JSBACH_S3_cVeg.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "cVeg")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$veg_gCm2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Soil (kgC/m2 -> gC/m2)
    input = nc_open("./JSBACH/JSBACH_S3_cSoil.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "cSoil")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$soil_gCm2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    ## define dimension
    lat_dimen <- ncdim_def( "lat", units="degree north (-90->90)", trendy$latitude[1,] )
    long_dimen <- ncdim_def( "lon", units="degree east (-180->180)", trendy$longitude[,1] )
    time_dimen <- ncdim_def( "time", units="", trendy$year)
    # model_dimen <- ncdim_def( "model", units="", c(1:length(trendy$model_list)))
    dimnchar <- ncdim_def("nchar", "", 1:18, create_dimvar=FALSE ) # Maximum number of characters used in string vector
    
    ## define output variable
    # var0 = ncvar_def("model_list", units = "-", longname = "Model name", 
    #                  dim=list(dimnchar,model_dimen), prec="char", compression = 9)
    var0 = ncvar_def("area", units = "m2", longname = paste("Pixel area",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    var1 = ncvar_def("land_fraction", units = "1", longname = paste("Fraction of pixel which is land",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    var2 = ncvar_def("vegetated_cover", units = "1", longname = paste("Pixels determined to have vegetation in CARDAMOM grid",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    
    # Define the output file name
    output_name = paste(output_dir,output_prefix,i,"_trendy_output",output_suffix,".nc", sep="")
    # Delete if the file currently exists
    if (file.exists(output_name)) {file.remove(output_name)}
    # Create the empty file space
    new_file=nc_create(filename=output_name, vars=list(var0,var1,var2), force_v4 = TRUE)
    # Load first variable into the file
    # # Model list
    # ncvar_put(new_file, var0, trendy$model_list)
    # Grid area 
    ncvar_put(new_file, var0, trendy$area_m2)
    # Land fraction
    ncvar_put(new_file, var1, trendy$land_fraction)
    # Land fraction
    ncvar_put(new_file, var2, trendy$vegetated_cover)
    # Close the existing file to ensure its written to file
    nc_close(new_file)
    
    ## Re-open the file so that we can add to it a variable at a time
    new_file <- nc_open( output_name, write=TRUE )
    
    # Leaf area index
    var_new  = ncvar_def("LAI", unit="m2.m-2", longname = "Leaf Area Index", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$lai_m2m2)
    
    # Gross Primary Productivity
    var_new  = ncvar_def("GPP", unit="gC.m-2.d-1", longname = "Gross Primary Productivity", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$gpp_gCm2day)
    
    # Autotrophic respiration
    var_new  = ncvar_def("Ra", unit="gC.m-2.d-1", longname = "Autotrophic respiration", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$ra_gCm2day)
    
    # Heterotrophic respiration
    var_new  = ncvar_def("Rh", unit="gC.m-2.d-1", longname = "Heterotrophic respiration", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$rh_gCm2day)
    
    # Fire
    var_new  = ncvar_def("FIRE", unit="gC.m-2.d-1", longname = "Fire C emissions", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$fire_gCm2day)
    
    # NBP
    var_new  = ncvar_def("NBP", unit="gC.m-2.d-1", longname = "Net Biome Productivity", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$nbp_gCm2day)
    
    # NBE
    var_new  = ncvar_def("NBE", unit="gC.m-2.d-1", longname = "Net Biome Exchange", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$nbe_gCm2day)
    
    # Soil
    var_new  = ncvar_def("SOIL", unit="gC.m-2", longname = "Soil C stocks", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$soil_gCm2)
    
    # Vegetation
    var_new  = ncvar_def("TOT", unit="gC.m-2", longname = "Vegetation C stocks", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$veg_gCm2)
    
    # Close file
    nc_close(new_file)
    
  }
  if (i == "JULES") {
    ###
    ## JULES
    
    # model_counter = 6
    
    # LAI (m2m2)
    input = nc_open("./JULES/JULES_S3_lai.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "lai")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$lai_m2m2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # GPP (kgC/m2/s -> gC/m2/d)
    input = nc_open("./JULES/JULES_S3_gpp.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "gpp")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$gpp_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Ra (kgC/m2/s -> gC/m2/d)
    input = nc_open("./JULES/JULES_S3_ra.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "ra")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$ra_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Rh (kgC/m2/s -> gC/m2/d)
    input = nc_open("./JULES/JULES_S3_rh.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "rh")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$rh_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Fire (kgC/m2/s -> gC/m2/d)
    input = nc_open("./JULES/JULES_S3_fFire.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "fFire")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$fire_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # NBP (kgC/m2/s -> gC/m2/d)
    input = nc_open("./JULES/JULES_S3_nbp.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "nbp")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data_nbp = input_data * 1e3 * 86400 #save for NBE calculation
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$nbp_gCm2day[,,] = input_data$var
    
    # NBE (kgC/m2/s -> gC/m2/d)
    # Unit convertion
    input_data_nbp = input_data_nbp * -1
    # Begin regridding
    input_data_nbp = regrid_func(input_data_nbp,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$nbe_gCm2day[,,] = input_data_nbp$var
    # Tidy
    nc_close(input); rm(input_data) ; rm(input_data_nbp) ; gc()
    
    # Veg (kgC/m2 -> gC/m2)
    input = nc_open("./JULES/JULES_S3_cVeg.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "cVeg")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$veg_gCm2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Soil (kgC/m2 -> gC/m2)
    input = nc_open("./JULES/JULES_S3_cSoil.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "cSoil")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$soil_gCm2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Wood (kgC/m2 -> gC/m2)
    input = nc_open("./JULES/JULES_S3_cWood.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "cWood")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180) 
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$wood_gCm2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    ## define dimension
    lat_dimen <- ncdim_def( "lat", units="degree north (-90->90)", trendy$latitude[1,] )
    long_dimen <- ncdim_def( "lon", units="degree east (-180->180)", trendy$longitude[,1] )
    time_dimen <- ncdim_def( "time", units="", trendy$year)
    # model_dimen <- ncdim_def( "model", units="", c(1:length(trendy$model_list)))
    dimnchar <- ncdim_def("nchar", "", 1:18, create_dimvar=FALSE ) # Maximum number of characters used in string vector
    
    ## define output variable
    # var0 = ncvar_def("model_list", units = "-", longname = "Model name", 
    #                  dim=list(dimnchar,model_dimen), prec="char", compression = 9)
    var0 = ncvar_def("area", units = "m2", longname = paste("Pixel area",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    var1 = ncvar_def("land_fraction", units = "1", longname = paste("Fraction of pixel which is land",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    var2 = ncvar_def("vegetated_cover", units = "1", longname = paste("Pixels determined to have vegetation in CARDAMOM grid",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    
    # Define the output file name
    output_name = paste(output_dir,output_prefix,i,"_trendy_output",output_suffix,".nc", sep="")
    # Delete if the file currently exists
    if (file.exists(output_name)) {file.remove(output_name)}
    # Create the empty file space
    new_file=nc_create(filename=output_name, vars=list(var0,var1,var2), force_v4 = TRUE)
    # Load first variable into the file
    # # Model list
    # ncvar_put(new_file, var0, trendy$model_list)
    # Grid area 
    ncvar_put(new_file, var0, trendy$area_m2)
    # Land fraction
    ncvar_put(new_file, var1, trendy$land_fraction)
    # Land fraction
    ncvar_put(new_file, var2, trendy$vegetated_cover)
    # Close the existing file to ensure its written to file
    nc_close(new_file)
    
    ## Re-open the file so that we can add to it a variable at a time
    new_file <- nc_open( output_name, write=TRUE )
    
    # Leaf area index
    var_new  = ncvar_def("LAI", unit="m2.m-2", longname = "Leaf Area Index", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$lai_m2m2)
    
    # Gross Primary Productivity
    var_new  = ncvar_def("GPP", unit="gC.m-2.d-1", longname = "Gross Primary Productivity", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$gpp_gCm2day)
    
    # Autotrophic respiration
    var_new  = ncvar_def("Ra", unit="gC.m-2.d-1", longname = "Autotrophic respiration", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$ra_gCm2day)
    
    # Heterotrophic respiration
    var_new  = ncvar_def("Rh", unit="gC.m-2.d-1", longname = "Heterotrophic respiration", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$rh_gCm2day)
    
    # Fire
    var_new  = ncvar_def("FIRE", unit="gC.m-2.d-1", longname = "Fire C emissions", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$fire_gCm2day)
    
    # NBP
    var_new  = ncvar_def("NBP", unit="gC.m-2.d-1", longname = "Net Biome Productivity", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$nbp_gCm2day)
    
    # NBE
    var_new  = ncvar_def("NBE", unit="gC.m-2.d-1", longname = "Net Biome Exchange", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$nbe_gCm2day)
    
    # Soil
    var_new  = ncvar_def("SOIL", unit="gC.m-2", longname = "Soil C stocks", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$soil_gCm2)
    
    # Vegetation
    var_new  = ncvar_def("TOT", unit="gC.m-2", longname = "Vegetation C stocks", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$veg_gCm2)
    
    # Wood
    var_new  = ncvar_def("WOOD", unit="gC.m-2", longname = "Wood C stocks", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$wood_gCm2)
    
    # Close file
    nc_close(new_file)
    
  }
  if (i == "LPJ-GUESS") {
    ###
    ## LPJ-GUESS
    
    # model_counter = 7
    
    # LAI (m2m2)
    input = nc_open("./LPJ-GUESS/LPJ-GUESS_S3_lai.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "lai")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$lai_m2m2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # GPP (kgC/m2/s -> gC/m2/d)
    input = nc_open("./LPJ-GUESS/LPJ-GUESS_S3_gpp.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "gpp")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$gpp_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Ra (kgC/m2/s -> gC/m2/d)
    input = nc_open("./LPJ-GUESS/LPJ-GUESS_S3_ra.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "ra")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$ra_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    ## Rh (kgC/m2/s -> gC/m2/d)
    #input = nc_open("./LPJ-GUESS/LPJ-GUESS_S3_rh.nc")
    ##input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    #input_data = ncvar_get(input, "rh")
    #input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    ## Turn lat_in / long_in from vectors to arrays
    #input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    #input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    ## Check for lat / long in -90 / 90, -180 / 180 repectively
    #check_long = which(input_long > 180)
    #if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    ## Restrict time period, assume we want the last 22 years
    ## Assume that as all models should be ending in 2021 and monthly time step,
    ## we can extract the last 22 years
    #finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    #input_data = input_data[,,start:finish]
    ## Unit convertion
    #input_data = input_data * 1e3 * 86400
    ## Begin regridding
    #input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    ## Assign to output variable
    #trendy$rh_gCm2day[,,] = input_data$var
    ## Tidy
    #nc_close(input) ; rm(input_data) ; gc()
    
    # Fire (kgC/m2/s -> gC/m2/d)
    input = nc_open("./LPJ-GUESS/LPJ-GUESS_S3_fFire.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "fFire")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$fire_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # NBP (kgC/m2/s -> gC/m2/d)
    input = nc_open("./LPJ-GUESS/LPJ-GUESS_S3_nbp.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "nbp")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data_nbp = input_data * 1e3 * 86400 #save for NBE calculation
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$nbp_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # fLuc (kgC/m2/s -> gC/m2/d)
    input = nc_open("./LPJ-GUESS/LPJ-GUESS_S3_fLuc.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "fLuc")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180) 
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data_fluc = input_data * 1e3 * 86400
    # Tidy
    rm(input_data) ; gc()
    
    # NBE (kgC/m2/s -> gC/m2/d)
    # Unit convertion
    input_data = (input_data_nbp+input_data_fluc)*-1
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$nbe_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Veg (kgC/m2 -> gC/m2)
    input = nc_open("./LPJ-GUESS/LPJ-GUESS_S3_cVeg.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "cVeg")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$veg_gCm2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Soil (kgC/m2 -> gC/m2)
    input = nc_open("./LPJ-GUESS/LPJ-GUESS_S3_cSoil.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "cSoil")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$soil_gCm2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Wood (kgC/m2 -> gC/m2)
    input = nc_open("./LPJ-GUESS/LPJ-GUESS_S3_cWood.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "cWood")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180) 
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$wood_gCm2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    ## define dimension
    lat_dimen <- ncdim_def( "lat", units="degree north (-90->90)", trendy$latitude[1,] )
    long_dimen <- ncdim_def( "lon", units="degree east (-180->180)", trendy$longitude[,1] )
    time_dimen <- ncdim_def( "time", units="", trendy$year)
    # model_dimen <- ncdim_def( "model", units="", c(1:length(trendy$model_list)))
    dimnchar <- ncdim_def("nchar", "", 1:18, create_dimvar=FALSE ) # Maximum number of characters used in string vector
    
    ## define output variable
    # var0 = ncvar_def("model_list", units = "-", longname = "Model name", 
    #                  dim=list(dimnchar,model_dimen), prec="char", compression = 9)
    var0 = ncvar_def("area", units = "m2", longname = paste("Pixel area",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    var1 = ncvar_def("land_fraction", units = "1", longname = paste("Fraction of pixel which is land",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    var2 = ncvar_def("vegetated_cover", units = "1", longname = paste("Pixels determined to have vegetation in CARDAMOM grid",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    
    # Define the output file name
    output_name = paste(output_dir,output_prefix,i,"_trendy_output",output_suffix,".nc", sep="")
    # Delete if the file currently exists
    if (file.exists(output_name)) {file.remove(output_name)}
    # Create the empty file space
    new_file=nc_create(filename=output_name, vars=list(var0,var1,var2), force_v4 = TRUE)
    # Load first variable into the file
    # # Model list
    # ncvar_put(new_file, var0, trendy$model_list)
    # Grid area 
    ncvar_put(new_file, var0, trendy$area_m2)
    # Land fraction
    ncvar_put(new_file, var1, trendy$land_fraction)
    # Land fraction
    ncvar_put(new_file, var2, trendy$vegetated_cover)
    # Close the existing file to ensure its written to file
    nc_close(new_file)
    
    ## Re-open the file so that we can add to it a variable at a time
    new_file <- nc_open( output_name, write=TRUE )
    
    # Leaf area index
    var_new  = ncvar_def("LAI", unit="m2.m-2", longname = "Leaf Area Index", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$lai_m2m2)
    
    # Gross Primary Productivity
    var_new  = ncvar_def("GPP", unit="gC.m-2.d-1", longname = "Gross Primary Productivity", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$gpp_gCm2day)
    
    # Autotrophic respiration
    var_new  = ncvar_def("Ra", unit="gC.m-2.d-1", longname = "Autotrophic respiration", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$ra_gCm2day)
    
    # Heterotrophic respiration
    var_new  = ncvar_def("Rh", unit="gC.m-2.d-1", longname = "Heterotrophic respiration", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$rh_gCm2day)
    
    # Fire
    var_new  = ncvar_def("FIRE", unit="gC.m-2.d-1", longname = "Fire C emissions", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$fire_gCm2day)
    
    # NBP
    var_new  = ncvar_def("NBP", unit="gC.m-2.d-1", longname = "Net Biome Productivity", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$nbp_gCm2day)
    
    # NBE
    var_new  = ncvar_def("NBE", unit="gC.m-2.d-1", longname = "Net Biome Exchange", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$nbe_gCm2day)
    
    # Soil
    var_new  = ncvar_def("SOIL", unit="gC.m-2", longname = "Soil C stocks", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$soil_gCm2)
    
    # Vegetation
    var_new  = ncvar_def("TOT", unit="gC.m-2", longname = "Vegetation C stocks", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$veg_gCm2)
    
    # Wood
    var_new  = ncvar_def("WOOD", unit="gC.m-2", longname = "Wood C stocks", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$wood_gCm2)
    
    # Close file
    nc_close(new_file)
    
  }
  if (i == "LPX-Bern") {
    ###
    ## LPX-Bern
    
    # model_counter = 8
    
    # LAI (m2m2)
    input = nc_open("./LPX-Bern/LPX-Bern_S3_lai.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "lai")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$lai_m2m2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # GPP (kgC/m2/s -> gC/m2/d)
    input = nc_open("./LPX-Bern/LPX-Bern_S3_gpp.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "gpp")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$gpp_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Ra (kgC/m2/s -> gC/m2/d)
    input = nc_open("./LPX-Bern/LPX-Bern_S3_ra.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "ra")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$ra_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Rh (kgC/m2/s -> gC/m2/d)
    input = nc_open("./LPX-Bern/LPX-Bern_S3_rh.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "rh")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$rh_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Fire (kgC/m2/s -> gC/m2/d)
    input = nc_open("./LPX-Bern/LPX-Bern_S3_fFire.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "fFire")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$fire_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # NBP (kgC/m2/s -> gC/m2/d)
    input = nc_open("./LPX-Bern/LPX-Bern_S3_nbp.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "nbp")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data_nbp = input_data * 1e3 * 86400 #save for NBE calculation
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$nbp_gCm2day[,,] = input_data$var
    
    # NBE (kgC/m2/s -> gC/m2/d)
    # Unit convertion
    input_data_nbp = input_data_nbp * -1
    # Begin regridding
    input_data_nbp = regrid_func(input_data_nbp,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$nbe_gCm2day[,,] = input_data_nbp$var
    # Tidy
    nc_close(input) ; rm(input_data) ; rm(input_data_nbp) ; gc()
    
    # Veg (kgC/m2 -> gC/m2)
    input = nc_open("./LPX-Bern/LPX-Bern_S3_cVeg.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "cVeg")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$veg_gCm2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Soil (kgC/m2 -> gC/m2)
    input = nc_open("./LPX-Bern/LPX-Bern_S3_cSoil.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "cSoil")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$soil_gCm2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Wood (kgC/m2 -> gC/m2)
    input = nc_open("./LPX-Bern/LPX-Bern_S3_cWood.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "cWood")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180) 
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$wood_gCm2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    ## define dimension
    lat_dimen <- ncdim_def( "lat", units="degree north (-90->90)", trendy$latitude[1,] )
    long_dimen <- ncdim_def( "lon", units="degree east (-180->180)", trendy$longitude[,1] )
    time_dimen <- ncdim_def( "time", units="", trendy$year)
    # model_dimen <- ncdim_def( "model", units="", c(1:length(trendy$model_list)))
    dimnchar <- ncdim_def("nchar", "", 1:18, create_dimvar=FALSE ) # Maximum number of characters used in string vector
    
    ## define output variable
    # var0 = ncvar_def("model_list", units = "-", longname = "Model name", 
    #                  dim=list(dimnchar,model_dimen), prec="char", compression = 9)
    var0 = ncvar_def("area", units = "m2", longname = paste("Pixel area",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    var1 = ncvar_def("land_fraction", units = "1", longname = paste("Fraction of pixel which is land",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    var2 = ncvar_def("vegetated_cover", units = "1", longname = paste("Pixels determined to have vegetation in CARDAMOM grid",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    
    # Define the output file name
    output_name = paste(output_dir,output_prefix,i,"_trendy_output",output_suffix,".nc", sep="")
    # Delete if the file currently exists
    if (file.exists(output_name)) {file.remove(output_name)}
    # Create the empty file space
    new_file=nc_create(filename=output_name, vars=list(var0,var1,var2), force_v4 = TRUE)
    # Load first variable into the file
    # # Model list
    # ncvar_put(new_file, var0, trendy$model_list)
    # Grid area 
    ncvar_put(new_file, var0, trendy$area_m2)
    # Land fraction
    ncvar_put(new_file, var1, trendy$land_fraction)
    # Land fraction
    ncvar_put(new_file, var2, trendy$vegetated_cover)
    # Close the existing file to ensure its written to file
    nc_close(new_file)
    
    ## Re-open the file so that we can add to it a variable at a time
    new_file <- nc_open( output_name, write=TRUE )
    
    # Leaf area index
    var_new  = ncvar_def("LAI", unit="m2.m-2", longname = "Leaf Area Index", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$lai_m2m2)
    
    # Gross Primary Productivity
    var_new  = ncvar_def("GPP", unit="gC.m-2.d-1", longname = "Gross Primary Productivity", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$gpp_gCm2day)
    
    # Autotrophic respiration
    var_new  = ncvar_def("Ra", unit="gC.m-2.d-1", longname = "Autotrophic respiration", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$ra_gCm2day)
    
    # Heterotrophic respiration
    var_new  = ncvar_def("Rh", unit="gC.m-2.d-1", longname = "Heterotrophic respiration", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$rh_gCm2day)
    
    # Fire
    var_new  = ncvar_def("FIRE", unit="gC.m-2.d-1", longname = "Fire C emissions", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$fire_gCm2day)
    
    # NBP
    var_new  = ncvar_def("NBP", unit="gC.m-2.d-1", longname = "Net Biome Productivity", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$nbp_gCm2day)
    
    # NBE
    var_new  = ncvar_def("NBE", unit="gC.m-2.d-1", longname = "Net Biome Exchange", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$nbe_gCm2day)
    
    # Soil
    var_new  = ncvar_def("SOIL", unit="gC.m-2", longname = "Soil C stocks", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$soil_gCm2)
    
    # Vegetation
    var_new  = ncvar_def("TOT", unit="gC.m-2", longname = "Vegetation C stocks", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$veg_gCm2)
    
    # Wood
    var_new  = ncvar_def("WOOD", unit="gC.m-2", longname = "Wood C stocks", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$wood_gCm2)
    
    # Close file
    nc_close(new_file)
    
  }
  if (i == "OCRHIDEE") {
    ###
    ## ORCHIDEE
    
    # model_counter = 9
    
    # LAI (m2m2)
    input = nc_open("./ORCHIDEE/ORCHIDEE_S3_lai.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "lai")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$lai_m2m2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # GPP (kgC/m2/s -> gC/m2/d)
    input = nc_open("./ORCHIDEE/ORCHIDEE_S3_gpp.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "gpp")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$gpp_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Ra (kgC/m2/s -> gC/m2/d)
    input = nc_open("./ORCHIDEE/ORCHIDEE_S3_ra.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "ra")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$ra_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Rh (kgC/m2/s -> gC/m2/d)
    input = nc_open("./ORCHIDEE/ORCHIDEE_S3_rh.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "rh")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$rh_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Fire (kgC/m2/s -> gC/m2/d)
    input = nc_open("./ORCHIDEE/ORCHIDEE_S3_fFire.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "fFire")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$fire_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # NBP (kgC/m2/s -> gC/m2/d)
    input = nc_open("./ORCHIDEE/ORCHIDEE_S3_nbp.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "nbp")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data_nbp = input_data * 1e3 * 86400 #save for NBE calculation
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$nbp_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # fLuc (kgC/m2/s -> gC/m2/d)
    input = nc_open("./OCRHIDEE/OCRHIDEE_S3_fLuc.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "fLuc")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180) 
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data_fluc = input_data * 1e3 * 86400
    # Tidy
    rm(input_data) ; gc()
    
    # NBE (kgC/m2/s -> gC/m2/d)
    # Unit convertion
    input_data = (input_data_nbp+input_data_fluc)*-1
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$nbe_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Veg (kgC/m2 -> gC/m2)
    input = nc_open("./ORCHIDEE/ORCHIDEE_S3_cVeg.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "cVeg")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$veg_gCm2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Soil (kgC/m2 -> gC/m2)
    input = nc_open("./ORCHIDEE/ORCHIDEE_S3_cSoil.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "cSoil")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$soil_gCm2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Wood (kgC/m2 -> gC/m2)
    input = nc_open("./ORCHIDEE/ORCHIDEE_S3_cWood.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "cWood")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180) 
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$wood_gCm2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    ## define dimension
    lat_dimen <- ncdim_def( "lat", units="degree north (-90->90)", trendy$latitude[1,] )
    long_dimen <- ncdim_def( "lon", units="degree east (-180->180)", trendy$longitude[,1] )
    time_dimen <- ncdim_def( "time", units="", trendy$year)
    # model_dimen <- ncdim_def( "model", units="", c(1:length(trendy$model_list)))
    dimnchar <- ncdim_def("nchar", "", 1:18, create_dimvar=FALSE ) # Maximum number of characters used in string vector
    
    ## define output variable
    # var0 = ncvar_def("model_list", units = "-", longname = "Model name", 
    #                  dim=list(dimnchar,model_dimen), prec="char", compression = 9)
    var0 = ncvar_def("area", units = "m2", longname = paste("Pixel area",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    var1 = ncvar_def("land_fraction", units = "1", longname = paste("Fraction of pixel which is land",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    var2 = ncvar_def("vegetated_cover", units = "1", longname = paste("Pixels determined to have vegetation in CARDAMOM grid",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    
    # Define the output file name
    output_name = paste(output_dir,output_prefix,i,"_trendy_output",output_suffix,".nc", sep="")
    # Delete if the file currently exists
    if (file.exists(output_name)) {file.remove(output_name)}
    # Create the empty file space
    new_file=nc_create(filename=output_name, vars=list(var0,var1,var2), force_v4 = TRUE)
    # Load first variable into the file
    # # Model list
    # ncvar_put(new_file, var0, trendy$model_list)
    # Grid area 
    ncvar_put(new_file, var0, trendy$area_m2)
    # Land fraction
    ncvar_put(new_file, var1, trendy$land_fraction)
    # Land fraction
    ncvar_put(new_file, var2, trendy$vegetated_cover)
    # Close the existing file to ensure its written to file
    nc_close(new_file)
    
    ## Re-open the file so that we can add to it a variable at a time
    new_file <- nc_open( output_name, write=TRUE )
    
    # Leaf area index
    var_new  = ncvar_def("LAI", unit="m2.m-2", longname = "Leaf Area Index", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$lai_m2m2)
    
    # Gross Primary Productivity
    var_new  = ncvar_def("GPP", unit="gC.m-2.d-1", longname = "Gross Primary Productivity", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$gpp_gCm2day)
    
    # Autotrophic respiration
    var_new  = ncvar_def("Ra", unit="gC.m-2.d-1", longname = "Autotrophic respiration", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$ra_gCm2day)
    
    # Heterotrophic respiration
    var_new  = ncvar_def("Rh", unit="gC.m-2.d-1", longname = "Heterotrophic respiration", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$rh_gCm2day)
    
    # Fire
    var_new  = ncvar_def("FIRE", unit="gC.m-2.d-1", longname = "Fire C emissions", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$fire_gCm2day)
    
    # NBP
    var_new  = ncvar_def("NBP", unit="gC.m-2.d-1", longname = "Net Biome Productivity", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$nbp_gCm2day)
    
    # NBE
    var_new  = ncvar_def("NBE", unit="gC.m-2.d-1", longname = "Net Biome Exchange", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$nbe_gCm2day)
    
    # Soil
    var_new  = ncvar_def("SOIL", unit="gC.m-2", longname = "Soil C stocks", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$soil_gCm2)
    
    # Vegetation
    var_new  = ncvar_def("TOT", unit="gC.m-2", longname = "Vegetation C stocks", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$veg_gCm2)
    
    # Wood
    var_new  = ncvar_def("WOOD", unit="gC.m-2", longname = "Wood C stocks", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$wood_gCm2)
    
    # Close file
    nc_close(new_file)
    
  }
  if (i == "SDGVM") {
    ###
    ## SDGVM
    
    # model_counter = 10
    
    # LAI (m2m2)
    input = nc_open("./SDGVM/SDGVM_S3_lai.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "lai")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$lai_m2m2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # GPP (kgC/m2/s -> gC/m2/d)
    input = nc_open("./SDGVM/SDGVM_S3_gpp.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "gpp")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$gpp_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Ra (kgC/m2/s -> gC/m2/d)
    input = nc_open("./SDGVM/SDGVM_S3_ra.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "ra")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$ra_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Rh (kgC/m2/s -> gC/m2/d)
    input = nc_open("./SDGVM/SDGVM_S3_rh.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "rh")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$rh_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Fire (kgC/m2/s -> gC/m2/d)
    input = nc_open("./SDGVM/SDGVM_S3_fFire.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "fFire")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$fire_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # NBP (kgC/m2/s -> gC/m2/d)
    input = nc_open("./SDGVM/SDGVM_S3_nbp.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "nbp")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data_nbp = input_data * 1e3 * 86400 #save for NBE calculation
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$nbp_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # fLuc (kgC/m2/s -> gC/m2/d)
    input = nc_open("./SDGVM/SDGVM_S3_fLuc.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "fLuc")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180) 
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data_fluc = input_data * 1e3 * 86400
    # Tidy
    rm(input_data) ; gc()
    
    # NBE (kgC/m2/s -> gC/m2/d)
    # Unit convertion
    input_data = (input_data_nbp+input_data_fluc)*-1
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$nbe_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Veg (kgC/m2 -> gC/m2)
    input = nc_open("./SDGVM/SDGVM_S3_cVeg.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "cVeg")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$veg_gCm2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Soil (kgC/m2 -> gC/m2)
    input = nc_open("./SDGVM/SDGVM_S3_cSoil.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "cSoil")
    input_lat = ncvar_get(input,"latitude") ; input_long = ncvar_get(input, "longitude")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$soil_gCm2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    ## define dimension
    lat_dimen <- ncdim_def( "lat", units="degree north (-90->90)", trendy$latitude[1,] )
    long_dimen <- ncdim_def( "lon", units="degree east (-180->180)", trendy$longitude[,1] )
    time_dimen <- ncdim_def( "time", units="", trendy$year)
    # model_dimen <- ncdim_def( "model", units="", c(1:length(trendy$model_list)))
    dimnchar <- ncdim_def("nchar", "", 1:18, create_dimvar=FALSE ) # Maximum number of characters used in string vector
    
    ## define output variable
    # var0 = ncvar_def("model_list", units = "-", longname = "Model name", 
    #                  dim=list(dimnchar,model_dimen), prec="char", compression = 9)
    var0 = ncvar_def("area", units = "m2", longname = paste("Pixel area",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    var1 = ncvar_def("land_fraction", units = "1", longname = paste("Fraction of pixel which is land",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    var2 = ncvar_def("vegetated_cover", units = "1", longname = paste("Pixels determined to have vegetation in CARDAMOM grid",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    
    # Define the output file name
    output_name = paste(output_dir,output_prefix,i,"_trendy_output",output_suffix,".nc", sep="")
    # Delete if the file currently exists
    if (file.exists(output_name)) {file.remove(output_name)}
    # Create the empty file space
    new_file=nc_create(filename=output_name, vars=list(var0,var1,var2), force_v4 = TRUE)
    # Load first variable into the file
    # # Model list
    # ncvar_put(new_file, var0, trendy$model_list)
    # Grid area 
    ncvar_put(new_file, var0, trendy$area_m2)
    # Land fraction
    ncvar_put(new_file, var1, trendy$land_fraction)
    # Land fraction
    ncvar_put(new_file, var2, trendy$vegetated_cover)
    # Close the existing file to ensure its written to file
    nc_close(new_file)
    
    ## Re-open the file so that we can add to it a variable at a time
    new_file <- nc_open( output_name, write=TRUE )
    
    # Leaf area index
    var_new  = ncvar_def("LAI", unit="m2.m-2", longname = "Leaf Area Index", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$lai_m2m2)
    
    # Gross Primary Productivity
    var_new  = ncvar_def("GPP", unit="gC.m-2.d-1", longname = "Gross Primary Productivity", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$gpp_gCm2day)
    
    # Autotrophic respiration
    var_new  = ncvar_def("Ra", unit="gC.m-2.d-1", longname = "Autotrophic respiration", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$ra_gCm2day)
    
    # Heterotrophic respiration
    var_new  = ncvar_def("Rh", unit="gC.m-2.d-1", longname = "Heterotrophic respiration", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$rh_gCm2day)
    
    # Fire
    var_new  = ncvar_def("FIRE", unit="gC.m-2.d-1", longname = "Fire C emissions", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$fire_gCm2day)
    
    # NBP
    var_new  = ncvar_def("NBP", unit="gC.m-2.d-1", longname = "Net Biome Productivity", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$nbp_gCm2day)
    
    # NBE
    var_new  = ncvar_def("NBE", unit="gC.m-2.d-1", longname = "Net Biome Exchange", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$nbe_gCm2day)
    
    # Soil
    var_new  = ncvar_def("SOIL", unit="gC.m-2", longname = "Soil C stocks", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$soil_gCm2)
    
    # Vegetation
    var_new  = ncvar_def("TOT", unit="gC.m-2", longname = "Vegetation C stocks", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$veg_gCm2)
    
    # Close file
    nc_close(new_file)
    
  }
  if (i == "VISIT") {
    ###
    ## VISIT
    
    # model_counter = 11
    
    # LAI (m2m2)
    input = nc_open("./VISIT/VISIT_S3_lai.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "lai")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$lai_m2m2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # GPP (kgC/m2/s -> gC/m2/d)
    input = nc_open("./VISIT/VISIT_S3_gpp.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "gpp")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$gpp_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Ra (kgC/m2/s -> gC/m2/d)
    input = nc_open("./VISIT/VISIT_S3_ra.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "ra")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$ra_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Rh (kgC/m2/s -> gC/m2/d)
    input = nc_open("./VISIT/VISIT_S3_rh.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "rh")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$rh_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Fire (kgC/m2/s -> gC/m2/d)
    input = nc_open("./VISIT/VISIT_S3_fFire.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "fFire")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$fire_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # NBP (kgC/m2/s -> gC/m2/d)
    input = nc_open("./VISIT/VISIT_S3_nbp.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "nbp")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data_nbp = input_data * 1e3 * 86400 #save for NBE calculation
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$nbp_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # fLuc (kgC/m2/s -> gC/m2/d)
    input = nc_open("./VISIT/VISIT_S3_fLuc.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "fLuc")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180) 
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data_fluc = input_data * 1e3 * 86400
    # Tidy
    rm(input_data) ; gc()
    
    # NBE (kgC/m2/s -> gC/m2/d)
    # Unit convertion
    input_data = (input_data_nbp+input_data_fluc)*-1
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$nbe_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Veg (kgC/m2 -> gC/m2)
    input = nc_open("./VISIT/VISIT_S3_cVeg.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "cVeg")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$veg_gCm2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Soil (kgC/m2 -> gC/m2)
    input = nc_open("./VISIT/VISIT_S3_cSoil.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "cSoil")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$soil_gCm2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    ## define dimension
    lat_dimen <- ncdim_def( "lat", units="degree north (-90->90)", trendy$latitude[1,] )
    long_dimen <- ncdim_def( "lon", units="degree east (-180->180)", trendy$longitude[,1] )
    time_dimen <- ncdim_def( "time", units="", trendy$year)
    # model_dimen <- ncdim_def( "model", units="", c(1:length(trendy$model_list)))
    dimnchar <- ncdim_def("nchar", "", 1:18, create_dimvar=FALSE ) # Maximum number of characters used in string vector
    
    ## define output variable
    # var0 = ncvar_def("model_list", units = "-", longname = "Model name", 
    #                  dim=list(dimnchar,model_dimen), prec="char", compression = 9)
    var0 = ncvar_def("area", units = "m2", longname = paste("Pixel area",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    var1 = ncvar_def("land_fraction", units = "1", longname = paste("Fraction of pixel which is land",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    var2 = ncvar_def("vegetated_cover", units = "1", longname = paste("Pixels determined to have vegetation in CARDAMOM grid",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    
    # Define the output file name
    output_name = paste(output_dir,output_prefix,i,"_trendy_output",output_suffix,".nc", sep="")
    # Delete if the file currently exists
    if (file.exists(output_name)) {file.remove(output_name)}
    # Create the empty file space
    new_file=nc_create(filename=output_name, vars=list(var0,var1,var2), force_v4 = TRUE)
    # Load first variable into the file
    # # Model list
    # ncvar_put(new_file, var0, trendy$model_list)
    # Grid area 
    ncvar_put(new_file, var0, trendy$area_m2)
    # Land fraction
    ncvar_put(new_file, var1, trendy$land_fraction)
    # Land fraction
    ncvar_put(new_file, var2, trendy$vegetated_cover)
    # Close the existing file to ensure its written to file
    nc_close(new_file)
    
    ## Re-open the file so that we can add to it a variable at a time
    new_file <- nc_open( output_name, write=TRUE )
    
    # Leaf area index
    var_new  = ncvar_def("LAI", unit="m2.m-2", longname = "Leaf Area Index", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$lai_m2m2)
    
    # Gross Primary Productivity
    var_new  = ncvar_def("GPP", unit="gC.m-2.d-1", longname = "Gross Primary Productivity", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$gpp_gCm2day)
    
    # Autotrophic respiration
    var_new  = ncvar_def("Ra", unit="gC.m-2.d-1", longname = "Autotrophic respiration", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$ra_gCm2day)
    
    # Heterotrophic respiration
    var_new  = ncvar_def("Rh", unit="gC.m-2.d-1", longname = "Heterotrophic respiration", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$rh_gCm2day)
    
    # Fire
    var_new  = ncvar_def("FIRE", unit="gC.m-2.d-1", longname = "Fire C emissions", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$fire_gCm2day)
    
    # NBP
    var_new  = ncvar_def("NBP", unit="gC.m-2.d-1", longname = "Net Biome Productivity", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$nbp_gCm2day)
    
    # NBE
    var_new  = ncvar_def("NBE", unit="gC.m-2.d-1", longname = "Net Biome Exchange", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$nbe_gCm2day)
    
    # Soil
    var_new  = ncvar_def("SOIL", unit="gC.m-2", longname = "Soil C stocks", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$soil_gCm2)
    
    # Vegetation
    var_new  = ncvar_def("TOT", unit="gC.m-2", longname = "Vegetation C stocks", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$veg_gCm2)
    
    # Close file
    nc_close(new_file)
    
  }
  if (i == "VISIT-NIES") {
    ###
    ## VISIT-NIES
    
    # model_counter = 12
    
    # LAI (m2m2)
    input = nc_open("./VISIT-NIES/VISIT-NIES_S3_lai.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "lai")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$lai_m2m2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # GPP (kgC/m2/s -> gC/m2/d)
    input = nc_open("./VISIT-NIES/VISIT-NIES_S3_gpp.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "gpp")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$gpp_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Ra (kgC/m2/s -> gC/m2/d)
    input = nc_open("./VISIT-NIES/VISIT-NIES_S3_ra.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "ra")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$ra_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Rh (kgC/m2/s -> gC/m2/d)
    input = nc_open("./VISIT-NIES/VISIT-NIES_S3_rh.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "rh")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$rh_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Fire (kgC/m2/s -> gC/m2/d)
    input = nc_open("./VISIT-NIES/VISIT-NIES_S3_fFire.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "fFire")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$fire_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # NBP (kgC/m2/s -> gC/m2/d)
    input = nc_open("./VISIT-NIES/VISIT-NIES_S3_nbp.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "nbp")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data_nbp = input_data * 1e3 * 86400 #save for NBE calculation
    input_data = input_data * 1e3 * 86400
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$nbp_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # fLuc (kgC/m2/s -> gC/m2/d)
    input = nc_open("./VISIT-NIES/VISIT-NIES_S3_fLuc.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "fLuc")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180) 
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data_fluc = input_data * 1e3 * 86400
    # Tidy
    rm(input_data) ; gc()
    
    # NBE (kgC/m2/s -> gC/m2/d)
    # Unit convertion
    input_data = (input_data_nbp+input_data_fluc)*-1
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$nbe_gCm2day[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Veg (kgC/m2 -> gC/m2)
    input = nc_open("./VISIT-NIES/VISIT-NIES_S3_cVeg.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "cVeg")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$veg_gCm2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Soil (kgC/m2 -> gC/m2)
    input = nc_open("./VISIT-NIES/VISIT-NIES_S3_cSoil.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "cSoil")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180)
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$soil_gCm2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    # Wood (kgC/m2 -> gC/m2)
    input = nc_open("./VISIT-NIES/VISIT-NIES_S3_cWood.nc")
    #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
    input_data = ncvar_get(input, "cWood")
    input_lat = ncvar_get(input,"lat") ; input_long = ncvar_get(input, "lon")
    # Turn lat_in / long_in from vectors to arrays
    input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
    input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
    # Check for lat / long in -90 / 90, -180 / 180 repectively
    check_long = which(input_long > 180) 
    if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
    # Restrict time period, assume we want the last 22 years
    # Assume that as all models should be ending in 2021 and monthly time step,
    # we can extract the last 22 years
    finish = dim(input_data)[3] ; start = (finish - (22*12))+1
    input_data = input_data[,,start:finish]
    # Unit convertion
    input_data = input_data * 1e3 
    # Begin regridding
    input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
    # Assign to output variable
    trendy$wood_gCm2[,,] = input_data$var
    # Tidy
    nc_close(input) ; rm(input_data) ; gc()
    
    ## define dimension
    lat_dimen <- ncdim_def( "lat", units="degree north (-90->90)", trendy$latitude[1,] )
    long_dimen <- ncdim_def( "lon", units="degree east (-180->180)", trendy$longitude[,1] )
    time_dimen <- ncdim_def( "time", units="", trendy$year)
    # model_dimen <- ncdim_def( "model", units="", c(1:length(trendy$model_list)))
    dimnchar <- ncdim_def("nchar", "", 1:18, create_dimvar=FALSE ) # Maximum number of characters used in string vector
    
    ## define output variable
    # var0 = ncvar_def("model_list", units = "-", longname = "Model name", 
    #                  dim=list(dimnchar,model_dimen), prec="char", compression = 9)
    var0 = ncvar_def("area", units = "m2", longname = paste("Pixel area",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    var1 = ncvar_def("land_fraction", units = "1", longname = paste("Fraction of pixel which is land",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    var2 = ncvar_def("vegetated_cover", units = "1", longname = paste("Pixels determined to have vegetation in CARDAMOM grid",sep=""), 
                     dim=list(long_dimen,lat_dimen), missval = NA, prec="double", compression = 9)
    
    # Define the output file name
    output_name = paste(output_dir,output_prefix,i,"_trendy_output",output_suffix,".nc", sep="")
    # Delete if the file currently exists
    if (file.exists(output_name)) {file.remove(output_name)}
    # Create the empty file space
    new_file=nc_create(filename=output_name, vars=list(var0,var1,var2), force_v4 = TRUE)
    # Load first variable into the file
    # # Model list
    # ncvar_put(new_file, var0, trendy$model_list)
    # Grid area 
    ncvar_put(new_file, var0, trendy$area_m2)
    # Land fraction
    ncvar_put(new_file, var1, trendy$land_fraction)
    # Land fraction
    ncvar_put(new_file, var2, trendy$vegetated_cover)
    # Close the existing file to ensure its written to file
    nc_close(new_file)
    
    ## Re-open the file so that we can add to it a variable at a time
    new_file <- nc_open( output_name, write=TRUE )
    
    # Leaf area index
    var_new  = ncvar_def("LAI", unit="m2.m-2", longname = "Leaf Area Index", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$lai_m2m2)
    
    # Gross Primary Productivity
    var_new  = ncvar_def("GPP", unit="gC.m-2.d-1", longname = "Gross Primary Productivity", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$gpp_gCm2day)
    
    # Autotrophic respiration
    var_new  = ncvar_def("Ra", unit="gC.m-2.d-1", longname = "Autotrophic respiration", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$ra_gCm2day)
    
    # Heterotrophic respiration
    var_new  = ncvar_def("Rh", unit="gC.m-2.d-1", longname = "Heterotrophic respiration", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$rh_gCm2day)
    
    # Fire
    var_new  = ncvar_def("FIRE", unit="gC.m-2.d-1", longname = "Fire C emissions", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$fire_gCm2day)
    
    # NBP
    var_new  = ncvar_def("NBP", unit="gC.m-2.d-1", longname = "Net Biome Productivity", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$nbp_gCm2day)
    
    # NBE
    var_new  = ncvar_def("NBE", unit="gC.m-2.d-1", longname = "Net Biome Exchange", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$nbe_gCm2day)
    
    # Soil
    var_new  = ncvar_def("SOIL", unit="gC.m-2", longname = "Soil C stocks", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$soil_gCm2)
    
    # Vegetation
    var_new  = ncvar_def("TOT", unit="gC.m-2", longname = "Vegetation C stocks", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$veg_gCm2)
    
    # Wood
    var_new  = ncvar_def("WOOD", unit="gC.m-2", longname = "Wood C stocks", dim=list(long_dimen,lat_dimen,time_dimen), missval = NA, prec="double",compression = 9)
    new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
    ncvar_put(new_file, var_new,  trendy$wood_gCm2)
    
    # Close file
    nc_close(new_file)
  }
}











