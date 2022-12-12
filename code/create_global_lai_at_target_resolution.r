
###
## Script to extract Copernicus 300 m LAI information, restrict to a target area and then aggragte to 4x5 degree
###

# Specify the resolution in degrees (WGS-84) you want to final output at
res_to_aggregate = c(0.125,0.125)

library(compiler) ; library(fields) ; library(ncdf4) ; library(raster) 
#source("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/R_functions/function_closest2d.r")
# Make a list of all the files we need to process
list_of_files = list.files("/disk/scratch/local.2/copernicus/LAI_1km_linked", full.names=TRUE, recursive = TRUE)
list_of_files = list_of_files[grepl("c_gls_LAI",list_of_files)]
list_of_files = list_of_files[grepl(".nc",list_of_files)]
#list_of_files = list_of_files[(length(list_of_files)*0.8):length(list_of_files)]
# Determine what our output file name will be
outfilename = unlist(list_of_files)
outfilename = gsub("c_gls_LAI-RT0","c_gls_LAI",outfilename) # remove this flag as it makes life easier to read into CARDAMOM later
outfilename = gsub("c_gls_LAI-RT1","c_gls_LAI",outfilename)
outfilename = gsub("c_gls_LAI-RT2","c_gls_LAI",outfilename)
outfilename = gsub("c_gls_LAI-RT3","c_gls_LAI",outfilename)
outfilename = gsub("c_gls_LAI-RT4","c_gls_LAI",outfilename)
outfilename = gsub("c_gls_LAI-RT5","c_gls_LAI",outfilename)
outfilename = gsub("c_gls_LAI-RT6","c_gls_LAI",outfilename)
outfilename = strsplit(outfilename, "/")
outfilename = unlist(outfilename) ; outfilename = outfilename[seq(7,length(outfilename),7)]
outfilename = paste("/disk/scratch/local.2/copernicus/LAI_0.125deg/",outfilename, sep="")

# Determine the i,j range which is covered but the desired area
# NOTE: I assume here that all files to be read have the same spatial grid underlying them
data1 = nc_open(list_of_files[1])
lat_orig = ncvar_get(data1, "lat") ; long_orig = ncvar_get(data1, "lon")
lat_orig = lat_orig[length(lat_orig):1]
nc_close(data1)
# determine the original resolution of the dataset
orig_resolution = abs(diff(lat_orig[1:2]))
radius = res_to_aggregate / orig_resolution

# Extract the current founds for the dataset
lat_bounds = range(lat_orig) 
long_bounds = range(long_orig)

# Create raster with user specified lat / long information
#dat1 = list(x = long_orig, y = lat_orig, z = lai)
#r=raster(dat1)

# Begin loop through each file
for (n in seq(1, length(list_of_files))) {

     # If file does not exist we will create it
     if (file.exists(outfilename[n]) == FALSE) {
         # open the next file
         data1 = nc_open(list_of_files[n])
         # read in the leaf area index...
         lai = ncvar_get(data1, "LAI")
         # convert to raster and aggregate
         lai = raster(t(lai), xmn=long_bounds[1], xmx=long_bounds[2], ymn=lat_bounds[1], ymx=lat_bounds[2])
         lai = aggregate(lai, radius, fun = mean, na.rm=TRUE) ; gc()
         # extract lat / long information
         lat_out = coordinates(lai)
         long_out = lat_out[,1] ; lat_out = lat_out[,2]
         # read in the rmse for the lai...
         lai_rmse = ncvar_get(data1, "RMSE")
         lai_rmse = raster(t(lai_rmse), xmn=long_bounds[1], xmx=long_bounds[2], ymn=lat_bounds[1], ymx=lat_bounds[2])
         lai_rmse = aggregate(lai_rmse, radius, fun = mean, na.rm=TRUE) ; gc()
         # tidy up 
         nc_close(data1)

         # re-structure back into array
         dims = dim(lai)
         lai = (array(as.vector(lai), dim=c(dims[2],dims[1])))         ; lai_rmse = (array(as.vector(lai_rmse), dim=c(dims[2],dims[1])))
         lat_out = (array(as.vector(lat_out), dim=c(dims[2],dims[1]))) ; long_out = (array(as.vector(long_out), dim=c(dims[2],dims[1])))
         # Create lat / long vector
         lat_out = lat_out[1,] ; long_out = long_out[,1]

         # write out the files to netcdf 
    
         # define dimensions
         lat_axis = ncdim_def("lat_axis", "", vals=c(1:dim(lai)[2]))
         long_axis = ncdim_def("long_axis", "", vals=c(1:dim(lai)[1]))

         # define variable
         var0 = ncvar_def(name="LAI",units="m2/m2",list(long_axis,lat_axis),-9999,longname="Leaf Area Index (m2/m2)", prec="double", compression = 9)
         var1 = ncvar_def(name="RMSE",units="m2/m2",list(long_axis,lat_axis),-9999,longname="Root Mean Square Error", prec="double", compression = 9)
         var2 = ncvar_def(name="lat",units="-",list(lat_axis),-9999,longname="Latitude (-90/90)", prec="double", compression = 9)
         var3 = ncvar_def(name="lon",units="-",list(long_axis),-9999,longname="Longitude (-180/180)", prec="double", compression = 9)
         # create a list of the variables to be added
         varlist = list(var0,var1,var2,var3)
         # create ncdf file with these variables and dimensions
         outfile = nc_create(filename=outfilename[n], varlist, verbose=FALSE, force_v4 = TRUE)
         # now add the variables to the file
         ncvar_put(outfile, var0, lai)
         ncvar_put(outfile, var1, lai_rmse)
         ncvar_put(outfile, var2, lat_out)
         ncvar_put(outfile, var3, long_out)
         # now close the file
         nc_close(outfile)

         rm(lai,lai_rmse) ; gc()
     } # does output file already exist

} # loop through each file



