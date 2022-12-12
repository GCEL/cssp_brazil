
###
## Script to extract Copernicus 300 m LAI information, restrict to a target area and then aggragte to 4x5 degree
###

# Specify the resolution in degrees (WGS-84) you want to final output at
res_to_aggregate = c(1,1)

library(compiler) ; library(fields) ; library(ncdf4) ; library(raster)
target = brick("/exports/csce/datastore/geos/users/cnwobi/ILAMB_beta_devel/RAINFOR_leeds_run/nbe_analysis/geoschem_CARDAMOM_Amazon_1deg_monthly_2001_updated_2019.nc",varname='LAI')

#source("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/R_functions/function_closest2d.r")
# Make a list of all the files we need to process
list_of_files = list.files("/disk/scratch/local.2/copernicus/LAI_1km_linked", full.names=TRUE, recursive = TRUE)
list_of_files = list_of_files[grepl("c_gls_LAI_20",list_of_files)]
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

# determine the original resolution of the dataset
orig_resolution = abs(diff(lat_orig[1:2]))
radius = res_to_aggregate / orig_resolution

# Extract the current founds for the dataset
lat_bounds = range(lat_orig) 
long_bounds = range(long_orig)

lai_list <- list()
lai_unc_list <- list()

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
         lai = aggregate(lai, radius, fun = mean, na.rm=TRUE) 
         lai = resample(lai, target, method = "bilinear", na.rm=TRUE); gc()
         # extract lat / long information
         lat_out = coordinates(lai)
         long_out = lat_out[,1] ; lat_out = lat_out[,2]
         # read in the rmse for the lai...
         lai_rmse = ncvar_get(data1, "RMSE")
         lai_rmse = raster(t(lai_rmse), xmn=long_bounds[1], xmx=long_bounds[2], ymn=lat_bounds[1], ymx=lat_bounds[2])
         lai_rmse = aggregate(lai_rmse, radius, fun = mean, na.rm=TRUE) 
         lai_rmse = resample(lai_rmse, target, method = "bilinear", na.rm=TRUE); gc()
         
         lai_list[n]<- lai
         lai_unc_list[n]<- lai_rmse
         
         # tidy up 
         nc_close(data1)
         
         rm(lai,lai_rmse) ; gc()
     } # does output file already exist
  
} # loop through each file
# return(lai_list)
# return (lai_unc_list)



