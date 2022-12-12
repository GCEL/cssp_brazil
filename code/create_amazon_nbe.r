
###
## Script to extract Liang's GEOSCHEM net biome exchange of CO2 (NBE; gC/m2/day) information and aggregate to target resolution (gC/m2/day)
## This script includes a specific modification to extract multiple years from the same file but to output them in a more CARDAMOM fieldly format
###

# Specify the resolution in degrees (WGS-84) you want to final output at
target_resolution = c(1,1)
# Specify the input directory
indir = "/exports/csce/datastore/geos/groups/gcel/AtmosphericInversions/GEOSCHEM/oco2_XCO2_2x2.5_to_1x1/LiangRaw/"
# Specify the output directory
outdir = "/exports/csce/datastore/geos/groups/gcel/AtmosphericInversions/GEOSCHEM/global_1deg_v2/"
# Desired year and month range
years_to_do = c(2015:2019) ; months = rep(NA, 12) ; for (m in seq(1, length(months))){ months[m] = sprintf('%02i',m) } 
# Convert into fill timeseries of desired
datetime_wanted = paste(rep(years_to_do, each = length(months)),months,sep="")
# Load libraries needed
library(compiler) ; library(fields) ; library(ncdf4) ; library(raster) 
target = brick("/exports/csce/datastore/geos/groups/gcel/cssp_rainfor_amazon_brazil/rainfor_leeds_data/AbovegroundBiomass_Mg_perHa_111km.tif")

# List the available files
list_of_files = list.files(indir, full.names=TRUE, recursive = TRUE)
# Restrict to those with the correct file prefix
list_of_files = list_of_files[grepl("oco2_v10r_odiac_flux_",list_of_files)]
# Restrict to those with correct file ending
list_of_files = list_of_files[grepl("\\.nc$",list_of_files)]
# Determine what our output file name will be
#outfilename = gsub(indir,outdir,list_of_files)
outfilename = paste(outdir,"geoschem_nbe_",datetime_wanted,".nc", sep="")
# outfilename = paste(outdir,"geoschem_nbe_",datetime_wanted,".tiff", sep="")

nbe_list <- list()
nbe_unc_list <- list()

# Determine the i,j range which is covered but the desired area
# NOTE: I assume here that all files to be read have the same spatial grid underlying them
data1 = nc_open(list_of_files[1])
lat_orig = ncvar_get(data1, "lat") ; long_orig = ncvar_get(data1, "lon")
lat_orig = lat_orig[length(lat_orig):1]
nc_close(data1)
# determine the original resolution of the dataset
orig_resolution = abs(diff(lat_orig[1:2]))

# Extract the current founds for the dataset
lat_bounds = range(lat_orig) 
long_bounds = range(long_orig)

# Create raster with user specified lat / long information
#dat1 = list(x = long_orig, y = lat_orig, z = nbe)
#r=raster(dat1)

# Begin loop through each file
for (n in seq(1, length(list_of_files))) {

     # Determine which time steps we want are present in this file
     # open the next file
     data1 = nc_open(list_of_files[n])     
     # Extract timing information
     datetime = ncvar_get(data1, "start_date")
     # Convert to format compatible with the desired
     datetime = paste(datetime[1,],sprintf('%02i',datetime[2,]), sep="")
     # set day
     day_out = 15
     # Load NBE estimate (gC/m2/day)
     nbe_all = ncvar_get(data1, "flux")
     # Load NBE uncertainty estimate (gC/m2/day)
     nbe_unc_all = ncvar_get(data1, "error")

     # Filter all missing data flags -999 to NA
     filter = which(nbe_all == -999 | nbe_unc_all == -999)
     nbe_all[filter] = NA ; nbe_unc_all[filter] = NA     
     
     # tidy up 
     nc_close(data1)

     # Loop through our desired year and output from this file the overlaping bits
     for (s in seq(1, length(outfilename))) {
         
          # Is the current year found in this file?
          wanted_loc = which(datetime == datetime_wanted[s])
          
          # If file does not exist we will create it
          if (length(wanted_loc) > 0 & file.exists(outfilename[s]) == FALSE) {

              # Extract the correct timestep from the database
              nbe = nbe_all[,,wanted_loc]
              # convert to raster assuming crs is WGS 84, espg: 4326
              nbe = raster(t(nbe), xmn=long_bounds[1]-(0.5*orig_resolution), xmx=long_bounds[2]+(0.5*orig_resolution), 
                                   ymn=lat_bounds[1]-(0.5*orig_resolution), ymx=lat_bounds[2]+(0.5*orig_resolution), 
                                   crs = ("+init=epsg:4326")); nbe<-flip(nbe, 2)
              nbe = resample(nbe, target, method = "bilinear", na.rm=TRUE)

              # extract the correct timestep from the database
              nbe_unc = nbe_unc_all[,,wanted_loc]
              # convert to raster assuming crs is WGS 84, espg: 4326
              nbe_unc = raster(t(nbe_unc), xmn=long_bounds[1]-(0.5*orig_resolution), xmx=long_bounds[2]+(0.5*orig_resolution), 
                                           ymn=lat_bounds[1]-(0.5*orig_resolution), ymx=lat_bounds[2]+(0.5*orig_resolution), 
                                           crs = ("+init=epsg:4326")) ; nbe<-flip(nbe, 2)
              # Aggregate to target resolution
              nbe_unc = resample(nbe_unc, target, method = "bilinear", na.rm=TRUE)
              # add to list
              nbe_list[s]<- nbe
              nbe_unc_list[s]<- nbe_unc

              # create ncdf file with these variables and dimensions
              # outfile = nc_create(filename=outfilename[s], varlist, verbose=FALSE, force_v4 = TRUE)
              # writeRaster(nbe, outfilename[s], overwrite=TRUE,
              #             format="CDF",     varname="NBE", varunit="g.m-2.d-1",force_v4=TRUE,
              #             longname="Net Biome Exchange of CO2", xname="lon",   yname="lat", zname='time', zunit='days since 2015-1-15 00:00:00')
              # writeRaster(nbe, outfilename[s], overwrite=TRUE,format="GTiff")
              } # does output file already exist
          
     } # Loop through the desired output files
     return(nbe_list)
     return (nbe_unc_list)
} # loop through each input file



