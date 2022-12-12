

source("M:/CARDAMOM/CARDAMOM/R_functions/read_binary_file_format.r")

# Load libraries needed
library(compiler) ; library(fields) ; library(ncdf4) ; library(raster)  ; library(zoo)

# load the CARDAMOM files
load("M:/CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/no_woody_data_copernicus/infofile.RData")

# lai = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,228))
lai_unc = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,228))

# Loop through all sites
for (i in seq(1, length(PROJECT$sites))) {
  
  # Determine the correct site number for the current location
  n = PROJECT$sites[i]
  
  drivers = read_binary_file_format(paste("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/no_woody_data_copernicus/DATA/no_woody_data_copernicus_",n,".bin",sep=""))
  
  # determine the lat / long location within the grid
  slot_j=as.numeric(n)/PROJECT$long_dim
  slot_i=as.numeric(n)-(floor(slot_j)*PROJECT$long_dim)
  
  if(slot_i == 0) {slot_i = PROJECT$long_dim} ; slot_j=ceiling(slot_j)
  
  # save for later
  grid_output$i_location[n] = slot_i ; grid_output$j_location[n] = slot_j
  
 # loop through parameters + likelihood
  # lai[slot_i,slot_j,] = drivers$obs[,3]
  lai_unc[slot_i,slot_j,] = drivers$obs[,4]
  
}
lai_unc_b <- brick(lai_unc)
# lai_unc_b <- t(lai_unc_b)
# lai_unc_b <- flip(lai_unc_b,2)
# target = brick("R://ILAMB_beta_devel/RAINFOR_leeds_run/nbe_analysis/geoschem_CARDAMOM_Amazon_1deg_monthly_2001_updated_2019.nc",varname='NBE')
# lai_unc_b = resample(lai_unc_b, target, method = "bilinear", na.rm=TRUE)

lai_unc_b <- t(lai_unc_b);lai_unc_b <- flip(lai_unc_b,2);extent(lai_unc_b)<-extent(target);crs(lai_unc_b)<- "+init=epsg:4326"
