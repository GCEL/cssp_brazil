load("M:/CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/nbe_data_alone/RESULTS_PROCESSED/nbe_data_alone_stock_flux.RData")

grid_output$mean_gpp_gCm2day
print(tbl_df(grid_output$mean_gpp_gCm2day[,,1]), n=40)

extract_check_grid <- function (model_variant,region,region_name,variable_name,reference,tp){
  for (i in model_variant) {
    if (i=='raw_GEOSCHEM'){
      data <- nbe_stack_15_19_monthly
    } else {
      data <- brick(paste(nbe_prefix,model_variant[1],nbe_suffix,sep=""),varname=variable_name)
    }
    a <- data$X2015.01.01
    cropped <- crop(a, extent(region))
    masked1 <- mask(cropped, region)
    masked2 <- masked1 > -Inf
    poly <- rasterToPolygons(masked2, dissolve=TRUE)
    data_region <- mask(data,poly)
    data_region_1519<-data_region[[which(getZ(data_region) >= as.Date("2015-01-01") & getZ(data_region) <= as.Date("2019-12-01"))]]
    df_data_region_1519 <- stackApply(data_region_1519, indices =  rep(1,nlayers(data_region_1519)), fun = "mean")
    df_data_region_1519 <- as.matrix(t(df_data_region_1519))
    
  }
  return(df_data_region_1519)
}


View(extract_check_grid(nbe_model_variants[2],nbe_pixel_data[[1]],nbe_pixel_names[1],nbe_var[1],reference_nbe_data,nbe_time_period))

# nw [10,58]; sw [13,45] ; ec [23,53] ; bs [32,45] ; gs [26,62]