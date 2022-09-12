#All codes being Used
#data prep pre-cursors----
compare_var_names <- c('WOOD','NPP_wood_flx','OUTPUT_wood_flx')
amazonia_subset <- shapefile("R:/cssp_brazil/cssp_brazil_R/data/amazonia_subset.shp")
mod_var <- c('norainfor','biomass','productivity','mortality')
time_period <- c('01/09','10/16')
prefix <- 'R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/DATA/benchmark/Amazon_subset_'
suffix <- '_1deg_monthly_2001_updated_2019.nc'
# paste(prefix,mod_var[1],suffix,sep="")

#functions for extraction----
extract_var_mean_biom <- function (region,data,reference){
  data <- crop(data, extent(reference))
  a <- data$X2001.01.01
  cropped <- crop(a, extent(region))
  masked1 <- mask(cropped, region)
  masked2 <- masked1 > -Inf
  poly <- rasterToPolygons(masked2, dissolve=TRUE)
  data_region <- mask(data,poly)
  data_2010<-data_region[[which(getZ(data_region) >= as.Date("2010-01-01") & getZ(data_region) <= as.Date("2010-12-01"))]]
  data_2010_mean <- stackApply(data_2010, indices =  rep(1,nlayers(data_2010)), fun = "mean")
  return(data_2010_mean)
}
extract_var_mean <- function (region,data,reference){
  data <- crop(data, extent(reference))
  a <- data$X2001.01.01
  cropped <- crop(a, extent(region))
  masked1 <- mask(cropped, region)
  masked2 <- masked1 > -Inf
  poly <- rasterToPolygons(masked2, dissolve=TRUE)
  data_region <- mask(data,poly)
  data_0110<-data_region[[which(getZ(data_region) >= as.Date("2001-01-01") & getZ(data_region) <= as.Date("2009-12-01"))]]
  data_1016<-data_region[[which(getZ(data_region) >= as.Date("2010-01-01") & getZ(data_region) <= as.Date("2016-12-01"))]]
  data_0110_mean <- stackApply(data_0110, indices =  rep(1,nlayers(data_0110)), fun = "mean")
  data_1016_mean <- stackApply(data_1016, indices =  rep(1,nlayers(data_1016)), fun = "mean")
  return(c(data_0110_mean,data_1016_mean))
}

extract_var_mean_biom_subset <- function (data){
  data_2010<-data[[which(getZ(data) >= as.Date("2010-01-01") & getZ(data) <= as.Date("2010-12-01"))]]
  data_2010_mean <- stackApply(data_2010, indices =  rep(1,nlayers(data_2010)), fun = "mean")
  return(data_2010_mean)
}
extract_var_mean_subset <- function (data){
  data_0109<-data[[which(getZ(data) >= as.Date("2001-01-01") & getZ(data) <= as.Date("2009-12-01"))]]
  data_1016<-data[[which(getZ(data) >= as.Date("2010-01-01") & getZ(data) <= as.Date("2016-12-01"))]]
  data_0109_mean <- stackApply(data_0109, indices =  rep(1,nlayers(data_0109)), fun = "mean")
  data_1016_mean <- stackApply(data_1016, indices =  rep(1,nlayers(data_1016)), fun = "mean")
  return(c(data_0109_mean,data_1016_mean))
}

#resample CARD_RAINFOR and RAINFOR data and plot 1:1
res_df_merge_plot_cwood <- function (benchmark,model,model_type){
  model_res <- resample(model,benchmark)
  model_res_df <- as.data.frame(model_res, xy=TRUE)
  benchmark_df <- as.data.frame(benchmark, xy=TRUE)
  model_benchmark_df <- merge(benchmark_df, model_res_df, by=c("x","y"))
  names(model_benchmark_df)<-c('x','y','rainfor','cardamom')
  model_benchmark_df <- model_benchmark_df[!is.na(model_benchmark_df$rainfor)&!is.na(model_benchmark_df$cardamom),]
  if (model_type=='rainfor') {
    plot(model_benchmark_df$cardamom, model_benchmark_df$rainfor, main="AGB comparison 2001-2010",
         xlab="CARDAMOM_RAINFOR Wood biomass g.m-2 ", ylab="RAINFOR Wood biomass g.m-2", pch=19, xlim=c(0,25000), ylim=c(0,25000))
    abline(coef = c(0,1),col='red', lwd=3)
  }
  else {
    plot(model_benchmark_df$cardamom, model_benchmark_df$rainfor, main="AGB comparison 2001-2010",
         xlab="CARDAMOM_NORAINFOR Wood biomass g.m-2 ", ylab="RAINFOR Wood biomass g.m-2", pch=19, xlim=c(0,25000), ylim=c(0,25000))
    abline(coef = c(0,1),col='red', lwd=3)
  }
}

res_df_merge_plot_nppwood <- function (benchmark,model,model_type){
  model_res <- resample(model,benchmark)
  model_res_df <- as.data.frame(model_res, xy=TRUE)
  benchmark_df <- as.data.frame(benchmark, xy=TRUE)
  model_benchmark_df <- merge(benchmark_df, model_res_df, by=c("x","y"))
  names(model_benchmark_df)<-c('x','y','rainfor','cardamom')
  model_benchmark_df <- model_benchmark_df[!is.na(model_benchmark_df$rainfor)&!is.na(model_benchmark_df$cardamom),]
  if (model_type=='rainfor') {
    plot(model_benchmark_df$cardamom, model_benchmark_df$rainfor, main="Woody NPP comparison 2000/1-2009",
         xlab="CARDAMOM_RAINFOR Woody NPP g.m-2.d-1 ", ylab="RAINFOR Wood productivity g.m-2.d-1", pch=19, xlim=c(0,4), ylim=c(0,4))
    abline(coef = c(0,1),col='red', lwd=3)
  }
  else {
    plot(model_benchmark_df$cardamom, model_benchmark_df$rainfor, main="Woody NPP comparison 2000/1-2009",
         xlab="CARDAMOM_NORAINFOR Woody NPP g.m-2.d-1 ", ylab="RAINFOR Wood productivity g.m-2.d-1", pch=19, xlim=c(0,4), ylim=c(0,4))
    abline(coef = c(0,1),col='red', lwd=3)
  }
}

res_df_merge_plot_nppwood_2 <- function (benchmark,model,model_type){
  model_res <- resample(model,benchmark)
  model_res_df <- as.data.frame(model_res, xy=TRUE)
  benchmark_df <- as.data.frame(benchmark, xy=TRUE)
  model_benchmark_df <- merge(benchmark_df, model_res_df, by=c("x","y"))
  names(model_benchmark_df)<-c('x','y','rainfor','cardamom')
  model_benchmark_df <- model_benchmark_df[!is.na(model_benchmark_df$rainfor)&!is.na(model_benchmark_df$cardamom),]
  if (model_type=='rainfor') {
    plot(model_benchmark_df$cardamom, model_benchmark_df$rainfor, main="Woody NPP comparison 2010-2016",
         xlab="CARDAMOM_RAINFOR Woody NPP g.m-2.d-1 ", ylab="RAINFOR Wood productivity g.m-2.d-1", pch=19, xlim=c(0,4), ylim=c(0,4))
    abline(coef = c(0,1),col='red', lwd=3)
  }
  else {
    plot(model_benchmark_df$cardamom, model_benchmark_df$rainfor, main="Woody NPP comparison 2010-2016",
         xlab="CARDAMOM_NORAINFOR Woody NPP g.m-2.d-1 ", ylab="RAINFOR Wood productivity g.m-2.d-1", pch=19, xlim=c(0,4), ylim=c(0,4))
    abline(coef = c(0,1),col='red', lwd=3)
  }
}

res_df_merge_plot_outputwood <- function (benchmark,model,model_type){
  model_res <- resample(model,benchmark)
  model_res_df <- as.data.frame(model_res, xy=TRUE)
  benchmark_df <- as.data.frame(benchmark, xy=TRUE)
  model_benchmark_df <- merge(benchmark_df, model_res_df, by=c("x","y"))
  names(model_benchmark_df)<-c('x','y','rainfor','cardamom')
  model_benchmark_df <- model_benchmark_df[!is.na(model_benchmark_df$rainfor)&!is.na(model_benchmark_df$cardamom),]
  if (model_type=='rainfor') {
    plot(model_benchmark_df$cardamom, model_benchmark_df$rainfor, main="Woody Carbon losses comparison 2000/1-2009",
         xlab="CARDAMOM_RAINFOR Output Wood g.m-2.d-1 ", ylab="RAINFOR Wood mortality g.m-2.d-1", pch=19, xlim=c(0,4), ylim=c(0,4))
    abline(coef = c(0,1),col='red', lwd=3)
  }
  else {
    plot(model_benchmark_df$cardamom, model_benchmark_df$rainfor, main="Woody Carbon losses comparison 2000/1-2009",
         xlab="CARDAMOM_NORAINFOR Output Wood g.m-2.d-1 ", ylab="RAINFOR Wood mortality g.m-2.d-1", pch=19, xlim=c(0,4), ylim=c(0,4))
    abline(coef = c(0,1),col='red', lwd=3)
  }
}
res_df_merge_plot_outputwood_2 <- function (benchmark,model,model_type){
  model_res <- resample(model,benchmark)
  model_res_df <- as.data.frame(model_res, xy=TRUE)
  benchmark_df <- as.data.frame(benchmark, xy=TRUE)
  model_benchmark_df <- merge(benchmark_df, model_res_df, by=c("x","y"))
  names(model_benchmark_df)<-c('x','y','rainfor','cardamom')
  model_benchmark_df <- model_benchmark_df[!is.na(model_benchmark_df$rainfor)&!is.na(model_benchmark_df$cardamom),]
  if (model_type=='rainfor') {
    plot(model_benchmark_df$cardamom, model_benchmark_df$rainfor, main="Woody Carbon losses comparison 2010-2016",
         xlab="CARDAMOM_RAINFOR Output Wood g.m-2.d-1 ", ylab="RAINFOR Wood mortality g.m-2.d-1", pch=19, xlim=c(0,4), ylim=c(0,4))
    abline(coef = c(0,1),col='red', lwd=3)
  }
  else {
    plot(model_benchmark_df$cardamom, model_benchmark_df$rainfor, main="Woody Carbon losses comparison 2010-2016",
         xlab="CARDAMOM_NORAINFOR Output Wood g.m-2.d-1 ", ylab="RAINFOR Wood mortality g.m-2.d-1", pch=19, xlim=c(0,4), ylim=c(0,4))
    abline(coef = c(0,1),col='red', lwd=3)
  }
}

########################################################################
####new code to include variable, model and region######################
########################################################################
#extract subset of wood, npp and mortality from rainfor data
extract_subset <- function (region,reference){
  masked1 <- mask(reference, region)
  masked2 <- masked1 > -Inf
  poly <- rasterToPolygons(masked2, dissolve=TRUE)
  data_region <- mask(reference,poly)
  return(data_region)
}

#extract subset of wood, npp and mortality from models
extract_biomass_subset_all_models <- function (region,cardamom_var,model_variant,reference){
     for (i in model_variant) {
         if (i=='norainfor') {
             data <- brick(paste(prefix,i,suffix,sep=""),varname=cardamom_var)
             data <- crop(data, extent(reference))
             a <- data$X2001.01.01
             cropped <- crop(a, extent(region))
             masked1 <- mask(cropped, region)
             masked2 <- masked1 > -Inf
             poly <- rasterToPolygons(masked2, dissolve=TRUE)
             data_region <- mask(data,poly)
             data_2010<-data_region[[which(getZ(data_region) >= as.Date("2010-01-01") & getZ(data_region) <= as.Date("2010-12-01"))]]
             data_2010_mean <- stackApply(data_2010, indices =  rep(1,nlayers(data_2010)), fun = "mean")
             }
         else {
             data <- brick(paste(prefix,i,suffix,sep=""),varname=cardamom_var)
             data_2010<-data[[which(getZ(data) >= as.Date("2010-01-01") & getZ(data) <= as.Date("2010-12-01"))]]
             data_2010_mean <- stackApply(data_2010, indices =  rep(1,nlayers(data_2010)), fun = "mean")
           }
        }
     return(data_2010_mean)
}

extract_biomass_subset_all_models(amazonia_subset,compare_var_names[1],mod_var[1],biomass_amazon_gCm2)

collate_all_models_cwood <- function (region,cardamom_var,model_variant,reference){
  for (i in model_variant) {
    if (i=='norainfor') {
      cwood_norainfor <- extract_biomass_subset_all_models(region,cardamom_var,i,reference)
    }
    else if (i=='biomass') {
      cwood_biomass <- extract_biomass_subset_all_models(region,cardamom_var,i,reference)
    }
    else if (i=='productivity') {
      cwood_productivity <- extract_biomass_subset_all_models(region,cardamom_var,i,reference)
    }
    else if (i=='mortality') {
      cwood_mortality <- extract_biomass_subset_all_models(region,cardamom_var,i,reference)
    }
  }
  cwood_list <-c(cwood_norainfor,cwood_biomass,cwood_productivity,cwood_mortality)
  names(cwood_list) <- model_variant
  return(cwood_list)
}
collate_all_models_cwood(amazonia_subset,compare_var_names[1],mod_var,biomass_amazon_gCm2)

res_df_merge_plot_cwood <- function (benchmark,model,model_variant){
  model_res <- resample(model,benchmark)
  model_res_df <- as.data.frame(model_res, xy=TRUE)
  benchmark_df <- as.data.frame(benchmark, xy=TRUE)
  model_benchmark_df <- merge(benchmark_df, model_res_df, by=c("x","y"))
  names(model_benchmark_df)<-c('x','y','rainfor','cardamom')
  model_benchmark_df <- model_benchmark_df[!is.na(model_benchmark_df$rainfor)&!is.na(model_benchmark_df$cardamom),]
  if (model_variant=='norainfor') {
    plot(model_benchmark_df$cardamom, model_benchmark_df$rainfor, main="AGB comparison 2001-2010",
         xlab="CARDAMOM_RAINFOR Wood biomass g.m-2 ", ylab="RAINFOR Wood biomass g.m-2", pch=19, xlim=c(0,25000), ylim=c(0,25000))
    abline(coef = c(0,1),col='red', lwd=3)
  }
  if (model_variant=='biomass') {
    plot(model_benchmark_df$cardamom, model_benchmark_df$rainfor, main="AGB comparison 2001-2010",
         xlab="CARDAMOM_RAINFOR Wood biomass g.m-2 ", ylab="RAINFOR Wood biomass g.m-2", pch=19, xlim=c(0,25000), ylim=c(0,25000))
    abline(coef = c(0,1),col='red', lwd=3)
  }
  if (model_variant=='rainfor') {
    plot(model_benchmark_df$cardamom, model_benchmark_df$rainfor, main="AGB comparison 2001-2010",
         xlab="CARDAMOM_RAINFOR Wood biomass g.m-2 ", ylab="RAINFOR Wood biomass g.m-2", pch=19, xlim=c(0,25000), ylim=c(0,25000))
    abline(coef = c(0,1),col='red', lwd=3)
  }
  else {
    plot(model_benchmark_df$cardamom, model_benchmark_df$rainfor, main="AGB comparison 2001-2010",
         xlab="CARDAMOM_NORAINFOR Wood biomass g.m-2 ", ylab="RAINFOR Wood biomass g.m-2", pch=19, xlim=c(0,25000), ylim=c(0,25000))
    abline(coef = c(0,1),col='red', lwd=3)
  }
}

extract_npp_mort_subset_models <- function (region,cardamom_var,model_variant,reference,tp){
  for (i in model_variant) {
    if (i=='norainfor') {
      data <- brick(paste(prefix,i,suffix,sep=""),varname=cardamom_var)
      data <- crop(data, extent(reference))
      a <- data$X2001.01.01
      cropped <- crop(a, extent(region))
      masked1 <- mask(cropped, region)
      masked2 <- masked1 > -Inf
      poly <- rasterToPolygons(masked2, dissolve=TRUE)
      data_region <- mask(data,poly)
      data_0109<-data_region[[which(getZ(data_region) >= as.Date("2001-01-01") & getZ(data_region) <= as.Date("2009-12-01"))]]
      data_1016<-data_region[[which(getZ(data_region) >= as.Date("2010-01-01") & getZ(data_region) <= as.Date("2016-12-01"))]]
      data_0109_mean <- stackApply(data_0109, indices =  rep(1,nlayers(data_0109)), fun = "mean")
      data_1016_mean <- stackApply(data_1016, indices =  rep(1,nlayers(data_1016)), fun = "mean")
      }
    else {
      data <- brick(paste(prefix,model_variant[1],suffix,sep=""),varname=cardamom_var)
      data_0109<-data[[which(getZ(data) >= as.Date("2001-01-01") & getZ(data) <= as.Date("2009-12-01"))]]
      data_1016<-data[[which(getZ(data) >= as.Date("2010-01-01") & getZ(data) <= as.Date("2016-12-01"))]]
      data_0109_mean <- stackApply(data_0109, indices =  rep(1,nlayers(data_0109)), fun = "mean")
      data_1016_mean <- stackApply(data_1016, indices =  rep(1,nlayers(data_1016)), fun = "mean")
    }
    }
  return(c(data_0109_mean,data_1016_mean))
}

extract_npp_mort_subset_models(amazonia_subset,compare_var_names[2],mod_var[1],woodprod_00_09_gCm2d,time_period)

collate_all_models_npp_mort <- function (region,cardamom_var,model_variant,reference,tp){
  for (i in model_variant) {
        if (i=='norainfor') {
          nppwood_norainfor <- extract_npp_mort_subset_models(region,cardamom_var,i,reference,j)
        }
        else if (i=='biomass') {
          nppwood_biomass <- extract_npp_mort_subset_models(region,cardamom_var,i,reference,j)
        }
        else if (i=='productivity') {
          nppwood_productivity <- extract_npp_mort_subset_models(region,cardamom_var,i,reference,j)
        }
        else if (i=='mortality') {
          nppwood_mortality <- extract_npp_mort_subset_models(region,cardamom_var,i,reference,j)
        }
  }
  for (j in tp) {
    if (j=='01/09'){
      nppwood_list <-c(nppwood_norainfor[1],nppwood_biomass[1],nppwood_productivity[1],nppwood_mortality[1])
      names(nppwood_list) <- model_variant
      return(nppwood_list)
    }
    else if (j=='10/16') {
      nppwood_list <-c(nppwood_norainfor[2],nppwood_biomass[2],nppwood_productivity[2],nppwood_mortality[2])
      names(nppwood_list) <- model_variant
      return(nppwood_list)
    }
  }
}
collate_all_models_npp_mort(amazonia_subset,compare_var_names[2],mod_var,woodprod_00_09_gCm2d,time_period[2])

######extract for each variable#########

reg_card_var_2 <- function(s,s2,cardamom_variables){
  for (i in cardamom_variables){
    if (i=='WOOD') {
      cwood_data <- extract_var_mean_biom_subset(s,s2,i)}
    else if (i=='LAI') {
      lai_reg_data <- amazon_reg_to_card_var_2(s,s2,i)}
    else if (i=='CICA') {
      cica_reg_data <- amazon_reg_to_card_var_2(s,s2,i)}
    else if (i=='GSDSR') {
      gsdsr_reg_data <- amazon_reg_to_card_var_2(s,s2,i)}
    else if (i=='APAR') {
      apar_reg_data <- amazon_reg_to_card_var_2(s,s2,i)}
    else{
      et_reg_data <- amazon_reg_to_card_var_2(s,s2,i)
    }
  }
  
  region_name <- as.factor(rep(s2,nrow(gpp_reg_data)))
  no_dp<-nrow(gpp_reg_data)/228
  date<- as.Date(rep(seq(from = as.Date("2001-01-01"), to = as.Date("2019-12-01"), by = 'month'),no_dp))
  month<-format(date,format="%b")
  year <- format(date,format="%Y")
  
  return(na.omit(cbind(gpp_reg_data,lai_reg_data,cica_reg_data,gsdsr_reg_data,apar_reg_data,et_reg_data,region_name,date,month,year)))
}
  
  