#######################################################
###############all extraction code#####################
#######################################################
######function to extract country amazon regions#######
amazon_reg_to_card_var_2 <- function(s,s2,t) {
  r <- brick('R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/DATA/benchmark/NoRainfor_Amazon_1deg_monthly_2001_updated_2019_compare.nc',varname=t)
  a <- r$X2001.01.01
  cropped <- crop(a, extent(s))
  masked1 <- mask(cropped, s)
  masked2 <- masked1 > -Inf
  poly <- rasterToPolygons(masked2, dissolve=TRUE)
  card_region <- mask(r,poly)
  df_card_region <- as.data.frame(card_region,xy=T)
  df_card_region<-as.data.frame(t(df_card_region))
  df_card_region <- data.frame(y=unlist(df_card_region))
  names(df_card_region) <- paste(tolower(t))
  
  return(df_card_region)
}
# gpp_amazonia <- amazon_reg_to_card_var_2(reg_data[[6]],reg_names[6],var_names[1])

amazon_to_card_var_2 <- function(s2,t) {
  card_region <- brick('R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/DATA/benchmark/NoRainfor_Amazon_1deg_monthly_2001_updated_2019_compare.nc',varname=t)
  df_card_region <- as.data.frame(card_region,xy=F)
  df_card_region<-as.data.frame(t(df_card_region))
  df_card_region <- data.frame(y=unlist(df_card_region))
  names(df_card_region) <- paste(tolower(t))
  
  return(df_card_region)
}

amazonia_reg_to_card_regression <- function(s,s2,t) {
  r <- brick('R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/DATA/benchmark/NoRainfor_Amazon_1deg_monthly_2001_updated_2019_compare.nc',varname=t)
  a <- r$X2001.01.01
  cropped <- crop(a, extent(s))
  masked1 <- mask(cropped, s)
  masked2 <- masked1 > -Inf
  poly <- rasterToPolygons(masked2, dissolve=TRUE)
  card_region <- mask(r,poly)
  df_card <- as.data.frame(card_region,xy=T)
  df_card_region <- as.matrix(as.data.frame(card_region,xy=F))
  x_cord <- matrix(df_card$x, nrow=nrow(df_card_region), ncol=ncol(df_card_region))
  y_cord <- matrix(df_card$y, nrow=nrow(df_card_region), ncol=ncol(df_card_region))
  cell <- matrix(1:nrow(df_card_region), nrow=nrow(df_card_region), ncol=ncol(df_card_region))
  t_date <- matrix(names(as.data.frame(df_card_region)), nrow=nrow(df_card_region), ncol=ncol(df_card_region), byrow=TRUE)
  df_card_region <- data.frame(x=c(x_cord),y=c(y_cord), t=c(t_date), cell=factor(sprintf("C%04d",c(cell))),variable=c(df_card_region))
  names(df_card_region) <- c('x','y','date','cell_no',paste(tolower(t)))
  
  return(df_card_region)
}
# summary(amazonia_reg_to_card_regression(reg_data[[3]],reg_names[3],var_names[1]))

######extract for each variable#########
reg_card_var_2 <- function(s,s2,cardamom_variables){
  for (i in cardamom_variables){
    if (i=='GPP') {
      gpp_reg_data <- amazon_reg_to_card_var_2(s,s2,i)}
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

card_var_2 <- function(s2,cardamom_variables){
  for (i in cardamom_variables){
    if (i=='GPP') {
      gpp_amazon_data <- amazon_to_card_var_2(s2,i)}
    else if (i=='LAI') {
      lai_amazon_data <- amazon_to_card_var_2(s2,i)}
    else{
      cica_amazon_data <- amazon_to_card_var_2(s2,i)
    }
  }
  region_name <- as.factor(rep(s2,nrow(gpp_amazon_data)))
  no_dp<-nrow(gpp_amazon_data)/228
  date<- as.Date(rep(seq(from = as.Date("2001-01-01"), to = as.Date("2019-12-01"), by = 'month'),no_dp))
  month<-format(date,format="%b")
  year <- format(date,format="%Y")
  
  return(na.omit(cbind(gpp_amazon_data,lai_amazon_data,cica_amazon_data,region_name,date,month,year)))
}

reg_card_var_regression <- function(s,s2,cardamom_variables){
  for (i in cardamom_variables){
    if (i=='GPP') {
      gpp_reg_data <- amazonia_reg_to_card_regression(s,s2,i)}
    else if (i=='LAI') {
      lai_reg_data <- amazonia_reg_to_card_regression(s,s2,i)}
    else if (i=='CICA') {
      cica_reg_data <- amazonia_reg_to_card_regression(s,s2,i)}
    else if (i=='GSDSR') {
      gsdsr_reg_data <- amazonia_reg_to_card_regression(s,s2,i)}
    else if (i=='APAR') {
      apar_reg_data <- amazonia_reg_to_card_regression(s,s2,i)}
    else{
      et_reg_data <- amazonia_reg_to_card_regression(s,s2,i)
    }
  }
  region_name <- as.factor(rep(s2,nrow(gpp_reg_data)))
  df_card_regres <- cbind(gpp_reg_data,lai_reg_data,cica_reg_data,gsdsr_reg_data,apar_reg_data,et_reg_data)
  return(df_card_regres[,c(1:5,10,15,20,25,30)])
  #return(na.omit(cbind(gpp_reg_data,lai_reg_data,cica_reg_data,gsdsr_reg_data,apar_reg_data,et_reg_data)))
  #return(cbind(gpp_reg_data,lai_reg_data,cica_reg_data,gsdsr_reg_data,apar_reg_data,et_reg_data))
  }
# summary(reg_card_var_regression(reg_data[[6]],reg_names[6],var_names))

######Run for all regions and all variables###########
all_reg_run <- function(s,s2,cardamom_variables) {
  for (i in s2) {
    if(i=='amazon_nw'){
      nw_data <- reg_card_var(reg_data[[1]],i,var_names)
    }
    else if (i=='amazon_sw') {
      sw_data <- reg_card_var(reg_data[[2]],i,var_names)
    }
    else if (i=='amazon_ec') {
      ec_data <- reg_card_var(reg_data[[3]],i,var_names)
    }
    else if (i=='amazon_bs') {
      bs_data <- reg_card_var(reg_data[[4]],i,var_names)
    }
    else if (i=='amazon_gs') {
      gs_data <- reg_card_var(reg_data[[5]],i,var_names)
    }
    else{
      amazon_data <- card_var(i,var_names)
    }
  }
  return(rbind(nw_data,sw_data,ec_data,bs_data,gs_data,amazon_data))
}

all_reg_run_2 <- function(s,s2,cardamom_variables) {
  for (i in s2) {
    if(i=='amazon_nw'){
      nw_data <- reg_card_var_2(reg_data[[1]],i,var_names)
    }
    else if (i=='amazon_sw') {
      sw_data <- reg_card_var_2(reg_data[[2]],i,var_names)
    }
    else if (i=='amazon_ec') {
      ec_data <- reg_card_var_2(reg_data[[3]],i,var_names)
    }
    else if (i=='amazon_bs') {
      bs_data <- reg_card_var_2(reg_data[[4]],i,var_names)
    }
    else if (i=='amazon_gs') {
      gs_data <- reg_card_var_2(reg_data[[5]],i,var_names)
    }
    else{
      amazon_data <- reg_card_var_2(reg_data[[6]],i,var_names)
    }
  }
  return(rbind(nw_data,sw_data,ec_data,bs_data,gs_data,amazon_data))
}

all_reg_var_regres_run <- function(s,s2,cardamom_variables) {
  for (i in s2) {
    if(i=='amazon_nw'){
      nw_data <- reg_card_var_regression(reg_data[[1]],i,var_names)
    }
    else if (i=='amazon_sw') {
      sw_data <- reg_card_var_regression(reg_data[[2]],i,var_names)
    }
    else if (i=='amazon_ec') {
      ec_data <- reg_card_var_regression(reg_data[[3]],i,var_names)
    }
    else if (i=='amazon_bs') {
      bs_data <- reg_card_var_regression(reg_data[[4]],i,var_names)
    }
    else if (i=='amazon_gs') {
      gs_data <- reg_card_var_regression(reg_data[[5]],i,var_names)
    }
    else{
      amazon_data <- reg_card_var_regression(reg_data[[6]],i,var_names)
    }
  }
  return(rbind(nw_data,sw_data,ec_data,bs_data,gs_data))
}

##############################################
###########Pixel based analysis###############
##############################################
extractfun <- function(m) {
  tinfo <- summary(m)$sigma
  r2 <- summary(m)$adj.r.squared
  data.frame(Rsq = r2,Resid.se = tinfo)
}
run_all_var<- function(data,variables){
  for (i in variables){
    if (variables=='LAI'){
      split_data_lm <- dlply(na.omit(data), .(cell_no,x,y), function(s) lm(gpp~lai, data=s))
      gpp_lai_lm_summary<-ldply(split_data_lm, extractfun)
      pixels_gpp_lai_lm_rsq_raster <- rasterFromXYZ(gpp_lai_lm_summary[,2:4])
      pixels_gpp_lai_lm_rse_raster <- rasterFromXYZ(gpp_lai_lm_summary[,c(2:3,5)])
      crs(pixels_gpp_lai_lm_rsq_raster) <- '+proj=longlat +datum=WGS84 +no_defs'
      crs(pixels_gpp_lai_lm_rse_raster) <- '+proj=longlat +datum=WGS84 +no_defs'
      par(mfrow = c(1,2),oma = c(0, 0, 2, 0))
      plot(pixels_gpp_lai_lm_rsq_raster,main='R sq');plot(amazon_nw_poly,add=T);plot(amazon_sw_poly,add=T);plot(amazon_ec_poly,add=T);plot(amazon_bs_poly,add=T);plot(amazon_gs_poly,add=T)
      plot(pixels_gpp_lai_lm_rse_raster,main='Residual SE');plot(amazon_nw_poly,add=T);plot(amazon_sw_poly,add=T);plot(amazon_ec_poly,add=T);plot(amazon_bs_poly,add=T);plot(amazon_gs_poly,add=T)
      mtext(paste(i), line=-1, outer=TRUE, cex=1.5)
    }
    else if (variables=='CICA'){
      split_data_lm <- dlply(na.omit(data), .(cell_no,x,y), function(s) lm(gpp~cica, data=s))
      gpp_cica_lm_summary<-ldply(split_data_lm, extractfun)
      pixels_gpp_cica_lm_rsq_raster <- rasterFromXYZ(gpp_cica_lm_summary[,2:4])
      pixels_gpp_cica_lm_rse_raster <- rasterFromXYZ(gpp_cica_lm_summary[,c(2:3,5)])
      crs(pixels_gpp_cica_lm_rsq_raster) <- '+proj=longlat +datum=WGS84 +no_defs'
      crs(pixels_gpp_cica_lm_rse_raster) <- '+proj=longlat +datum=WGS84 +no_defs'
      par(mfrow = c(1,2),oma = c(0, 0, 2, 0))
      plot(pixels_gpp_cica_lm_rsq_raster,main='R sq');plot(amazon_nw_poly,add=T);plot(amazon_sw_poly,add=T);plot(amazon_ec_poly,add=T);plot(amazon_bs_poly,add=T);plot(amazon_gs_poly,add=T)
      plot(pixels_gpp_cica_lm_rse_raster,main='Residual SE');plot(amazon_nw_poly,add=T);plot(amazon_sw_poly,add=T);plot(amazon_ec_poly,add=T);plot(amazon_bs_poly,add=T);plot(amazon_gs_poly,add=T)
      mtext(paste(i), line=-1, outer=TRUE, cex=1.5)
    }
    else if (variables=='GSDSR'){
      split_data_lm <- dlply(na.omit(data), .(cell_no,x,y), function(s) lm(gpp~gsdsr, data=s))
      gpp_gsdsr_lm_summary<-ldply(split_data_lm, extractfun)
      pixels_gpp_gsdsr_lm_rsq_raster <- rasterFromXYZ(gpp_gsdsr_lm_summary[,2:4])
      pixels_gpp_gsdsr_lm_rse_raster <- rasterFromXYZ(gpp_gsdsr_lm_summary[,c(2:3,5)])
      crs(pixels_gpp_gsdsr_lm_rsq_raster) <- '+proj=longlat +datum=WGS84 +no_defs'
      crs(pixels_gpp_gsdsr_lm_rse_raster) <- '+proj=longlat +datum=WGS84 +no_defs'
      par(mfrow = c(1,2),oma = c(0, 0, 2, 0))
      plot(pixels_gpp_gsdsr_lm_rsq_raster,main='R sq');plot(amazon_nw_poly,add=T);plot(amazon_sw_poly,add=T);plot(amazon_ec_poly,add=T);plot(amazon_bs_poly,add=T);plot(amazon_gs_poly,add=T)
      plot(pixels_gpp_gsdsr_lm_rse_raster,main='Residual SE');plot(amazon_nw_poly,add=T);plot(amazon_sw_poly,add=T);plot(amazon_ec_poly,add=T);plot(amazon_bs_poly,add=T);plot(amazon_gs_poly,add=T)
      mtext(paste(i), line=-1, outer=TRUE, cex=1.5)
    }
    else if (variables=='APAR'){
      split_data_lm <- dlply(na.omit(data), .(cell_no,x,y), function(s) lm(gpp~apar, data=s))
      gpp_apar_lm_summary<-ldply(split_data_lm, extractfun)
      pixels_gpp_apar_lm_rsq_raster <- rasterFromXYZ(gpp_apar_lm_summary[,2:4])
      pixels_gpp_apar_lm_rse_raster <- rasterFromXYZ(gpp_apar_lm_summary[,c(2:3,5)])
      crs(pixels_gpp_apar_lm_rsq_raster) <- '+proj=longlat +datum=WGS84 +no_defs'
      crs(pixels_gpp_apar_lm_rse_raster) <- '+proj=longlat +datum=WGS84 +no_defs'
      par(mfrow = c(1,2),oma = c(0, 0, 2, 0))
      plot(pixels_gpp_apar_lm_rsq_raster,main='R sq');plot(amazon_nw_poly,add=T);plot(amazon_sw_poly,add=T);plot(amazon_ec_poly,add=T);plot(amazon_bs_poly,add=T);plot(amazon_gs_poly,add=T)
      plot(pixels_gpp_apar_lm_rse_raster,main='Residual SE');plot(amazon_nw_poly,add=T);plot(amazon_sw_poly,add=T);plot(amazon_ec_poly,add=T);plot(amazon_bs_poly,add=T);plot(amazon_gs_poly,add=T)
      mtext(paste(i), line=-1, outer=TRUE, cex=1.5)
    }
    else if (variables=='ET'){
      split_data_lm <- dlply(na.omit(data), .(cell_no,x,y), function(s) lm(gpp~et, data=s))
      gpp_et_lm_summary<-ldply(split_data_lm, extractfun)
      pixels_gpp_et_lm_rsq_raster <- rasterFromXYZ(gpp_et_lm_summary[,2:4])
      pixels_gpp_et_lm_rse_raster <- rasterFromXYZ(gpp_et_lm_summary[,c(2:3,5)])
      crs(pixels_gpp_et_lm_rsq_raster) <- '+proj=longlat +datum=WGS84 +no_defs'
      crs(pixels_gpp_et_lm_rse_raster) <- '+proj=longlat +datum=WGS84 +no_defs'
      par(mfrow = c(1,2),oma = c(0, 0, 2, 0))
      plot(pixels_gpp_et_lm_rsq_raster,main='R sq');plot(amazon_nw_poly,add=T);plot(amazon_sw_poly,add=T);plot(amazon_ec_poly,add=T);plot(amazon_bs_poly,add=T);plot(amazon_gs_poly,add=T)
      plot(pixels_gpp_et_lm_rse_raster,main='Residual SE');plot(amazon_nw_poly,add=T);plot(amazon_sw_poly,add=T);plot(amazon_ec_poly,add=T);plot(amazon_bs_poly,add=T);plot(amazon_gs_poly,add=T)
      mtext(paste(i), line=-1, outer=TRUE, cex=1.5)
    }
  }
}