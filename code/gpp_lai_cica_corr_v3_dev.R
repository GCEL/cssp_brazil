#####################################################################
##########################Compare GPP, LAI and CiCa##################
#####################################################################
######Load required packages#######
library(ncdf4); library(raster); library(dplyr); library(ggplot2);library(ggpubr);library(quantreg);library(ggpp);library(rgeos);library(ggpmisc);library(rgdal)
#end
######Read in shape files ###########
amazonia_nw <- shapefile("./data/amazon_nw.shp");amazonia_sw <- shapefile("./data/amazon_sw.shp");amazonia_ec <- shapefile("./data/amazon_ec.shp")
amazonia_bs <- shapefile("./data/amazon_bs.shp");amazonia_gs <- shapefile("./data/amazon_gs.shp")
amazonia<- biomass_amazon_mask_pol

var_names <- c('GPP','LAI','CICA','GSDSR','APAR','ET')
reg_data <- list(amazonia_nw,amazonia_sw,amazonia_ec,amazonia_bs,amazonia_gs,amazonia)
reg_names <- c('amazon_nw','amazon_sw','amazon_ec','amazon_bs','amazon_gs','amazonia')
#names(reg_data) <- reg_names

########extract subset for analysis----
amazonia_subset <- shapefile("./data/amazonia_subset.shp")

extract_subset <- function (region,reference){
  masked1 <- mask(reference, region)
  masked2 <- masked1 > -Inf
  poly <- rasterToPolygons(masked2, dissolve=TRUE)
  data_region <- mask(reference,poly)
  return(data_region)
}

amazonia_subset_raster<-extract_subset(amazonia_subset,biomass_amazon_gCm2)
#plot(amazonia_subset_raster)
writeRaster(amazonia_subset_raster, filename = "G://cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/amazonia_subset.tif", format = "GTiff")

#done
######function to extract country amazon regions#######
amazon_reg_to_card_var_2 <- function(s,s2,t) {
  r <- brick('R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/DATA/benchmark/NoRainfor_Amazon_1deg_monthly_2001_updated_2019_compare.nc',varname=t)
  a <- r$X2001.01.01
  cropped <- crop(a, extent(s))
  masked1 <- mask(cropped, s)
  masked2 <- masked1 > -Inf
  poly <- rasterToPolygons(masked2, dissolve=TRUE)
  card_region <- mask(r,poly)
  df_card_region <- as.data.frame(card_region,xy=F)
  df_card_region<-as.data.frame(t(df_card_region))
  df_card_region <- data.frame(y=unlist(df_card_region))
  names(df_card_region) <- paste(tolower(t))
  
  return(df_card_region)
}
#cica_nw <- amazon_reg_to_card_var_2(reg_data[[1]],reg_names[1],var_names[1])

amazon_to_card_var_2 <- function(s2,t) {
  card_region <- brick('R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/DATA/benchmark/NoRainfor_Amazon_1deg_monthly_2001_updated_2019_compare.nc',varname=t)
  df_card_region <- as.data.frame(card_region,xy=F)
  df_card_region<-as.data.frame(t(df_card_region))
  df_card_region <- data.frame(y=unlist(df_card_region))
  names(df_card_region) <- paste(tolower(t))
  
  return(df_card_region)
}
#cica_amazon <- amazon_to_card_var_2(reg_names[6],var_names[3])

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
#nw_data <- reg_card_var_2(reg_data[[1]],reg_names[1],var_names)
#nrow(nw_data)

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
#amazon_data <- card_var_2(reg_names[6],var_names)

######Run for all regions and all variables###########
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

######Run codes for extraction#####
all_pixels_compare_data <- all_reg_run_2(reg_data,reg_names,var_names)
rownames(all_pixels_compare_data) <- 1:nrow(all_pixels_compare_data)
all_pixels_compare_data_zero <- all_pixels_compare_data[apply(all_pixels_compare_data, 1, function(row) all(row !=0 )), ]  # Remove zero-rows
row_sub <- apply(all_pixels_compare_data_zero, 1, function(row) all(row !=0 ))
all_pixels_compare_data_zero <- all_pixels_compare_data_zero[row_sub,]

all_pixels_compare_data_zero[all_pixels_compare_data_zero==0] <- NA
all_pixels_compare_data_zero<-na.omit(all_pixels_compare_data_zero)

summary(all_pixels_compare_data);names(all_pixels_compare_data);head(all_pixels_compare_data)
summary(all_pixels_compare_data_zero);names(all_pixels_compare_data_zero);head(all_pixels_compare_data_zero);str(all_pixels_compare_data_zero)
all_pixels_compare_data_zero$month_f<-factor(all_pixels_compare_data_zero$month,levels=c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'))
#all_pixels_compare_data_zero$region_name_f<-factor(all_pixels_compare_data_zero$region_name,levels=c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'))
saveRDS(all_pixels_compare_data_zero,"./data/extracted_variables.RDS")
#done
######Time Series Plot######
#monthly mean
par(mfrow=c(3,1))
boxplot(all_pixels_compare_data$lai~all_pixels_compare_data$month,ylab='LAI',xlab='Month')
boxplot(all_pixels_compare_data$gpp~all_pixels_compare_data$month,ylab='GPP',xlab='Month',ylim=c(6.25,10))
boxplot(all_pixels_compare_data$cica~all_pixels_compare_data$month,ylab='CiCa',xlab='Month',ylim=c(0.6,0.8))
mtext("North West",side=3,line=-2,outer=TRUE, font=2)

boxplot(all_pixels_compare_data$lai_amazon_sw~all_pixels_compare_data$month,ylab='LAI',xlab='Month')
boxplot(all_pixels_compare_data$gpp_amazon_sw~all_pixels_compare_data$month,ylab='GPP',xlab='Month',ylim=c(6.25,10))
boxplot(all_pixels_compare_data$cica_amazon_sw~all_pixels_compare_data$month,ylab='CiCa',xlab='Month',ylim=c(0.6,0.8))
mtext("South West",side=3,line=-2,outer=TRUE, font=2)

boxplot(all_pixels_compare_data$lai_amazon_ec~all_pixels_compare_data$month,ylab='LAI',xlab='Month')
boxplot(all_pixels_compare_data$gpp_amazon_ec~all_pixels_compare_data$month,ylab='GPP',xlab='Month',ylim=c(6.25,10))
boxplot(all_pixels_compare_data$cica_amazon_ec~all_pixels_compare_data$month,ylab='CiCa',xlab='Month',ylim=c(0.6,0.8))
mtext("East Central",side=3,line=-2,outer=TRUE, font=2)

boxplot(all_pixels_compare_data$lai_amazon_bs~all_pixels_compare_data$month,ylab='LAI',xlab='Month')
boxplot(all_pixels_compare_data$gpp_amazon_bs~all_pixels_compare_data$month,ylab='GPP',xlab='Month',ylim=c(6.25,10))
boxplot(all_pixels_compare_data$cica_amazon_bs~all_pixels_compare_data$month,ylab='CiCa',xlab='Month',ylim=c(0.6,0.8))
mtext("Brazil Shield",side=3,line=-2,outer=TRUE, font=2)

boxplot(all_pixels_compare_data$lai_amazon_gs~all_pixels_compare_data$month,ylab='LAI',xlab='Month')
boxplot(all_pixels_compare_data$gpp_amazon_gs~all_pixels_compare_data$month,ylab='GPP',xlab='Month',ylim=c(6.25,10))
boxplot(all_pixels_compare_data$cica_amazon_gs~all_pixels_compare_data$month,ylab='CiCa',xlab='Month',ylim=c(0.6,0.8))
mtext("Guiana Shield",side=3,line=-2,outer=TRUE, font=2)

boxplot(all_pixels_compare_data$lai_amazonia~all_pixels_compare_data$month,ylab='LAI',xlab='Month')
boxplot(all_pixels_compare_data$gpp_amazonia~all_pixels_compare_data$month,ylab='GPP',xlab='Month',ylim=c(6.25,10))
boxplot(all_pixels_compare_data$cica_amazonia~all_pixels_compare_data$month,ylab='CiCa',xlab='Month',ylim=c(0.6,0.8))
mtext("Amazonia",side=3,line=-2,outer=TRUE, font=2)


gpplai_amazon_nw<-summary(lm(gpp ~ lai, data = all_pixels_compare_data))
rsq_gpplai_amazon_nw<-paste0("R.sp = ", round(gpplai$r.squared,2))
gppcica_amazon_nw<-summary(lm(gpp ~ cica, data = all_pixels_compare_data))
rsq_gppcica_amazon_nw<-paste0("R.sp = ", round(gppcica$r.squared,2))
      
gpplai_amazon_sw<-summary(lm(gpp_amazon_sw ~ lai_amazon_sw, data = all_pixels_compare_data))
rsq_gpplai_amazon_sw<-paste0("R.sp = ", round(gpplai_amazon_sw$r.squared,2))
gppcica_amazon_sw<-summary(lm(gpp_amazon_sw ~ cica_amazon_sw, data = all_pixels_compare_data))
rsq_gppcica_amazon_sw<-paste0("R.sp = ", round(gppcica_amazon_sw$r.squared,2))

gpplai_amazon_ec<-summary(lm(gpp_amazon_ec ~ lai_amazon_ec, data = all_pixels_compare_data))
rsq_gpplai_amazon_ec<-paste0("R.sp = ", round(gpplai_amazon_ec$r.squared,2))
gppcica_amazon_ec<-summary(lm(gpp_amazon_ec ~ cica_amazon_ec, data = all_pixels_compare_data))
rsq_gppcica_amazon_ec<-paste0("R.sp = ", round(gppcica_amazon_ec$r.squared,2))
    
gpplai_amazon_bs<-summary(lm(gpp_amazon_bs ~ lai_amazon_bs, data = all_pixels_compare_data))
rsq_gpplai_amazon_bs<-paste0("R.sp = ", round(gpplai_amazon_bs$r.squared,2))
gppcica_amazon_bs<-summary(lm(gpp_amazon_bs ~ cica_amazon_bs, data = all_pixels_compare_data))
rsq_gppcica_amazon_bs<-paste0("R.sp = ", round(gppcica_amazon_bs$r.squared,2))
    
gpplai_amazon_gs<-summary(lm(gpp_amazon_gs ~ lai_amazon_gs, data = all_pixels_compare_data))
rsq_gpplai_amazon_gs<-paste0("R.sp = ", round(gpplai_amazon_gs$r.squared,2))
gppcica_amazon_gs<-summary(lm(gpp_amazon_gs ~ cica_amazon_gs, data = all_pixels_compare_data))
rsq_gppcica_amazon_gs<-paste0("R.sp = ", round(gppcica_amazon_gs$r.squared,2))
    
gpplai_amazonia<-summary(lm(gpp_amazonia ~ lai_amazonia, data = all_pixels_compare_data))
rsq_gpplai_amazonia<-paste0("R.sp = ", round(gpplai_amazonia$r.squared,2))
gppcica_amazonia<-summary(lm(gpp_amazonia ~ cica_amazonia, data = all_pixels_compare_data))
rsq_gppcica_amazonia<-paste0("R.sp = ", round(gppcica_amazonia$r.squared,2))
    
#rsq_amazon_rltnshps<-rsq_relationships(all_pixels_compare_data,reg_names)

# gpplai<-summary(lm(gpp_mean ~ lai_mean, data = all_pixels_compare_data))
# rsq_gpplai<-paste0("R.sp = ", round(gpplai$r.squared,2))
# gppcica<-summary(lm(gpp_mean ~ cica_mean, data = all_pixels_compare_data))
# rsq_gppcica<-paste0("R.sp = ", round(gppcica$r.squared,2))

all_pixels_compare_data<-as.data.frame(all_pixels_compare_data)
par(mfrow=c(1,2))
plot(all_pixels_compare_data$gpp~all_pixels_compare_data$lai,xlab='LAI', ylab='GPP', main=rsq_gpplai_amazon_nw)
plot(all_pixels_compare_data$gpp~all_pixels_compare_data$cica, main=rsq_gppcica_amazon_nw,ylab='GPP',xlab='CiCa')
mtext("North West",side=3,line=-3,outer=TRUE, font=2,cex=1.5)

plot(all_pixels_compare_data$gpp_amazon_sw~all_pixels_compare_data$lai_amazon_sw, main=rsq_gpplai_amazon_sw,ylab='GPP',xlab='LAI')
plot(all_pixels_compare_data$gpp_amazon_sw~all_pixels_compare_data$cica_amazon_sw, main=rsq_gppcica_amazon_sw,ylab='GPP',xlab='CiCa')
mtext("South West",side=3,line=-3,outer=TRUE, font=2,cex=1.5)

plot(all_pixels_compare_data$gpp_amazon_ec~all_pixels_compare_data$lai_amazon_ec, main=rsq_gpplai_amazon_ec,ylab='GPP',xlab='LAI')
plot(all_pixels_compare_data$gpp_amazon_ec~all_pixels_compare_data$cica_amazon_ec, main=rsq_gppcica_amazon_ec,ylab='GPP',xlab='CiCa')
mtext("East Central",side=3,line=-3,outer=TRUE, font=2,cex=1.5)

plot(all_pixels_compare_data$gpp_amazon_bs~all_pixels_compare_data$lai_amazon_bs, main=rsq_gpplai_amazon_bs,ylab='GPP',xlab='LAI')
plot(all_pixels_compare_data$gpp_amazon_bs~all_pixels_compare_data$cica_amazon_bs, main=rsq_gppcica_amazon_bs,ylab='GPP',xlab='CiCa')
mtext("Brazil Shield",side=3,line=-3,outer=TRUE, font=2,cex=1.5)

plot(all_pixels_compare_data$gpp_amazon_gs~all_pixels_compare_data$lai_amazon_gs, main=rsq_gpplai_amazon_gs,ylab='GPP',xlab='LAI')
plot(all_pixels_compare_data$gpp_amazon_gs~all_pixels_compare_data$cica_amazon_gs, main=rsq_gppcica_amazon_gs,ylab='GPP',xlab='CiCa')
mtext("Guiana Shield",side=3,line=-3,outer=TRUE, font=2,cex=1.5)

plot(all_pixels_compare_data$gpp_amazonia~all_pixels_compare_data$lai_amazonia, main=rsq_gpplai_amazonia,ylab='GPP',xlab='LAI')
plot(all_pixels_compare_data$gpp_amazonia~all_pixels_compare_data$cica_amazonia, main=rsq_gppcica_amazonia,ylab='GPP',xlab='CiCa')
mtext("Amazonia",side=3,line=-3,outer=TRUE, font=2,cex=1.5)


#####multiple plots
formula <- y ~ x

all_pixels_compare_data %>%
  ggplot(aes(x=lai, 
             y=gpp, 
             shape = region_name))+
  geom_point(alpha = 0.2)+
  geom_smooth(method="lm", se = FALSE, fill = NA,fullrange=TRUE,aes(color=region_name))+
  theme_bw() + 
  theme(text = element_text(size = 16),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#ggsave("scatterplot_with_multiple_groups_ggplot2.png")

all_pixels_compare_data %>%
  ggplot(aes(x=cica, 
             y=gpp, 
             shape = region_name))+
  geom_point(alpha = 0.2)+
  geom_smooth(method="lm", se = FALSE, fill = NA,fullrange=TRUE,aes(color=region_name))+
  theme_bw() + 
  theme(text = element_text(size = 16),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#ggsave("scatterplot_with_multiple_groups_ggplot2.png")

dfg <- data.frame(x = c(1:100))
dfg$y <- 20 * c(0, 1) + 3 * dfg$x + rnorm(100, sd = 40)
dfg$group <- factor(rep(c("A", "B"), 50))

ggplot(data = dfg, aes(x = x, y = y, colour = group)) +
  stat_poly_line() +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*"))) +
  geom_point()

###################new plots##################
all_pixels_compare_data_zero %>%
  ggplot(aes(x=lai, 
             y=gpp, 
             colour = region_name))+
  geom_point(alpha = 0.2)+
  labs(title="Plot of GPP against LAI",x="LAI", y = "GPP (gC/m2/day)")+
  #geom_smooth(method="lm", se = FALSE, fill = NA,fullrange=TRUE,aes(color=region_name))+
  stat_poly_line()+
  stat_poly_eq(aes(label = paste(after_stat(rr.label), sep = "*\", \"*")))+
  geom_point()+
  theme_bw() + 
  theme(text = element_text(size = 16),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

all_pixels_compare_data_zero %>%
  ggplot(aes(x=cica, 
             y=gpp, 
             colour = region_name))+
  geom_point(alpha = 0.2)+
  labs(title="Plot of GPP against CiCa",x="CiCa", y = "GPP (gC/m2/day)")+
  #geom_smooth(method="lm", se = FALSE, fill = NA,fullrange=TRUE,aes(color=region_name))+
  stat_poly_line()+
  stat_poly_eq(aes(label = paste(after_stat(rr.label), sep = "*\", \"*")))+
  geom_point()+
  theme_bw() + 
  theme(text = element_text(size = 16),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#done
###################new plots and Facet Wrapping##################
all_pixels_compare_data_zero %>% 
  group_by(month,region_name)%>%
  ggplot(aes(x=lai, 
             y=gpp, 
             colour = region_name))+
  geom_point(alpha = 0.2)+
  labs(title="Plot of GPP against LAI",x="LAI", y = "GPP (gC/m2/day)")+
  stat_poly_line()+
  stat_poly_eq(aes(label = paste(after_stat(rr.label), sep = "*\", \"*")))+
  geom_point()+
  theme_bw() + 
  theme(text = element_text(size = 16),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

all_pixels_compare_data_zero %>%
  ggplot(aes(x=cica, 
             y=gpp, 
             colour = region_name))+
  geom_point(alpha = 0.2)+
  labs(title="Plot of GPP against CiCa",x="CiCa", y = "GPP (gC/m2/day)")+
  stat_poly_line()+
  stat_poly_eq(aes(label = paste(after_stat(rr.label), sep = "*\", \"*")))+
  geom_point()+
  theme_bw() + 
  theme(text = element_text(size = 16),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


####new plots#####
ggplot(all_pixels_compare_data_zero, aes(lai, gpp, colour = region_name)) +
  geom_point()+
  labs(title="Plot of Gross Primary Productivity against Leaf Area Index",x="LAI (m2/m2)", y = "GPP (gC/m2/day)") +
  stat_poly_line()+
  stat_poly_eq(aes(label = paste(after_stat(rr.label), sep = "*\", \"*"))) +
  facet_wrap(.~month_f)+
  theme_bw() + 
  theme(text = element_text(size = 16),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggplot(all_pixels_compare_data_zero, aes(cica, gpp, colour = region_name)) +
  geom_point()+
  labs(title="Plot of Gross Primary Productivity against CiCa",x="CiCa", y = "GPP (gC/m2/day)") +
  stat_poly_line()+
  stat_poly_eq(aes(label = paste(after_stat(rr.label), sep = "*\", \"*"))) +
  facet_wrap(.~month_f)+
  theme_bw() + 
  theme(text = element_text(size = 16),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggplot(all_pixels_compare_data_zero, aes(gsdsr, gpp, colour = region_name)) +
  geom_point()+
  labs(title="Plot of Gross Primary Productivity against gs demand:supply",x="GSDSR", y = "GPP (gC/m2/day)") +
  stat_poly_line()+
  stat_poly_eq(aes(label = paste(after_stat(rr.label), sep = "*\", \"*"))) +
  facet_wrap(.~month_f)+
  theme_bw() + 
  theme(text = element_text(size = 16),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggplot(all_pixels_compare_data_zero, aes(apar, gpp, colour = region_name)) +
  geom_point()+
  labs(title="Plot of Gross Primary Productivity against Absorbed Photosynthetically Active Radiation",x="APAR (MJ/m2/day)", y = "GPP (gC/m2/day)") +
  stat_poly_line()+
  stat_poly_eq(aes(label = paste(after_stat(rr.label), sep = "*\", \"*"))) +
  facet_wrap(.~month_f)+
  theme_bw() + 
  theme(text = element_text(size = 16),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggplot(all_pixels_compare_data_zero, aes(et, gpp, colour = region_name)) +
  geom_point()+
  labs(title="Plot of Gross Primary Productivity against Evapotranspiration",x="ET (kgH20/m2/day)", y = "GPP (gC/m2/day)") +
  stat_poly_line()+
  stat_poly_eq(aes(label = paste(after_stat(rr.label), sep = "*\", \"*"))) +
  facet_wrap(.~month_f)+
  theme_bw() + 
  theme(text = element_text(size = 16),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggplot(all_pixels_compare_data_zero, aes(gsdsr, apar, colour = region_name)) +
  geom_point()+
  labs(title="Plot of Absorbed Photosynthetically Active Radiation against gs demand:supply",x="GSDSR", y = "APAR (MJ/m2/day)") +
  stat_poly_line()+
  stat_poly_eq(aes(label = paste(after_stat(rr.label), sep = "*\", \"*"))) +
  facet_wrap(.~month_f)+
  theme_bw() + 
  theme(text = element_text(size = 16),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggplot(all_pixels_compare_data_zero, aes(et, gsdsr, colour = region_name)) +
  geom_point()+
  labs(title="Plot of gs demand:supply against Evapotranspiration",x="ET (kgH20/m2/day)", y = "GSDSR") +
  stat_poly_line()+
  stat_poly_eq(aes(label = paste(after_stat(rr.label), sep = "*\", \"*"))) +
  facet_wrap(.~month_f)+
  theme_bw() + 
  theme(text = element_text(size = 16),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggplot(all_pixels_compare_data_zero, aes(et, apar, colour = region_name)) +
  geom_point()+
  labs(title="Plot of Absorbed Photosynthetically Active Radiation against Evapotranspiration",x="ET (kgH20/m2/day)", y = "APAR (MJ/m2/day)") +
  stat_poly_line()+
  stat_poly_eq(aes(label = paste(after_stat(rr.label), sep = "*\", \"*"))) +
  facet_wrap(.~month_f)+
  theme_bw() + 
  theme(text = element_text(size = 16),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggplot(all_pixels_compare_data_zero[all_pixels_compare_data$region_name=='amazonia',], aes(et, apar)) +
  geom_hex(bins = 20)+
  labs(title="Plot of Absorbed Photosynthetically Active Radiation against Evapotranspiration",x="ET (kgH20/m2/day)", y = "APAR (MJ/m2/day)")+
  scale_fill_viridis_c()+
  facet_wrap(.~month_f)+
  theme_bw() + 
  theme(text = element_text(size = 16),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
