#####################################################################
##########################Compare GPP, LAI and CiCa##################
#####################################################################
library(ncdf4); library(raster); library(dplyr); library(ggplot2);library(ggpubr);library(quantreg);library(ggpp);library(rgdal)

######Read in shape files ###########
amazonia_nw <- shapefile("./amazon_nw.shp");amazonia_sw <- shapefile("./amazon_sw.shp");amazonia_ec <- shapefile("./amazon_ec.shp")
amazonia_bs <- shapefile("./amazon_bs.shp");amazonia_gs <- shapefile("./amazon_gs.shp")

var_names <- c('GPP','LAI','CICA','GSDSR','APAR','ET')
reg_data <- list(amazonia_nw,amazonia_sw,amazonia_ec,amazonia_bs,amazonia_gs)
reg_names <- c('amazon_nw','amazon_sw','amazon_ec','amazon_bs','amazon_gs','amazonia')
names(reg_data) <- reg_names

# function to extract country amazon regions
amazon_reg_to_card_var <- function(s,s2,t) {
  r <- stack('R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/DATA/benchmark/NoRainfor_Amazon_1deg_monthly_2001_updated_2019_compare.nc',varname=t)
  a <- r$X2001.01.01
  cropped <- crop(a, extent(s))
  masked1 <- mask(cropped, s)
  masked2 <- masked1 > -Inf
  poly <- rasterToPolygons(masked2, dissolve=TRUE)
  card_region <- mask(r,poly)
  df_card_region <- as.data.frame(card_region,xy=F)
  df_card_region_cm <- colMeans(df_card_region,na.rm=T)
  df_card_region_cm <- as.data.frame(df_card_region_cm)
  names(df_card_region_cm) <- paste(tolower(t))
  return(df_card_region_cm)
}
#cica_nw <- amazon_reg_to_card_var(reg_data[[1]],reg_names[1],var_names[3])

amazon_to_card_var <- function(s2,t) {
  card_region <- stack('R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/DATA/benchmark/NoRainfor_Amazon_1deg_monthly_2001_updated_2019_compare.nc',varname=t)
  df_card_region <- as.data.frame(card_region,xy=F)
  df_card_region_cm <- colMeans(df_card_region,na.rm=T)
  df_card_region_cm <- as.data.frame(df_card_region_cm)
  names(df_card_region_cm) <- paste(tolower(t))
  
  return(df_card_region_cm)
}
#cica_amazon <- amazon_to_card_var('amazon',var_names[3])

reg_card_var <- function(s,s2,cardamom_variables){
  for (i in cardamom_variables){
    if (i=='GPP') {
      gpp_reg_data <- amazon_reg_to_card_var(s,s2,i)}
    else if (i=='LAI') {
      lai_reg_data <- amazon_reg_to_card_var(s,s2,i)}
    else if (i=='CICA') {
      cica_reg_data <- amazon_reg_to_card_var(s,s2,i)}
    else if (i=='GSDR') {
      gsdr_reg_data <- amazon_reg_to_card_var(s,s2,i)}
    else if (i=='APAR') {
      apar_reg_data <- amazon_reg_to_card_var(s,s2,i)}
    else{
      et_reg_data <- amazon_reg_to_card_var(s,s2,i)
    }
  }
  region_name <- as.factor(rep(s2,nrow(gpp_reg_data)))
  return(cbind(gpp_reg_data,lai_reg_data,cica_reg_data,region_name))
}
#nw_data <- reg_card_var(reg_data[[1]],reg_names[1],var_names)
card_var <- function(s2,cardamom_variables){
  for (i in cardamom_variables){
    if (i=='GPP') {
      gpp_amazon_data <- amazon_to_card_var(s2,i)}
    else if (i=='LAI') {
      lai_amazon_data <- amazon_to_card_var(s2,i)}
    else{
      cica_amazon_data <- amazon_to_card_var(s2,i)
    }
  }
  region_name <- as.factor(rep(s2,nrow(gpp_amazon_data)))
  return(cbind(gpp_amazon_data,lai_amazon_data,cica_amazon_data,region_name))
}
#amazon_data <- card_var(reg_names[6],var_names)

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

all_compare_data <- all_reg_run(reg_data,reg_names,var_names)
names(all_compare_data);head(all_compare_data)

new_all_compare_data <- all_reg_run(reg_data,reg_names,var_names)
names(new_all_compare_data);head(new_all_compare_data);str(new_all_compare_data)

####Time Series Plot######
values = seq(from = as.Date("2001-01-01"), to = as.Date("2019-12-01"), by = 'month')
all_compare_data$date <- as.Date(values)
all_compare_data$month <- format(values,format="%b")
all_compare_data$year <- format(values,format="%Y")
#all_compare_data$mnth <- as.factor(format(values,format="%m"))


all_comparison_data_ts <- ts(all_compare_data, start=c(2001, 1), end=c(2019, 12), frequency=12)
plot.ts(all_comparison_data_ts[,c(1,2,3)],main = "North West", ylab = "")
plot.ts(all_comparison_data_ts[,c(4,5,6)],main = "South West", ylab = "")
plot.ts(all_comparison_data_ts[,c(7,8,9)],main = "East Central", ylab = "")
plot.ts(all_comparison_data_ts[,c(10,11,12)],main = "Brazil Shield", ylab = "")
plot.ts(all_comparison_data_ts[,c(13,14,15)],main = "Guiana Shield", ylab = "")
plot.ts(all_comparison_data_ts[,c(16,17,18)],main = "Amazonia", ylab = "")

#monthly mean
par(mfrow=c(3,1))
boxplot(all_compare_data$lai_amazon_nw~all_compare_data$month,ylab='LAI',xlab='Month')
boxplot(all_compare_data$gpp_amazon_nw~all_compare_data$month,ylab='GPP',xlab='Month',ylim=c(6.25,10))
boxplot(all_compare_data$cica_amazon_nw~all_compare_data$month,ylab='CiCa',xlab='Month',ylim=c(0.6,0.8))
mtext("North West",side=3,line=-2,outer=TRUE, font=2)

boxplot(all_compare_data$lai_amazon_sw~all_compare_data$month,ylab='LAI',xlab='Month')
boxplot(all_compare_data$gpp_amazon_sw~all_compare_data$month,ylab='GPP',xlab='Month',ylim=c(6.25,10))
boxplot(all_compare_data$cica_amazon_sw~all_compare_data$month,ylab='CiCa',xlab='Month',ylim=c(0.6,0.8))
mtext("South West",side=3,line=-2,outer=TRUE, font=2)

boxplot(all_compare_data$lai_amazon_ec~all_compare_data$month,ylab='LAI',xlab='Month')
boxplot(all_compare_data$gpp_amazon_ec~all_compare_data$month,ylab='GPP',xlab='Month',ylim=c(6.25,10))
boxplot(all_compare_data$cica_amazon_ec~all_compare_data$month,ylab='CiCa',xlab='Month',ylim=c(0.6,0.8))
mtext("East Central",side=3,line=-2,outer=TRUE, font=2)

boxplot(all_compare_data$lai_amazon_bs~all_compare_data$month,ylab='LAI',xlab='Month')
boxplot(all_compare_data$gpp_amazon_bs~all_compare_data$month,ylab='GPP',xlab='Month',ylim=c(6.25,10))
boxplot(all_compare_data$cica_amazon_bs~all_compare_data$month,ylab='CiCa',xlab='Month',ylim=c(0.6,0.8))
mtext("Brazil Shield",side=3,line=-2,outer=TRUE, font=2)

boxplot(all_compare_data$lai_amazon_gs~all_compare_data$month,ylab='LAI',xlab='Month')
boxplot(all_compare_data$gpp_amazon_gs~all_compare_data$month,ylab='GPP',xlab='Month',ylim=c(6.25,10))
boxplot(all_compare_data$cica_amazon_gs~all_compare_data$month,ylab='CiCa',xlab='Month',ylim=c(0.6,0.8))
mtext("Guiana Shield",side=3,line=-2,outer=TRUE, font=2)

boxplot(all_compare_data$lai_amazonia~all_compare_data$month,ylab='LAI',xlab='Month')
boxplot(all_compare_data$gpp_amazonia~all_compare_data$month,ylab='GPP',xlab='Month',ylim=c(6.25,10))
boxplot(all_compare_data$cica_amazonia~all_compare_data$month,ylab='CiCa',xlab='Month',ylim=c(0.6,0.8))
mtext("Amazonia",side=3,line=-2,outer=TRUE, font=2)


gpplai_amazon_nw<-summary(lm(gpp_amazon_nw ~ lai_amazon_nw, data = all_compare_data))
rsq_gpplai_amazon_nw<-paste0("R.sp = ", round(gpplai_amazon_nw$r.squared,2))
gppcica_amazon_nw<-summary(lm(gpp_amazon_nw ~ cica_amazon_nw, data = all_compare_data))
rsq_gppcica_amazon_nw<-paste0("R.sp = ", round(gppcica_amazon_nw$r.squared,2))
      
gpplai_amazon_sw<-summary(lm(gpp_amazon_sw ~ lai_amazon_sw, data = all_compare_data))
rsq_gpplai_amazon_sw<-paste0("R.sp = ", round(gpplai_amazon_sw$r.squared,2))
gppcica_amazon_sw<-summary(lm(gpp_amazon_sw ~ cica_amazon_sw, data = all_compare_data))
rsq_gppcica_amazon_sw<-paste0("R.sp = ", round(gppcica_amazon_sw$r.squared,2))

gpplai_amazon_ec<-summary(lm(gpp_amazon_ec ~ lai_amazon_ec, data = all_compare_data))
rsq_gpplai_amazon_ec<-paste0("R.sp = ", round(gpplai_amazon_ec$r.squared,2))
gppcica_amazon_ec<-summary(lm(gpp_amazon_ec ~ cica_amazon_ec, data = all_compare_data))
rsq_gppcica_amazon_ec<-paste0("R.sp = ", round(gppcica_amazon_ec$r.squared,2))
    
gpplai_amazon_bs<-summary(lm(gpp_amazon_bs ~ lai_amazon_bs, data = all_compare_data))
rsq_gpplai_amazon_bs<-paste0("R.sp = ", round(gpplai_amazon_bs$r.squared,2))
gppcica_amazon_bs<-summary(lm(gpp_amazon_bs ~ cica_amazon_bs, data = all_compare_data))
rsq_gppcica_amazon_bs<-paste0("R.sp = ", round(gppcica_amazon_bs$r.squared,2))
    
gpplai_amazon_gs<-summary(lm(gpp_amazon_gs ~ lai_amazon_gs, data = all_compare_data))
rsq_gpplai_amazon_gs<-paste0("R.sp = ", round(gpplai_amazon_gs$r.squared,2))
gppcica_amazon_gs<-summary(lm(gpp_amazon_gs ~ cica_amazon_gs, data = all_compare_data))
rsq_gppcica_amazon_gs<-paste0("R.sp = ", round(gppcica_amazon_gs$r.squared,2))
    
gpplai_amazonia<-summary(lm(gpp_amazonia ~ lai_amazonia, data = all_compare_data))
rsq_gpplai_amazonia<-paste0("R.sp = ", round(gpplai_amazonia$r.squared,2))
gppcica_amazonia<-summary(lm(gpp_amazonia ~ cica_amazonia, data = all_compare_data))
rsq_gppcica_amazonia<-paste0("R.sp = ", round(gppcica_amazonia$r.squared,2))
    
#rsq_amazon_rltnshps<-rsq_relationships(all_compare_data,reg_names)

# gpplai<-summary(lm(gpp_mean ~ lai_mean, data = all_compare_data))
# rsq_gpplai<-paste0("R.sp = ", round(gpplai$r.squared,2))
# gppcica<-summary(lm(gpp_mean ~ cica_mean, data = all_compare_data))
# rsq_gppcica<-paste0("R.sp = ", round(gppcica$r.squared,2))

all_compare_data<-as.data.frame(all_compare_data)
par(mfrow=c(1,2))
plot(all_compare_data$gpp_amazon_nw~all_compare_data$lai_amazon_nw,xlab='LAI', ylab='GPP', main=rsq_gpplai_amazon_nw)
plot(all_compare_data$gpp_amazon_nw~all_compare_data$cica_amazon_nw, main=rsq_gppcica_amazon_nw,ylab='GPP',xlab='CiCa')
mtext("North West",side=3,line=-3,outer=TRUE, font=2,cex=1.5)

plot(all_compare_data$gpp_amazon_sw~all_compare_data$lai_amazon_sw, main=rsq_gpplai_amazon_sw,ylab='GPP',xlab='LAI')
plot(all_compare_data$gpp_amazon_sw~all_compare_data$cica_amazon_sw, main=rsq_gppcica_amazon_sw,ylab='GPP',xlab='CiCa')
mtext("South West",side=3,line=-3,outer=TRUE, font=2,cex=1.5)

plot(all_compare_data$gpp_amazon_ec~all_compare_data$lai_amazon_ec, main=rsq_gpplai_amazon_ec,ylab='GPP',xlab='LAI')
plot(all_compare_data$gpp_amazon_ec~all_compare_data$cica_amazon_ec, main=rsq_gppcica_amazon_ec,ylab='GPP',xlab='CiCa')
mtext("East Central",side=3,line=-3,outer=TRUE, font=2,cex=1.5)

plot(all_compare_data$gpp_amazon_bs~all_compare_data$lai_amazon_bs, main=rsq_gpplai_amazon_bs,ylab='GPP',xlab='LAI')
plot(all_compare_data$gpp_amazon_bs~all_compare_data$cica_amazon_bs, main=rsq_gppcica_amazon_bs,ylab='GPP',xlab='CiCa')
mtext("Brazil Shield",side=3,line=-3,outer=TRUE, font=2,cex=1.5)

plot(all_compare_data$gpp_amazon_gs~all_compare_data$lai_amazon_gs, main=rsq_gpplai_amazon_gs,ylab='GPP',xlab='LAI')
plot(all_compare_data$gpp_amazon_gs~all_compare_data$cica_amazon_gs, main=rsq_gppcica_amazon_gs,ylab='GPP',xlab='CiCa')
mtext("Guiana Shield",side=3,line=-3,outer=TRUE, font=2,cex=1.5)

plot(all_compare_data$gpp_amazonia~all_compare_data$lai_amazonia, main=rsq_gpplai_amazonia,ylab='GPP',xlab='LAI')
plot(all_compare_data$gpp_amazonia~all_compare_data$cica_amazonia, main=rsq_gppcica_amazonia,ylab='GPP',xlab='CiCa')
mtext("Amazonia",side=3,line=-3,outer=TRUE, font=2,cex=1.5)


#####multiple plots
formula <- y ~ x

new_all_compare_data %>%
  ggplot(aes(x=lai, 
             y=gpp, 
             shape = region_name))+
  geom_point(alpha = 0.2)+
  geom_smooth(method="lm", se = FALSE, fill = NA,fullrange=TRUE,aes(color=region_name))+
  theme_bw() + 
  theme(text = element_text(size = 16),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#ggsave("scatterplot_with_multiple_groups_ggplot2.png")

new_all_compare_data %>%
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

new_all_compare_data %>%
  ggplot(aes(x=lai, 
             y=gpp, 
             colour = region_name))+
  labs(title="Plot of GPP against LAI",x="LAI", y = "GPP (gC/m2/day)")+
  #geom_smooth(method="lm", se = FALSE, fill = NA,fullrange=TRUE,aes(color=region_name))+
  stat_poly_line()+
  stat_poly_eq(aes(label = paste(after_stat(rr.label), sep = "*\", \"*")))+
  geom_point(alpha = 0.2)+
  theme_bw() + 
  theme(text = element_text(size = 16),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

new_all_compare_data %>%
  ggplot(aes(x=cica, 
             y=gpp, 
             colour = region_name))+
  labs(title="Plot of GPP against CiCa",x="CiCa", y = "GPP (gC/m2/day)")+
  #geom_smooth(method="lm", se = FALSE, fill = NA,fullrange=TRUE,aes(color=region_name))+
  stat_poly_line()+
  stat_poly_eq(aes(label = paste(after_stat(rr.label), sep = "*\", \"*")))+
  geom_point(alpha = 0.2)+
  theme_bw() + 
  theme(text = element_text(size = 16),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
