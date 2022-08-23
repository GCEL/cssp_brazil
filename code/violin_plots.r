#####################################################################
##########################Compare GPP, LAI and CiCa##################
#####################################################################
library(ncdf4); library(raster); library(ggplot2)

#####extract variables from different models########
amazon_biomass_gCm2_01_16_monthly
amazon_woodprod_gCm2d_01_16_monthly 
amazon_biommort_gCm2d_01_16_monthly 

amazon_woodprod_gCm2d_01_16_monthly_2 <- mask(amazon_woodprod_gCm2d_01_16_monthly, amazon_biomass_gCm2_01_16_monthly)
amazon_biommort_gCm2d_01_16_monthly_2 <- mask(amazon_biommort_gCm2d_01_16_monthly, amazon_biomass_gCm2_01_16_monthly)

cardamom_wb<- stack('R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/DATA/benchmark/Amazon_1deg_monthly_2001_2019_cwood_2010.nc',varname='WOOD')
cardamom_wp<- stack('R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/DATA/benchmark/Amazon_1deg_monthly_2001_2019_cwood_2010.nc',varname='NPP_wood_flx')
cardamom_wm<- stack('R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/DATA/benchmark/Amazon_1deg_monthly_2001_2019_cwood_2010.nc',varname='OUTPUT_wood_flx')

inland_wb<- stack('R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/MODELS/INLAND/INLAND_monthly_1x1_Amazon_2001_2019.nc',varname='WOOD')
inland_wp<- stack('R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/MODELS/INLAND/INLAND_monthly_1x1_Amazon_2001_2019.nc',varname='NPP_wood_flx')
inland_wm<- stack('R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/MODELS/INLAND/INLAND_monthly_1x1_Amazon_2001_2019.nc',varname='OUTPUT_wood_flx')

jules_wb<- stack('R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/MODELS/JULES/JULES_monthly_1x1_Amazon_2001_2019.nc',varname='WOOD')
jules_wp<- stack('R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/MODELS/JULES/JULES_monthly_1x1_Amazon_2001_2019.nc',varname='NPP_wood_flx')
jules_wm<- stack('R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/MODELS/JULES/JULES_monthly_1x1_Amazon_2001_2019.nc',varname='OUTPUT_wood_flx')


model_var_names <- c('woody biomass','woody productivity','woody mortality')
model_data <- list(inland_wb,inland_wp,inland_wm,jules_wb,jules_wp,jules_wm,cardamom_wb,cardamom_wp,cardamom_wm,amazon_biomass_gCm2_01_16_monthly,amazon_woodprod_gCm2d_01_16_monthly_2,amazon_biommort_gCm2d_01_16_monthly_2)
model_names <- c('inland','jules','cardamom_dalec','rainfor')

model_raster_mean_over_time_fun <- function(mr){
  mr_01_16 <- mr[[which(getZ(mr) >= as.Date("2001-01-01") & getZ(mr) <= as.Date("2016-12-01"))]]
  mr_01_16_mean <- stackApply(mr_01_16, indices =  rep(1,nlayers(mr_01_16)), fun = "mean")
}
#extract raster images of two periods of interest
inland_wb_01_16 <- inland_wb[[which(getZ(inland_wb) >= as.Date("2001-01-01") & getZ(inland_wb) <= as.Date("2016-12-01"))]]

#calculate mean over time
inland_wb_01_16_mean <- stackApply(inland_wb_01_16, indices =  rep(1,nlayers(inland_wb_01_16)), fun = "mean")
plot(inland_wb_01_16_mean)

inland_wb_01_16_mean<-model_raster_mean_over_time_fun(inland_wb)

model_amazon_crop_fun <- function(p,mr) {
  model_extent <- extent(mr) #use south america extent from INLAND SA
  r <- raster(model_extent)
  res(r)<- res(mr) #use resolution to mirror other masks in ilamb/ could be 1x1
  values(r) <- 1
  crs(r) <- "+proj=longlat +datum=WGS84 +no_defs"
  r_amazon <- mask(r, p, updatevalue=0)#create mask from polygon of RAINFOR amazon region
  r_amazon[r_amazon < 1] <- NA
  r_amazon <- r_amazon> -Inf
  r_amazon_pol <- rasterToPolygons(r_amazon, dissolve = TRUE)
  mr_crop <- crop(mr, extent(r_amazon_pol))
  mr_amazonia <- mask(mr_crop, r_amazon_pol)
}


mr_mot_amazonia_df_fun <- function(mr,modelname,variablename) {
  mr_01_16 <- mr[[which(getZ(mr) >= as.Date("2001-01-01") & getZ(mr) <= as.Date("2016-12-01"))]]
  mr_01_16_mean <- stackApply(mr_01_16, indices =  rep(1,nlayers(mr_01_16)), fun = "mean")
  model_extent <- extent(mr_01_16_mean) #use south america extent from INLAND SA
  r <- raster(model_extent)
  res(r)<- res(mr_01_16_mean) #use resolution to mirror other masks in ilamb/ could be 1x1
  values(r) <- 1
  crs(r) <- "+proj=longlat +datum=WGS84 +no_defs"
  r_amazon <- mask(r, amazonia_poly, updatevalue=0)#create mask from polygon of RAINFOR amazon region
  r_amazon[r_amazon < 1] <- NA
  r_amazon <- r_amazon> -Inf
  r_amazon_pol <- rasterToPolygons(r_amazon, dissolve = TRUE)
  mr_01_16_mean_crop <- crop(mr_01_16_mean, extent(r_amazon_pol))
  mr_01_16_mean_crop_amazonia <- mask(mr_01_16_mean_crop, r_amazon_pol)
  df_mr_01_16_mean_crop_amazonia <- as.data.frame(mr_01_16_mean_crop_amazonia,na.rm=T)
  names(df_mr_01_16_mean_crop_amazonia)<-paste(variablename)
  rownames(df_mr_01_16_mean_crop_amazonia) <- 1:nrow(df_mr_01_16_mean_crop_amazonia)
  df_mr_01_16_mean_crop_amazonia$variablename  <- df_mr_01_16_mean_crop_amazonia$index_1
  df_mr_01_16_mean_crop_amazonia$model_name <- as.factor(rep(modelname,nrow(df_mr_01_16_mean_crop_amazonia)))
  return(df_mr_01_16_mean_crop_amazonia)
}

benchmark_amazonia_df_fun <- function(mr,modelname,variablename) {
  mr_01_16_mean <- stackApply(mr, indices =  rep(1,nlayers(mr)), fun = "mean")
  df_mr_01_16_mean_crop_amazonia <- as.data.frame(mr_01_16_mean,xy=F,na.rm=T)
  names(df_mr_01_16_mean_crop_amazonia)<-paste(variablename)
  rownames(df_mr_01_16_mean_crop_amazonia) <- 1:nrow(df_mr_01_16_mean_crop_amazonia)
  df_mr_01_16_mean_crop_amazonia$variablename  <- df_mr_01_16_mean_crop_amazonia$index_1
  df_mr_01_16_mean_crop_amazonia$model_name <- as.factor(rep(modelname,nrow(df_mr_01_16_mean_crop_amazonia)))
  return(df_mr_01_16_mean_crop_amazonia)
}

df_inland_wb_01_16_mean_amazonia <- mr_mot_amazonia_df_fun(inland_model_data[[1]],model_names[1],model_var_names[1])
df_inland_wp_01_16_mean_amazonia <- mr_mot_amazonia_df_fun(inland_model_data[[2]],model_names[1],model_var_names[2])
df_inland_wm_01_16_mean_amazonia <- mr_mot_amazonia_df_fun(inland_model_data[[3]],model_names[1],model_var_names[3])

df_inland_wb_01_16_mean_amazonia<-df_inland_wb_01_16_mean_amazonia[,c(2,1,3,5)]

df_inland_01_16_mean_amazonia<-cbind(df_inland_wb_01_16_mean_amazonia,df_inland_wp_01_16_mean_amazonia,df_inland_wm_01_16_mean_amazonia)
df_inland_01_16_mean_amazonia<-df_inland_01_16_mean_amazonia[,c(2,1,3,5)]

df_benchmark_wb_01_16_mean_amazonia <- benchmark_amazonia_df_fun(model_data[[10]],model_names[4],model_var_names[1])
df_benchmark_wp_01_16_mean_amazonia <- benchmark_amazonia_df_fun(model_data[[11]],model_names[4],model_var_names[2])
df_benchmark_wm_01_16_mean_amazonia <- benchmark_amazonia_df_fun(model_data[[12]],model_names[4],model_var_names[3])

df_benchmark_wb_01_16_mean_amazonia<-df_benchmark_wb_01_16_mean_amazonia[,c(2,1,3,5)]

df_benchmark_01_16_mean_amazonia<-cbind(df_benchmark_wb_01_16_mean_amazonia,df_benchmark_wp_01_16_mean_amazonia,df_benchmark_wm_01_16_mean_amazonia)
df_benchmark_01_16_mean_amazonia<-df_benchmark_01_16_mean_amazonia[,c(2,1,3,5)]

all_model_vars_fun <- function(mr,modelname,variablename){
  for (i in variablename){
    if (i=='woody biomass') {
      wb_data <- mr_mot_amazonia_df_fun(mr[[1]],modelname,i)}
    else if (i=='woody productivity') {
      wp_data <- mr_mot_amazonia_df_fun(mr[[2]],modelname,i)}
    else{
      wm_data <- mr_mot_amazonia_df_fun(mr[[3]],modelname,i)
    }
  }
  model_all_var <- cbind(wb_data,wp_data,wm_data)
  return(model_all_var[,c(2,1,3,5)])
}

benchmark_vars_fun <- function(mr,modelname,variablename){
  for (i in variablename){
    if (i=='woody biomass') {
      wb_data <- benchmark_amazonia_df_fun(mr[[1]],modelname,i)}
    else if (i=='woody productivity') {
      wp_data <- benchmark_amazonia_df_fun(mr[[2]],modelname,i)}
    else{
      wm_data <- benchmark_amazonia_df_fun(mr[[3]],modelname,i)
    }
  }
  model_all_var <- cbind(wb_data,wp_data,wm_data)
  return(model_all_var[,c(2,1,3,5)])
}

df_inland_01_16_mean_amazonia <- all_model_vars_fun(model_data[1:3],model_names[1],model_var_names)
df_jules_01_16_mean_amazonia <- all_model_vars_fun(model_data[4:6],model_names[2],model_var_names)

rbind(df_inland_01_16_mean_amazonia,df_jules_01_16_mean_amazonia)

df_benchmark_01_16_mean_amazonia <- benchmark_vars_fun(model_data[10:12],model_names[4],model_var_names)


all_model_fun <- function(mr,modelname,variablename) {
  for (i in modelname) {
    if(i=='inland'){
      inland_data <- all_model_vars_fun(mr[1:3],i,variablename)
    }
    else if (i=='jules') {
      jules_data <- all_model_vars_fun(mr[4:6],i,variablename)
    }
    else if (i=='cardamom_dalec') {
      cardamom_data <- all_model_vars_fun(mr[7:9],i,variablename)
    }
    else{
      rainfor_data <- benchmark_vars_fun(mr[10:12],i,variablename)
    }
  }
  return(rbind(inland_data,jules_data,cardamom_data,rainfor_data))
}

df_models_01_16_mean_amazonia <- all_model_fun(model_data,model_names,model_var_names)
summary(df_models_01_16_mean_amazonia)
str(df_models_01_16_mean_amazonia)

names(df_models_01_16_mean_amazonia) <- c('model_name','woody_biomass','woody_productivity','woody_mortality')
df_models_01_16_mean_amazonia$woody_biomass[df_models_01_16_mean_amazonia$woody_biomass < 0] <- 0
df_models_01_16_mean_amazonia$woody_productivity[df_models_01_16_mean_amazonia$woody_productivity < 0] <- 0
df_models_01_16_mean_amazonia$woody_mortality[df_models_01_16_mean_amazonia$woody_mortality < 0] <- 0

df_models_01_16_mean_amazonia$model_name <- factor(df_models_01_16_mean_amazonia$model_name, levels = c('rainfor','cardamom_dalec','inland', 'jules'))

# Basic violin plot
ggplot(df_models_01_16_mean_amazonia, aes(x=dose, y=len)) + 
  geom_violin()

# Rotate the violin plot
p + coord_flip()
# Set trim argument to FALSE
ggplot(df_models_01_16_mean_amazonia, aes(x=model_name, y=woody_biomass, fill=model_name)) + 
  geom_violin(trim=FALSE)+ 
  labs(title="Plot of Woody Biomass by Model",x="Model", y = "Woody Biomass (gC/m2)")+
  #scale_x_discrete(limits=c('rainfor','cardamom_dalec','inland', 'jules'))+
  geom_boxplot(width=0.05)+
  scale_fill_brewer(palette="Blues") +
  guides(fill=guide_legend('Model Name')) +
  #scale_color_manual(values=c( "pink", "orange","purple","black")) +
  theme_bw() + 
  theme(text = element_text(size = 16),panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ 
  coord_flip()

ggplot(df_models_01_16_mean_amazonia, aes(x=model_name, y=woody_productivity, fill=model_name)) + 
  geom_violin(trim=FALSE)+ 
  labs(title="Plot of Woody Productivity by Model",x="Model", y = "Woody Productivity (gC/m2/day)")+
  #scale_x_discrete(limits=c('rainfor','cardamom_dalec','inland', 'jules'))+
  stat_summary(fun.y=median, geom="point", size=3, color="green")+
  geom_boxplot(width=0.05)+
  scale_fill_brewer(palette="Blues") +
  guides(fill=guide_legend('Model Name')) +
  #scale_color_manual(values=c( "pink", "orange","purple","black")) +
  theme_bw() + 
  theme(text = element_text(size = 16),panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ 
  coord_flip()

ggplot(df_models_01_16_mean_amazonia, aes(x=model_name, y=woody_mortality, fill=model_name)) + 
  geom_violin(trim=FALSE)+ 
  labs(title="Plot of Woody Mortality by Model",x="Model", y = "Woody Mortality (gC/m2/day)")+
  #scale_x_discrete(limits=c('rainfor','cardamom_dalec','inland', 'jules'))+
  stat_summary(fun.y=median, geom="point", size=3, color="red")+
  geom_boxplot(width=0.05)+
  scale_fill_brewer(palette="Blues") +
  guides(fill=guide_legend('Model Name')) +
  #scale_color_manual(values=c( "pink", "orange","purple","black")) +
  theme_bw() + 
  theme(text = element_text(size = 16),panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ 
  coord_flip()


