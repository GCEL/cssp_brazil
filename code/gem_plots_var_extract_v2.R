##########Masks for Creating Ecoregions##################
amazon_region_gem <- extent(biomass_C_amazon)
r_amazon_region_gem <- raster(amazon_region_gem,res=1)
values(r_amazon_region_gem) <- 0
r_amazon_region_gem <- rasterize(biomass_amazon_mask_pol,r_amazon_region_gem)
plot(r_amazon_region_gem)

allpahuayo_poly <- shapefile("./data/Allpahuayo.shp")
caxiuana_poly <- shapefile("./data/Caxiuana.shp")
kenia_poly <- shapefile("./data/Kenia.shp")
tambopata_poly <- shapefile("./data/Tambopata.shp")
tanguro_poly <- shapefile("./data/Tanguro.shp")

# plot(biomass_C_amazon);plot(amazon_nw_poly,add=T);plot(amazon_sw_poly,add=T);plot(amazon_ec_poly,add=T);plot(amazon_bs_poly,add=T);plot(amazon_gs_poly,add=T)

plot(biomass_amazon_C_ifl);plot(amazon_nw_poly,add=T);plot(amazon_sw_poly,add=T);plot(amazon_ec_poly,add=T);plot(amazon_bs_poly,add=T);plot(amazon_gs_poly,add=T);plot(allpahuayo_poly,add=T);plot(caxiuana_poly,add=T);plot(kenia_poly,add=T);plot(tambopata_poly,add=T);plot(tanguro_poly,add=T)


allpahuayo_pol_mask <- country_to_amazon_crop_fun(biomass_C_amazon,allpahuayo_poly)
caxiuana_pol_mask <- country_to_amazon_crop_fun(biomass_C_amazon,caxiuana_poly)
kenia_pol_mask <- country_to_amazon_crop_fun(biomass_C_amazon,kenia_poly)
tambopata_pol_mask <- country_to_amazon_crop_fun(biomass_C_amazon,tambopata_poly)
tanguro_pol_mask <- country_to_amazon_crop_fun(biomass_C_amazon,tanguro_poly)

r_allpahuayo <- mask(r_amazon_region_gem, allpahuayo_pol_mask, updatevalue=0)
r_caxiuana <- mask(r_allpahuayo, caxiuana_pol_mask, inverse=T, updatevalue=2)
r_kenia <- mask(r_caxiuana, tambopata_pol_mask, inverse=T, updatevalue=3)
r_tambopata <- mask(r_kenia, kenia_pol_mask, inverse=T, updatevalue=4)
r_amazonia_gem <- mask(r_tambopata, tanguro_pol_mask, inverse=T, updatevalue=5)
# plot(r_amazonia_gem)
# plot(flip(r_amazonia_gem,2))
# plot(t(r_amazonia_gem))
# 
# plot(allpahuayo_poly,add=T);plot(caxiuana_poly,add=T);plot(kenia_poly,add=T);plot(tambopata_poly,add=T);plot(tanguro_poly,add=T)
# plot(amazonia_extent,add=T)

amazonia_gem_matrix <- t(as.matrix(flip(r_amazonia_gem,2)))

amazonia_gem_matrices <- array(numeric(),c(37,32,3)) 
amazonia_gem_matrices[,,1]<-amazonia_gem_matrix
amazonia_gem_matrices[,,2]<-amazonia_gem_matrix
amazonia_gem_matrices[,,3]<-amazonia_gem_matrix

gem_plots_data <- read.csv("./data/gem_plot_data.csv", header=TRUE)

# test_gem <- amazonia_gem_matrices==1
# dim(test_gem)
gem_region <- c('allpahuayo','caxiuana','kenia','tambopata','tanguro')
gem_var <- c('gpp','cue','npp','mrt_wood','npp_wood','wood')
gem_prefix <- "M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/amazonia_"
gem_suffix_1 <- "_nomngt/infofile.RData"
gem_midfix <- "_nomngt/RESULTS_PROCESSED/amazonia_"
gem_suffix_2 <- "_nomngt_stock_flux.RData"

for (i in new_mod_var) {
  new_gem_plots_df <- data.frame(matrix(nrow=5,ncol=ncol(gem_plots_data)))
  colnames(new_gem_plots_df) <- colnames(gem_plots_data)
  
  if (i =='esa_cci_agb') {
    new_gem_plots_df$model_plot_name <- rep('C',nrow(new_gem_plots_df))
  }
  else if (i =='rainfor_biomass_annual') {
    new_gem_plots_df$model_plot_name <- rep('R0',nrow(new_gem_plots_df))
  }
  else if (i =='rainfor_biomass_productivity_2005') {
    new_gem_plots_df$model_plot_name <- rep('R1',nrow(new_gem_plots_df))
  }
  else {
    new_gem_plots_df$model_plot_name <- rep('R2',nrow(new_gem_plots_df))
  }
  new_gem_plots_df$gem_region <- gem_region
  
for (e in gem_plot_names) {
  # Load the processed site file
  load(paste(gem_prefix,i,gem_suffix_1,sep=""))
  load(paste(gem_prefix,i,gem_midfix,i,gem_suffix_2,sep=""))
  
  if (e =='allpahuayo') {
    gem_plot_mask <- amazonia_gem_matrices==1
  }
  else if (e =='caxiuana') {
    gem_plot_mask <- amazonia_gem_matrices==2
  }
  else if (e =='kenia') {
    gem_plot_mask <- amazonia_gem_matrices==3
  }
  else if (e =='tambopata') {
    gem_plot_mask <- amazonia_gem_matrices==4
  }
  else {
    gem_plot_mask <- amazonia_gem_matrices==5
  }
dp = 2
quantiles_wanted = c(1,4,7)

#Masking phase
# NATURAL FLUXES
grid_output$mean_gpp_gCm2day[,,quantiles_wanted][!gem_plot_mask]<-NA
grid_output$mean_cue[,,quantiles_wanted][!gem_plot_mask]<-NA
grid_output$mean_npp_gCm2day[,,quantiles_wanted][!gem_plot_mask]<-NA
grid_output$mean_alloc_wood_gCm2day[,,quantiles_wanted][!gem_plot_mask]<-NA
# STOCKS
grid_output$mean_wood_gCm2[,,quantiles_wanted][!gem_plot_mask]<-NA
# MEAN TRANSIENT TIMES
grid_output$MTT_wood_years[,,quantiles_wanted][!gem_plot_mask]<-NA

# Extract or calculate required derived values
# NOTE: unit conversion from gC/m2/day to MgC/ha/yr
# NATURAL FLUXES
gpp_MgChayr = format(round(apply(grid_output$mean_gpp_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
cue = format(round(apply(grid_output$mean_cue[,,quantiles_wanted],3,mean, na.rm=TRUE), digits = 2), nsmall = dp)
npp_MgChayr = format(round(apply(grid_output$mean_npp_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
npp_wood_MgChayr = format(round(apply(grid_output$mean_alloc_wood_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * (365.25/100), digits = dp), nsmall = dp)
# STOCKS
wood_MgCha = format(round(apply(grid_output$mean_wood_gCm2[,,quantiles_wanted],3,mean, na.rm=TRUE) / 100, digits = dp), nsmall = dp)
# MEAN TRANSIENT TIMES
mrt_wood_years = format(round(apply(grid_output$MTT_wood_years[,,quantiles_wanted],3,mean, na.rm=TRUE), digits = dp), nsmall = dp)

#collate all variables
all_var_means <- as.numeric(c(gpp_MgChayr[2],cue[2],npp_MgChayr[2],mtt_wood_years[2],npp_wood_MgChayr[2],wood_MgCha[2],
                              gpp_MgChayr[1],cue[1],npp_MgChayr[1],mtt_wood_years[1],npp_wood_MgChayr[1],wood_MgCha[1],
                              gpp_MgChayr[3],cue[3],npp_MgChayr[3],mtt_wood_years[3],npp_wood_MgChayr[3],wood_MgCha[3]))

#add to dataframe
new_gem_plots_df[new_gem_plots_df$gem_region==e,c(3:20)]<-all_var_means
}
  assign(paste("new_gem_plots_df_",i,sep=""),new_gem_plots_df)
}
gem_plots_model_data <- rbind(gem_plots_data,new_gem_plots_df_esa_cci_agb,new_gem_plots_df_rainfor_biomass_annual,new_gem_plots_df_rainfor_biomass_productivity_2005,new_gem_plots_df_rainfor_biomass_annual_productivity)
gem_plots_model_data$model_plot_name <- as.factor(gem_plots_model_data$model_plot_name)
gem_plots_model_data$gem_region <- as.factor(gem_plots_model_data$gem_region)
str(gem_plots_model_data)

for (v in gem_var) {
  if (v =='gpp') {
    plot_data <- gem_plots_model_data[,c(1,2,3,9,15)]
    names(plot_data) <- c('model_plot_name','gem_region','variable','variable_min','variable_max')
    y_axis<-bquote("GPP ("~ Mg~C~ ha^-1~year^-1~")")
  }
  else if (v =='cue') {
    plot_data <- gem_plots_model_data[,c(1,2,4,10,16)]
    names(plot_data) <- c('model_plot_name','gem_region','variable','variable_min','variable_max')
    y_axis<-bquote("Carbon Use Efficiency")
  }
  else if (v =='npp') {
    plot_data <- gem_plots_model_data[,c(1,2,5,11,17)]
    names(plot_data) <- c('model_plot_name','gem_region','variable','variable_min','variable_max')
    y_axis<-bquote("NPP ("~ Mg~C~ ha^-1~year^-1~")")
  }
  else if (v =='mrt_wood') {
    plot_data <- gem_plots_model_data[,c(1,2,6,12,18)]
    names(plot_data) <- c('model_plot_name','gem_region','variable','variable_min','variable_max')
    y_axis<-bquote("Total Coarse Wood Carbon Residence Time (years)")
  }
  else if (v =='npp_wood') {
    plot_data <- gem_plots_model_data[,c(1,2,7,13,19)]
    names(plot_data) <- c('model_plot_name','gem_region','variable','variable_min','variable_max')
    y_axis<-bquote("Total Coarse Wood Carbon NPP ("~ Mg~C~ ha^-1~year^-1~")")
  }
  else {
    plot_data <- gem_plots_model_data[,c(1,2,8,14,20)]
    names(plot_data) <- c('model_plot_name','gem_region','variable','variable_min','variable_max')
    y_axis<-bquote("Total Coarse Wood Carbon Biomass ("~Mg~C~ ha^-1~")")
  }
  plot_data
  y_axis
}
for (i in gem_region) {
    if (i =='allpahuayo') {
      p1 <- ggplot(data=plot_data[plot_data$gem_region==i,],aes(x = model_plot_name, y = variable))+
        geom_point(data=plot_data[plot_data$gem_region==i,],position=position_dodge(0.1)) +
        geom_errorbar(data=plot_data[plot_data$gem_region==i,],aes(ymin=variable_min, ymax=variable_max), width=.1, position=position_dodge(0.1)) +
        labs(title= "", x="Plot/Model Variant",y=y_axis) +
        theme_test()
    }
    else if (i =='caxiuana') {
      p2 <- ggplot(data=plot_data[plot_data$gem_region==i,],aes(x = model_plot_name, y = variable))+
        geom_point(data=plot_data[plot_data$gem_region==i,],position=position_dodge(0.1)) +
        geom_errorbar(data=plot_data[plot_data$gem_region==i,],aes(ymin=variable_min, ymax=variable_max), width=.1, position=position_dodge(0.1)) +
        labs(title= "", x="Plot/Model Variant",y=y_axis) +
        theme_test()
    }
    else if (i =='kenia') {
      p3 <- ggplot(data=plot_data[plot_data$gem_region==i,],aes(x = model_plot_name, y = variable))+
        geom_point(data=plot_data[plot_data$gem_region==i,],position=position_dodge(0.1)) +
        geom_errorbar(data=plot_data[plot_data$gem_region==i,],aes(ymin=variable_min, ymax=variable_max), width=.1, position=position_dodge(0.1)) +
        labs(title= "", x="Plot/Model Variant",y=y_axis) +
        theme_test()
    }
    else if (i =='tambopata') {
      p4 <- ggplot(data=plot_data[plot_data$gem_region==i,],aes(x = model_plot_name, y = variable))+
        geom_point(data=plot_data[plot_data$gem_region==i,],position=position_dodge(0.1)) +
        geom_errorbar(data=plot_data[plot_data$gem_region==i,],aes(ymin=variable_min, ymax=variable_max), width=.1, position=position_dodge(0.1)) +
        labs(title= "", x="Plot/Model Variant",y=y_axis) +
        theme_test()
    }
    else {
      p5 <- ggplot(data=plot_data[plot_data$gem_region==i,],aes(x = model_plot_name, y = variable))+
        geom_point(data=plot_data[plot_data$gem_region==i,],position=position_dodge(0.1)) +
        geom_errorbar(data=plot_data[plot_data$gem_region==i,],aes(ymin=variable_min, ymax=variable_max), width=.1, position=position_dodge(0.1)) +
        labs(title= "", x="Plot/Model Variant",y=y_axis) +
        theme_test()
    }
    print(ggarrange(p1, p2, p3, p4, p5, nrow = 2, ncol = 3))
  }


extract_plot_gem_var_df <- function (gem_data,gem_variables,gem_plots) {
  for (v in gem_variables) {
    if (v =='gpp') {
      plot_data <- gem_data[,c(1,2,3,9,15)]
      names(plot_data) <- c('model_plot_name','gem_region','variable','variable_min','variable_max')
      y_axis<-bquote("GPP ("~ Mg~C~ ha^-1~year^-1~")")
    }
    else if (v =='cue') {
      plot_data <- gem_data[,c(1,2,4,10,16)]
      names(plot_data) <- c('model_plot_name','gem_region','variable','variable_min','variable_max')
      y_axis<-bquote("Carbon Use Efficiency")
    }
    else if (v =='npp') {
      plot_data <- gem_data[,c(1,2,5,11,17)]
      names(plot_data) <- c('model_plot_name','gem_region','variable','variable_min','variable_max')
      y_axis<-bquote("NPP ("~ Mg~C~ ha^-1~year^-1~")")
    }
    else if (v =='mrt_wood') {
      plot_data <- gem_data[,c(1,2,6,12,18)]
      names(plot_data) <- c('model_plot_name','gem_region','variable','variable_min','variable_max')
      y_axis<-bquote("TCWC Residence Time (years)")
    }
    else if (v =='npp_wood') {
      plot_data <- gem_data[,c(1,2,7,13,19)]
      names(plot_data) <- c('model_plot_name','gem_region','variable','variable_min','variable_max')
      y_axis<-bquote("TCWC NPP ("~ Mg~C~ ha^-1~year^-1~")")
    }
    else {
      plot_data <- gem_data[,c(1,2,8,14,20)]
      names(plot_data) <- c('model_plot_name','gem_region','variable','variable_min','variable_max')
      y_axis<-bquote("TCWC Biomass ("~Mg~C~ ha^-1~")")
    }
    plot_data
    y_axis
  }
  for (i in gem_plots) {
    if (i =='allpahuayo') {
      p1 <- ggplot(data=plot_data[plot_data$gem_region==i,],aes(x = model_plot_name, y = variable))+
        ylim(c(min(plot_data$variable_min),max(plot_data$variable_max)))+
        geom_point(data=plot_data[plot_data$gem_region==i,],position=position_dodge(0.1),size=3) +
        geom_errorbar(data=plot_data[plot_data$gem_region==i,],aes(ymin=variable_min, ymax=variable_max), linewidth=1.5, width=.1, position=position_dodge(0.1)) +
        labs(title= "", x="Allpahuayo Plot/ CARDAMOM Model Variant",y=y_axis) + 
        theme_classic()+ 
        theme(text = element_text(size = 20))
    }
    else if (i =='caxiuana') {
      p2 <- ggplot(data=plot_data[plot_data$gem_region==i,],aes(x = model_plot_name, y = variable))+
        ylim(c(min(plot_data$variable_min),max(plot_data$variable_max)))+
        geom_point(data=plot_data[plot_data$gem_region==i,],position=position_dodge(0.1),size=3) +
        geom_errorbar(data=plot_data[plot_data$gem_region==i,],aes(ymin=variable_min, ymax=variable_max), linewidth=1.5, width=.1, position=position_dodge(0.1)) +
        labs(title= "", x="Caxiuana Plot/ CARDAMOM Model Variant",y="") + 
        theme_classic()+ 
        theme(text = element_text(size = 20))
    }
    else if (i =='kenia') {
      p3 <- ggplot(data=plot_data[plot_data$gem_region==i,],aes(x = model_plot_name, y = variable))+
        ylim(c(min(plot_data$variable_min),max(plot_data$variable_max)))+
        geom_point(data=plot_data[plot_data$gem_region==i,],position=position_dodge(0.1),size=3) +
        geom_errorbar(data=plot_data[plot_data$gem_region==i,],aes(ymin=variable_min, ymax=variable_max), linewidth=1.5, width=.1, position=position_dodge(0.1)) +
        labs(title= "", x="Kenia Plot/ CARDAMOM Model Variant",y=y_axis) + 
        theme_classic()+ 
        theme(text = element_text(size = 20))
    }
    else if (i =='tambopata') {
      p4 <- ggplot(data=plot_data[plot_data$gem_region==i,],aes(x = model_plot_name, y = variable))+
        ylim(c(min(plot_data$variable_min),max(plot_data$variable_max)))+
        geom_point(data=plot_data[plot_data$gem_region==i,],position=position_dodge(0.1),size=3) +
        geom_errorbar(data=plot_data[plot_data$gem_region==i,],aes(ymin=variable_min, ymax=variable_max), linewidth=1.5, width=.1, position=position_dodge(0.1)) +
        labs(title= "", x="Tambopata Plot/ CARDAMOM Model Variant",y="") + 
        theme_classic()+ 
        theme(text = element_text(size = 20))
    }
    # else {
    #   p5 <- ggplot(data=plot_data[plot_data$gem_region==i,],aes(x = model_plot_name, y = variable))+
    #     ylim(c(min(plot_data$variable_min),max(plot_data$variable_max)))+
    #     geom_point(data=plot_data[plot_data$gem_region==i,],position=position_dodge(0.1)) +
    #     geom_errorbar(data=plot_data[plot_data$gem_region==i,],aes(ymin=variable_min, ymax=variable_max), width=.1, position=position_dodge(0.1)) +
    #     labs(title= "", x="Tanguro Plot/ CARDAMOM Model Variant",y="") +
    #     theme_test()
    # }
    # print(ggarrange(p1, p2, p3, p4, p5, nrow = 2, ncol = 3))
    print(ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2))
  }
}
extract_plot_gem_var_df(gem_plots_model_data,gem_var[2],gem_region)


#   facet_wrap(.~gem_region) +


