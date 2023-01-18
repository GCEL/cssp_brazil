#All codes being Used
#install packages----
# library(ncdf4); library(raster); library(dplyr); library(ggplot2);library(ggpubr);library(quantreg);library(ggpp);library(rgeos);library(ggpmisc);library(rgdal);library(Metrics)

#done
#data prep pre-cursors----
# amazonia_subset <- shapefile("R:/cssp_brazil/cssp_brazil_R/data/amazonia_subset.shp")
# nbe_var_name <- c('NBE')
# mod_var_nbe <- c('esa_cci_biomass_only','biomass_annual_productivity','esa_cci_agb_nbe','biomass_annual_productivity_nbe')
# nbe_time_period <- c('2015_2019')
# prefix_nbe <- 'R://ILAMB_beta_devel/RAINFOR_leeds_run/cssp_brazil_amazon_run/DATA/benchmark/Amazon_subset_'
# suffix_nbe <- '_2000_updated_2019.nc'

# nbe_benchmark_data<- nbe_stack_15_19_monthly
# save(nbe_benchmark_data,file = "R:/cssp_brazil/cssp_brazil_R/data/nbe_benchmark_data.RData")
# plot(nbe_benchmark_data$X2015.01.01)

#done
########################################################################
####new code to include variable, model and region######################
########################################################################
#extract subset of wood, npp and mortality from rainfor benchmark data
extract_subset <- function (region,reference){
  target = raster(crs = ("+init=epsg:4326"), ext = extent(c(-79.25504, -43.39195, -20.42212, 10.45082)), resolution = c(0.9692725, 0.9647794))
  reference<-resample(reference, target, method = "bilinear", na.rm=TRUE)
  masked1 <- mask(reference, region)
  masked2 <- masked1 > -Inf
  poly <- rasterToPolygons(masked2, dissolve=TRUE)
  data_region <- mask(reference,poly)
  data_region_mean <- stackApply(data_region, indices =  rep(1,nlayers(data_region)), fun = "mean")
  return(data_region_mean)
}
# check_extract_subset<-extract_subset(amazonia_subset,nbe_benchmark_data)
# plot(check_extract_subset)

#done
###################################################
################NBE data######################
###################################################
#extract subset of wood, npp and mortality from models
extract_nbe_subset_all_models <- function (region,cardamom_var,model_variant,reference,tp){
  for (i in model_variant) {
    data <- brick(paste(prefix_nbe,model_variant[1],suffix_nbe,sep=""),varname=cardamom_var)
    data_1519<-data[[which(getZ(data) >= as.Date("2015-01-01") & getZ(data) <= as.Date("2019-12-01"))]]
    data_1519_mean <- stackApply(data_1519, indices =  rep(1,nlayers(data_1519)), fun = "mean")
  }
  return(data_1519_mean)
}

collate_all_models_nbe <- function (region,cardamom_var,model_variant,reference,tp){
  for (i in model_variant) {
    nbe_raster <- extract_nbe_subset_all_models(region,cardamom_var,i,reference,tp)
  }
  return(nbe_raster)
}
# , xlim=c(0,25000), ylim=c(0,25000)
res_df_merge_plot_nbe <- function (region,cardamom_var,model_variant,reference,tp){
  for (i in model_variant) {
    benchmark <- extract_subset(region,reference)
    df_benchmark <- as.data.frame(benchmark, xy=TRUE, na.rm=TRUE)
    if (i=='esa_cci_biomass_only') {
      model_res <- extract_nbe_subset_all_models(region,cardamom_var,i,reference,tp)
      # model_res <- resample(model,benchmark)
      model_res_df <- as.data.frame(model_res, xy=TRUE, na.rm=TRUE)
      model_benchmark_df <- merge(df_benchmark, model_res_df, by = 'row.names')
      model_benchmark_df <- model_benchmark_df[,c(5,6,4,7)]
      names(model_benchmark_df)<-c('x','y','benchmark',i)
      model_benchmark_df<-na.omit(model_benchmark_df)
      model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]),digits=2)
      assign(paste("plot", i, sep = "_"),
             plot(model_benchmark_df[,4], model_benchmark_df[,3], main="",xlab="Default CARDAMOM", ylab="GEOSCHEM benchmark", pch=19,
                  cex.lab=1.5, cex.axis=1.5, xlim=c(-0.7,0.7), ylim=c(-0.4,0.4))
             +abline(coef = c(0,1),col='red', lwd=3)
             +text(0.4, 0, paste('RMSE = ',model_rmse, sep=''),cex=1.5)
             )
      }
    else if (i=='biomass_annual_productivity') {
      model_res <- extract_nbe_subset_all_models(region,cardamom_var,i,reference,tp)
      # model_res <- resample(model,benchmark)
      model_res_df <- as.data.frame(model_res, xy=TRUE, na.rm=TRUE)
      model_benchmark_df <- merge(df_benchmark, model_res_df, by = 'row.names')
      model_benchmark_df <- model_benchmark_df[,c(5,6,4,7)]
      names(model_benchmark_df)<-c('x','y','benchmark',i)
      model_benchmark_df<-na.omit(model_benchmark_df)
      model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]),digits=2)
      assign(paste("plot", i, sep = "_"),
             plot(model_benchmark_df[,4], model_benchmark_df[,3], main="",xlab="New CARDAMOM Benchmark", ylab="GEOSCHEM benchmark", pch=19, 
                  cex.lab=1.5, cex.axis=1.5, xlim=c(-0.7,0.7), ylim=c(-0.4,0.4))
             +abline(coef = c(0,1),col='red', lwd=3)
             +text(0.4, 0, paste('RMSE = ',model_rmse, sep=''),cex=1.5)
             )
      }
    else if (i=='esa_cci_agb_nbe') {
      model_res <- extract_nbe_subset_all_models(region,cardamom_var,i,reference,tp)
      # model_res <- resample(model,benchmark)
      model_res_df <- as.data.frame(model_res, xy=TRUE, na.rm=TRUE)
      model_benchmark_df <- merge(df_benchmark, model_res_df, by = 'row.names')
      model_benchmark_df <- model_benchmark_df[,c(5,6,4,7)]
      names(model_benchmark_df)<-c('x','y','benchmark',i)
      model_benchmark_df<-na.omit(model_benchmark_df)
      model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]),digits=2)
      assign(paste("plot", i, sep = "_"),
             plot(model_benchmark_df[,4], model_benchmark_df[,3], main="",xlab="Default CARDAMOM with GEOSCHEM NBE", ylab="GEOSCHEM benchmark", pch=19, 
                  cex.lab=1.5, cex.axis=1.5, xlim=c(-0.7,0.7), ylim=c(-0.4,0.4))
             +abline(coef = c(0,1),col='red', lwd=3)
             +text(0.4, 0, paste('RMSE = ',model_rmse, sep=''),cex=1.5)
             )
      }
    else if (i=='biomass_annual_productivity_nbe') {
      model_res <- extract_nbe_subset_all_models(region,cardamom_var,i,reference,tp)
      # model_res <- resample(model,benchmark)
      model_res_df <- as.data.frame(model_res, xy=TRUE, na.rm=TRUE)
      model_benchmark_df <- merge(df_benchmark, model_res_df, by = 'row.names')
      model_benchmark_df <- model_benchmark_df[,c(5,6,4,7)]
      names(model_benchmark_df)<-c('x','y','benchmark',i)
      model_benchmark_df<-na.omit(model_benchmark_df)
      model_rmse<- round(rmse(model_benchmark_df[,3],model_benchmark_df[,4]),digits=2)
      assign(paste("plot", i, sep = "_"),
             plot(model_benchmark_df[,4], model_benchmark_df[,3], main="",xlab="New CARDAMOM Benchmark with GEOSCHEM NBE", ylab="GEOSCHEM benchmark", pch=19, 
                  cex.lab=1.5, cex.axis=1.5, xlim=c(-0.7,0.7), ylim=c(-0.4,0.4))
             +abline(coef = c(0,1),col='red', lwd=3)
             +text(0.4, 0, paste('RMSE = ',model_rmse, sep=''),cex=1.5)
             )
      }
    }
  par(mfrow = c(2,2),oma = c(0, 0, 2, 0))
  plot_esa_cci_biomass_only
  plot_biomass_annual_productivity
  plot_esa_cci_agb_nbe
  plot_biomass_annual_productivity_nbe
  mtext(expression("NBE Comparison 2015-2019 gC.m-2.d-1"), line=-2, side=3, outer=TRUE, cex=1.5)
}

# trial_s1_default<-extract_nbe_subset_all_models(amazonia_subset,nbe_var_name,mod_var_nbe[1],nbe_benchmark_data)
# collate_all_models_nbe(amazonia_subset,nbe_var_name,mod_var_nbe[4],nbe_benchmark_data)
# par(mfrow = c(2,4));res_df_merge_plot_cwood(amazonia_subset,compare_var_names[1],mod_var,biomass_amazon_gCm2)
# res_df_merge_plot_nbe(amazonia_subset,nbe_var_name,mod_var_nbe,nbe_benchmark_data)
#done
#######################################################
###############extract model parameters################
#######################################################
# prefix_new <- 'M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/'
# midfix <- '_subset/RESULTS_PROCESSED/'
# suffix_new <- '_subset_stock_flux.RData'
# 
# mod_var_nbe <- c('esa_cci_agb_only','Rainfor_woody_biomass_annual_productivity',
#                  'esa_cci_agb_nbe','Rainfor_woody_biomass_annual_productivity_nbe'
# )
# mod_var_abr_nbe <- c('default','benchmark','default_nbe','benchmark_nbe')

model_pars_list <- list()
# model_pars_list_sd <- list()

extract_parameters_all_models <- function (model_variant,model_variant_abr){
  for (i in seq_along(model_variant)) {
    load(paste(prefix_new,model_variant[i],midfix,model_variant[i],suffix_new,sep=""))
    data <- grid_output$parameters[,,,4]
    data_par_mean <- colMeans(data, dims = 2,na.rm=T)
    data_par_mean[34] <- model_variant_abr[i]
    model_pars_list[[i]]<-data_par_mean
  }
  return(model_pars_list)
}
plot_all_parameters <- function(data,parameter_names) {
  for (i in 2:ncol(data)){
    print(ggplot(data, aes(x = model_var,y=data[,i])) +
            geom_bar(stat = "identity", width=0.5)+
            labs(title=parameter_names[i], y = parameter_names[i])+
            scale_fill_grey(start = 0, end = 1) +
            theme_bw()+
            theme(text = element_text(size = 16),panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")))
  }
}
round_df <- function(x, digits) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}


# extract_parameters_all_models_sd <- function (model_variant,model_variant_abr){
#   for (i in seq_along(model_variant)) {
#     load(paste(prefix_new,model_variant[i],midfix,model_variant[i],suffix_new,sep=""))
#     data <- grid_output$parameters[,,,4]
#     # for (j in seq_along(data[37,32,])){
#     #   data_par_sd <- sd(data[37,32,j],na.rm=T)
#     # }
#     # data_par_sd[34] <- model_variant_abr[i]
#     # model_pars_list_sd[[i]]<-data_par_sd
#     data_par_sd <- sd(data, dims = 2,na.rm=T)
#     data_par_sd[34] <- model_variant_abr[i]
#     model_pars_list[[i]]<-data_par_mean
#   }
#   return(model_pars_list_sd)
# }

# model_pars_df_sd<-as.data.frame(t(as.data.frame(extract_parameters_all_models_sd(mod_var_new,mod_var_abr))), stringsAsFactors = FALSE)
# model_pars_df_sd[] <- lapply(model_pars_df_sd, type.convert, as.is = TRUE)

# par_names <- c('Lit2SOM (day_1)','GPP%Ra','NPP_fol_frac','NPP_root_frac','Leaf lifespan','TO Wood','TO Roots','Mineralise Lit','SOM2Rh',
#                'Temp fac','Canopy Eff (gC/m2leaf/day)','Max bud burst day','NPP_lab_frac','Labile_release_period','max leaf fall',
#                'leaf fall period','LMA (gC/m2)','Labile_initial','Foliar_initial','Roots_initial','Wood_initial','Lit_initial',
#                'SOM_initial','Soil Water_initial','Cwood_coarseR_frac','CoarseR Biomass 50% Max Depth','Max Root Depth',
#                'Res factor Burned Cstocks','CCF Fol','CCF Wood and FineR','CCF Soil','CCF Fol and FineR Lit','unknown','model_var')
# # summary(model_pars_df)
#
# param_names <- c('model_var','Lit2SOM (day_1)','GPP%Ra','NPP_fol_frac','NPP_root_frac','Leaf lifespan','TO Wood','TO Roots','Mineralise Lit','SOM2Rh',
#                'Temp fac','Canopy Eff (gC/m2leaf/day)','Max bud burst day','NPP_lab_frac','Labile_release_period','max leaf fall',
#                'leaf fall period','LMA (gC/m2)','Labile_initial','Foliar_initial','Roots_initial','Wood_initial','Lit_initial',
#                'SOM_initial','Soil Water_initial','Cwood_coarseR_frac','CoarseR Biomass 50% Max Depth','Max Root Depth',
#                'Res factor Burned Cstocks','CCF Fol','CCF Wood and FineR','CCF Soil','CCF Fol and FineR Lit')

# model_pars_df<-as.data.frame(t(data.frame(extract_parameters_all_models(mod_var_new,mod_var_abr))), stringsAsFactors = FALSE)
# model_pars_nbe_df<-as.data.frame(t(data.frame(extract_parameters_all_models(mod_var_nbe,mod_var_abr_nbe))), stringsAsFactors = FALSE)
# model_pars_nbe_df[] <- lapply(model_pars_nbe_df, type.convert, as.is = TRUE)
# rownames(model_pars_nbe_df)<-  seq_along(model_pars_nbe_df[,1])
# colnames(model_pars_nbe_df)<-  par_names
# model_pars_nbe_df$model_var<-as.factor(model_pars_nbe_df$model_var)
# model_pars_nbe_df<-model_pars_nbe_df[,c(34,1:32)]
# new_model_pars_df<-round_df(model_pars_df[,2:33], 3)
# rownames(new_model_pars_df)<-mod_var_new
# print(new_model_pars_df)
# plot_all_parameters(model_pars_df,param_names)
