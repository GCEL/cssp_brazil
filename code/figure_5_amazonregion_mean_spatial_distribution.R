#load R libraries
library(ncdf4); library(raster); library(dplyr); library(ggplot2);library(ggpubr);library(quantreg);library(ggpp);library(rgeos);library(ggpmisc);library(rgdal);library(Metrics)

cardamom_model_exp <- c('C','R0','R1','R2') #parameters for comparisons
ecoregion_names <- c('northwest','southwest','brazilshield','eastcentral','guyanashield') #ecoregion names
ci_names <- c('2pt5','mean','97pt5') #Confidence interval names
var_spadis_names <- c('gpp','cue','npp','nppwood','outputwood','wood','mrtwood') #CARDAMOM variables  needed
new_mod_var <- c("esa_cci_agb","rainfor_biomass_annual","rainfor_biomass_productivity_2005","rainfor_biomass_annual_productivity") ##CARDAMOM run variants

for (e in ecoregion_names) {###loop run for each ecoregion
  for (i in cardamom_model_exp) {#loop over all 4 CARDAMOM run variants found in amazonia_ifl_cardamom_runs zip folder
    if (i =='C') {
      load("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/amazonia_ifl_esa_cci_agb_nomngt//infofile.RData")
      load(paste("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/amazonia_ifl_esa_cci_agb_nomngt/RESULTS_PROCESSED/","amazonia_ifl_esa_cci_agb_nomngt","_stock_flux.RData",sep=""))
      print('CARDAMOM C')
    }
    else if (i =='R0') {
      load("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/amazonia_ifl_rainfor_biomass_annual_nomngt//infofile.RData")
      load(paste("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/amazonia_ifl_rainfor_biomass_annual_nomngt/RESULTS_PROCESSED/","amazonia_ifl_rainfor_biomass_annual_nomngt","_stock_flux.RData",sep=""))
      print('CARDAMOM R0')
    }
    else if (i =='R1') {
      load("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/amazonia_ifl_rainfor_biomass_productivity_2005_nomngt//infofile.RData")
      load(paste("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/amazonia_ifl_rainfor_biomass_productivity_2005_nomngt/RESULTS_PROCESSED/","amazonia_ifl_rainfor_biomass_productivity_2005_nomngt","_stock_flux.RData",sep=""))
      print('CARDAMOM R1')
    }
    else {
      load("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/amazonia_ifl_rainfor_biomass_annual_productivity_nomngt//infofile.RData")
      load(paste("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/amazonia_ifl_rainfor_biomass_annual_productivity_nomngt/RESULTS_PROCESSED/","amazonia_ifl_rainfor_biomass_annual_productivity_nomngt","_stock_flux.RData",sep=""))
      print('CARDAMOM R2')
    }
    
    model_reg_par_mean_spadist_df <- data.frame(matrix(nrow=35,ncol=4)) #create df to accomodate data
    colnames(model_reg_par_mean_spadist_df) <- cardamom_model_exp #column names= CARDAMOM variants
    model_reg_par_mean_spadist_df$vars <- rep(var_spadis_names,length(ecoregion_names)) #create column for variable name for each eco region
    model_reg_par_mean_spadist_df$region <- c(rep(ecoregion_names[1],length(var_spadis_names)), #create column for each ecoregion
                                              rep(ecoregion_names[2],length(var_spadis_names)),
                                              rep(ecoregion_names[3],length(var_spadis_names)),
                                              rep(ecoregion_names[4],length(var_spadis_names)),
                                              rep(ecoregion_names[5],length(var_spadis_names)))
  if (e =='northwest') {
    region_mask <- amazonia_ecoregions_matrix==1
    print(e)
  }
  else if (e =='southwest') {
    region_mask <- amazonia_ecoregions_matrix==2
    print(e)
  }
  else if (e =='brazilshield') {
    region_mask <- amazonia_ecoregions_matrix==3
    print(e)
  }
  else if (e =='eastcentral') {
    region_mask <- amazonia_ecoregions_matrix==4
    print(e)
  }
  else {
    region_mask <- amazonia_ecoregions_matrix==5
    print(e)
  }
#Masking phase
# NATURAL FLUXES
grid_output$mean_gpp_gCm2day[,,4][!region_mask]<-NA
grid_output$mean_cue[,,4][!region_mask]<-NA
grid_output$mean_npp_gCm2day[,,4][!region_mask]<-NA
grid_output$mean_alloc_wood_gCm2day[,,4][!region_mask]<-NA
grid_output$mean_wood_to_litter_gCm2day[,,4][!region_mask]<-NA
# STOCKS
wood <- grid_output$mean_wood_gCm2[,,4][!region_mask]<-NA
# MEAN TRANSIENT TIMES
mrtwood <- grid_output$MTT_wood_years[,,4][!region_mask]<-NA

#create empty df to populate for each eco-region
reg_par_mean_spadist_df <- data.frame(matrix(nrow=3,ncol=8))
colnames(reg_par_mean_spadist_df) <- c(var_spadis_names,'region')
rownames(reg_par_mean_spadist_df) <- ci_names
reg_par_mean_spadist_df$region <- rep(e,nrow(reg_par_mean_spadist_df))

#create spatial distribution (mean and confidence interval)
#gpp (t/ha/year)
gpp_list<-unlist(as.list(grid_output$mean_gpp_gCm2day[,,4]))* (365.25/100)
gpp_list<-gpp_list[!is.na(gpp_list)]
gpp_list_model <- lm (gpp_list~1)
gpp_list_ci<-confint(gpp_list_model,level=0.95)
gpp_list_mean<-mean(gpp_list)
gpp_spadis<- round(c(gpp_list_ci[,1],gpp_list_mean[1],gpp_list_ci[,2]),digits=2) #[2.5, mean. 97.5]
gpp_range <- paste(gpp_spadis[2]," (",gpp_spadis[1],"-",gpp_spadis[3],")", sep="")
reg_par_mean_spadist_df[,1] <- gpp_spadis


#cue
cue_list<-unlist(as.list(grid_output$mean_cue[,,4]))
cue_list<-cue_list[!is.na(cue_list)]
cue_list_model <- lm (cue_list~1)
cue_list_ci<-confint(cue_list_model,level=0.95)
cue_list_mean<-mean(cue_list)
cue_spadis<- round(c(cue_list_ci[,1],cue_list_mean[1],cue_list_ci[,2]),digits=3)#[2.5, mean. 97.5]
cue_range <- paste(cue_spadis[2]," (",cue_spadis[1],"-",cue_spadis[3],")", sep="")
reg_par_mean_spadist_df[,2] <- cue_spadis

#npp (t/ha/year)
npp_list<-unlist(as.list(grid_output$mean_npp_gCm2day[,,4]))* (365.25/100)
npp_list<-npp_list[!is.na(npp_list)]
npp_list_model <- lm (npp_list~1)
npp_list_ci<-confint(npp_list_model,level=0.95)
npp_list_mean<-mean(npp_list)
npp_spadis<- round(c(npp_list_ci[,1],npp_list_mean[1],npp_list_ci[,2]),digits=2)#[2.5, mean. 97.5]
npp_range <- paste(npp_spadis[2]," (",npp_spadis[1],"-",npp_spadis[3],")", sep="")
reg_par_mean_spadist_df[,3] <- npp_spadis

#nppwood (t/ha/year)
nppwood_list<-unlist(as.list(grid_output$mean_alloc_wood_gCm2day[,,4]))* (365.25/100)
nppwood_list<-nppwood_list[!is.na(nppwood_list)]
nppwood_list_model <- lm (nppwood_list~1)
nppwood_list_ci<-confint(nppwood_list_model,level=0.95)
nppwood_list_mean<-mean(nppwood_list)
nppwood_spadis<- round(c(nppwood_list_ci[,1],nppwood_list_mean[1],nppwood_list_ci[,2]),digits=2)#[2.5, mean. 97.5]
nppwood_range <- paste(nppwood_spadis[2]," (",nppwood_spadis[1],"-",nppwood_spadis[3],")", sep="")
reg_par_mean_spadist_df[,4] <- nppwood_spadis

#outputwood (t/ha/year)
outputwood_list<-unlist(as.list(grid_output$mean_wood_to_litter_gCm2day[,,4]))* (365.25/100)
outputwood_list<-outputwood_list[!is.na(outputwood_list)]
outputwood_list_model <- lm (outputwood_list~1)
outputwood_list_ci<-confint(outputwood_list_model,level=0.95)
outputwood_list_mean<-mean(outputwood_list)
outputwood_spadis<- round(c(outputwood_list_ci[,1],outputwood_list_mean[1],outputwood_list_ci[,2]),digits=2)#[2.5, mean. 97.5]
outputwood_range <- paste(outputwood_spadis[2]," (",outputwood_spadis[1],"-",outputwood_spadis[3],")", sep="")
reg_par_mean_spadist_df[,5] <- outputwood_spadis

#wood (t/ha)
wood_list<-unlist(as.list(grid_output$mean_wood_gCm2[,,4])) /100
wood_list<-wood_list[!is.na(wood_list)]
wood_list_model <- lm (wood_list~1)
wood_list_ci<-confint(wood_list_model,level=0.95)
wood_list_mean<-mean(wood_list)
wood_spadis<- round(c(wood_list_ci[,1],wood_list_mean[1],wood_list_ci[,2]),digits=2)#[2.5, mean. 97.5]
wood_range <- paste(wood_spadis[2]," (",wood_spadis[1],"-",wood_spadis[3],")", sep="")
reg_par_mean_spadist_df[,6] <- wood_spadis

#mrtwood (year)
mrtwood_list<-unlist(as.list(grid_output$MTT_wood_years[,,4]))
mrtwood_list<-mrtwood_list[!is.na(mrtwood_list)]
mrtwood_list_model <- lm (mrtwood_list~1)
mrtwood_list_ci<-confint(mrtwood_list_model,level=0.95)
mrtwood_list_mean<-mean(mrtwood_list)
mrtwood_spadis<- round(c(mrtwood_list_ci[,1],mrtwood_list_mean[1],mrtwood_list_ci[,2]),digits=2)#[2.5, mean. 97.5]
mrtwood_range <- paste(mrtwood_spadis[2]," (",mrtwood_spadis[1],"-",mrtwood_spadis[3],")", sep="")
reg_par_mean_spadist_df[,7] <- mrtwood_spadis

# print(reg_par_mean_spadist_df)
assign(paste(i,e,"_par_mean_spadist_df",sep="_"),reg_par_mean_spadist_df)
regional_var <- c(gpp_range,cue_range,npp_range,nppwood_range,outputwood_range,wood_range,mrtwood_range)
assign(paste(i,e,"par_mean_spadist_list",sep="_"),regional_var)
print(paste('CARDAMOM ',i, " ",e, " ",'done',sep=""))
}
print(paste('CARDAMOM ',i, " ",'done',sep=""))
}

model_reg_par_mean_spadist_df$C <- c(C_northwest_par_mean_spadist_list,
                                     C_southwest_par_mean_spadist_list,
                                     C_brazilshield_par_mean_spadist_list,
                                     C_eastcentral_par_mean_spadist_list,
                                     C_guyanashield_par_mean_spadist_list)
model_reg_par_mean_spadist_df$R0 <- c(R0_northwest_par_mean_spadist_list,
                                      R0_southwest_par_mean_spadist_list,
                                      R0_brazilshield_par_mean_spadist_list,
                                      R0_eastcentral_par_mean_spadist_list,
                                      R0_guyanashield_par_mean_spadist_list)
model_reg_par_mean_spadist_df$R1 <- c(R1_northwest_par_mean_spadist_list,
                                      R1_southwest_par_mean_spadist_list,
                                      R1_brazilshield_par_mean_spadist_list,
                                      R1_eastcentral_par_mean_spadist_list,
                                      R1_guyanashield_par_mean_spadist_list)
model_reg_par_mean_spadist_df$R2 <- c(R2_northwest_par_mean_spadist_list,
                                      R2_southwest_par_mean_spadist_list,
                                      R2_brazilshield_par_mean_spadist_list,
                                      R2_eastcentral_par_mean_spadist_list,
                                      R2_guyanashield_par_mean_spadist_list)

write.csv(model_reg_par_mean_spadist_df,"all_model_reg_par_mean_spadist_df.csv")




