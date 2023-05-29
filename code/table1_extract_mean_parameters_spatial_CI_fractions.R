
#install packages
library(ncdf4); library(raster); library(dplyr); library(ggplot2);library(ggpubr);library(quantreg);library(ggpp);library(rgeos);library(ggpmisc);library(rgdal);library(Metrics);library(gridExtra);library(grid)

#define terms
cardamom_model_exp <- c('C','R0','R1','R2') #representation of each CARDAMOM run variants
par_meanvi_names<-c('cue','nppfol_frac','npproots_frac','nppwood_frac','mrtwood','mrtroots') #parameters to extract
new_mod_var <- c("esa_cci_agb","rainfor_biomass_annual","rainfor_biomass_productivity_2005","rainfor_biomass_annual_productivity") #index of the CARDAMOM run variations

par_med_95ci_df <- data.frame(matrix(nrow=6,ncol=4)) #rows=parameters; columns=model variant names
colnames(par_med_95ci_df) <- cardamom_model_exp
rownames(par_med_95ci_df) <- par_meanvi_names

for (i in cardamom_model_exp) {#loop over all 4 CARDAMOM run variants found in amazonia_ifl_cardamom_runs zip folder
    if (i =='C') {
      load("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/amazonia_ifl_esa_cci_agb_nomngt//infofile.RData")
      load(paste("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/amazonia_ifl_esa_cci_agb_nomngt/RESULTS_PROCESSED/","amazonia_ifl_esa_cci_agb_nomngt","_stock_flux.RData",sep=""))
      print(paste('CARDAMOM',i,sep=" "))
    }
    else if (i =='R0') {
      load("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/amazonia_ifl_rainfor_biomass_annual_nomngt//infofile.RData")
      load(paste("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/amazonia_ifl_rainfor_biomass_annual_nomngt/RESULTS_PROCESSED/","amazonia_ifl_rainfor_biomass_annual_nomngt","_stock_flux.RData",sep=""))
      print(paste('CARDAMOM',i,sep=" "))
    }
    else if (i =='R1') {
      load("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/amazonia_ifl_rainfor_biomass_productivity_2005_nomngt//infofile.RData")
      load(paste("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/amazonia_ifl_rainfor_biomass_productivity_2005_nomngt/RESULTS_PROCESSED/","amazonia_ifl_rainfor_biomass_productivity_2005_nomngt","_stock_flux.RData",sep=""))
      print(paste('CARDAMOM',i,sep=" "))
    }
    else {
      load("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/amazonia_ifl_rainfor_biomass_annual_productivity_nomngt//infofile.RData")
      load(paste("M://CARDAMOM/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/amazonia_ifl_rainfor_biomass_annual_productivity_nomngt/RESULTS_PROCESSED/","amazonia_ifl_rainfor_biomass_annual_productivity_nomngt","_stock_flux.RData",sep=""))
      print(paste('CARDAMOM',i,sep=" "))
    }

  quantiles_wanted<- c(1,4,7) #quantiles of 2.5, 50 and 97.5%
  dp <-2 #decimal places
# Extract parameters from grid_output
cue <- format(round(apply(grid_output$mean_cue[,,quantiles_wanted],3,mean, na.rm=TRUE), digits = 2), nsmall = dp) 
nppfol <- format(round(apply(grid_output$NPP_foliage_fraction[,,quantiles_wanted],3,mean, na.rm=TRUE), digits = 2), nsmall = dp) 
npproots <- format(round(apply(grid_output$NPP_roots_fraction[,,quantiles_wanted],3,mean, na.rm=TRUE), digits = 2), nsmall = dp) 
nppwood <- format(round(apply(grid_output$NPP_wood_fraction[,,quantiles_wanted],3,mean, na.rm=TRUE), digits = 2), nsmall = dp)
mrtwood <- format(round(apply(grid_output$MTT_wood_years[,,quantiles_wanted],3,mean, na.rm=TRUE), digits = 2), nsmall = dp)
mrtroots <- format(round(apply(grid_output$MTT_roots_years[,,quantiles_wanted],3,mean, na.rm=TRUE), digits = 2), nsmall = dp)

#create spa dist
#cue
cue_range <- paste(cue[2]," (",cue[1],"-",cue[3],")", sep="")
par_med_95ci_df[1,i] <- cue_range

#nppfol
nppfol_range <- paste(nppfol[2]," (",nppfol[1],"-",nppfol[3],")", sep="")
par_med_95ci_df[2,i] <- nppfol_range

#npproots
npproots_range <- paste(npproots[2]," (",npproots[1],"-",npproots[3],")", sep="")
par_med_95ci_df[3,i] <- npproots_range

#nppwood
nppwood_range <- paste(nppwood[2]," (",nppwood[1],"-",nppwood[3],")", sep="")
par_med_95ci_df[4,i] <- nppwood_range

#mrtroots
mrtroots_range <- paste(mrtroots[2]," (",mrtroots[1],"-",mrtroots[3],")", sep="")
par_med_95ci_df[5,i] <- mrtroots_range

#mrtwood
mrtwood_range <- paste(mrtwood[2]," (",mrtwood[1],"-",mrtwood[3],")", sep="")
par_med_95ci_df[6,i] <- mrtwood_range

}

par_med_95ci_df
# write.csv(par_med_95ci_df,"par_med_95ci_df.csv")




