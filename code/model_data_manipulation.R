#########Check data for editing##########
nc_sa_card_data <- nc_open('R:/ILAMB_beta_devel/RAINFOR_leeds_run/DATA_BRAZIL/benchmark/CARDAMOM_Brazil_1x1_2001_2017_v1.0_wOW4.nc') #gcel
nc_sa_card_data
nc_close(nc_sa_card_data)
nc_sa_inland_data <- nc_open('R:/ILAMB_beta_devel/RAINFOR_leeds_run/3d_INLAND_monthly_1x1_Brazil.nc') #inland
nc_sa_inland_data
nc_close(nc_sa_inland_data)
nc_sa_jules_data <- nc_open('R:/ILAMB_beta_devel/RAINFOR_leeds_run/JULES_monthly_1x1_Brazil.nc') #jules
nc_sa_jules_data
nc_close(nc_sa_jules_data)
nc_sa_rainfor_data <- nc_open('R:/ILAMB_beta_devel/RAINFOR_leeds_run/correct_rainfor_outputwood_01_10_monthly.nc') #jules
nc_sa_rainfor_data
nc_close(nc_sa_rainfor_data)

nc_inland <- nc_open('R:/ILAMB_beta_devel/RAINFOR_leeds_run/MODELS_BRAZIL/JULES/JULES_monthly_1x1_Brazil.nc', write=TRUE )
#ncatt_put( nc_inland, "alwood", "units", "g.m-2.d-1")

old_varname <-'litC'
new_varname <- 'NPP_wood'

nc_inland <- ncvar_rename( nc_inland, old_varname, new_varname )
nc_close(nc_inland) #close nc file


nc_sa_card_data <- nc_open('R:/ILAMB_beta_devel/RAINFOR_leeds_run/DATA_BRAZIL/benchmark/CARDAMOM_Brazil_1x1_2001_2017_v1.0_wOW4.nc', write=TRUE) #gcel
nc_sa_card_data
ncatt_put(nc_sa_card_data, "outputwood", "units", "g.m-2.d-1")
ncatt_put(nc_sa_card_data, "outputwood_25pc", "units", "g.m-2.d-1")
ncatt_put(nc_sa_card_data, "outputwood_75pc", "units", "g.m-2.d-1")
nc_close(nc_sa_card_data)

nc_sa_inland_data <- nc_open('R:/ILAMB_beta_devel/RAINFOR_leeds_run/3d_INLAND_monthly_1x1_Brazil.nc', write=TRUE) #inland
ncatt_put(nc_sa_inland_data, "fallw", "units", "g.m-2.d-1")
old_varname <-'fallw'
new_varname <- 'outputwood'
nc_sa_inland_data <- ncvar_rename(nc_sa_inland_data, old_varname, new_varname )
nc_sa_inland_data
nc_close(nc_sa_inland_data)

nc_sa <- nc_open('G:/ILAMB_runs_output/CSSP_stippling/MODELS/JULES/JULES_monthly_1x1_Brazil.nc') #gcel
nc_sa
nc_close(nc_sa) #close nc file
sa_data <- stack('G:/ILAMB_runs_output/CSSP_stippling/MODELS/JULES/JULES_monthly_1x1_Brazil.nc',varname="npptot")
sa_data
plot(sa_data$X2001.01.15)
summary(sa_data$X2001.01.01)
#litC,alwood,NPP_wood,nppwood

#check masks
nc_trial <- nc_open('R:/ILAMB_beta_devel/RAINFOR_leeds_run/masks/ilamb-mask-FGUYAMAZON.nc') #gcel
nc_trial
nc_close(nc_trial)

bol_mask_data <- raster('R:/ILAMB_beta_devel/RAINFOR_leeds_run/masks/ilamb-mask-BRAZILAMAZON.nc')
bol_mask_data
plot(bol_mask_data)
