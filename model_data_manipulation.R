#########Check data for editing##########
nc_sa_card_data <- nc_open('R:/ILAMB_beta_devel/RAINFOR_leeds_run/DATA/benchmark/CARDAMOM_monthly_1x1_SAmerica.nc') #gcel
nc_sa_card_data
nc_close(nc_sa_card_data)
nc_sa_inland_data <- nc_open('R:/ILAMB_beta_devel/RAINFOR_leeds_run/MODELS/INLAND/INLAND_monthly_1x1_SAmerica.nc') #inland
nc_sa_inland_data
nc_close(nc_sa_inland_data)
nc_sa_jules_data <- nc_open('R:/ILAMB_beta_devel/RAINFOR_leeds_run/MODELS/JULES/JULES_monthly_1x1_SAmerica.nc') #jules
nc_sa_jules_data
nc_close(nc_sa_jules_data)

nc_sa <- nc_open('G:/ILAMB_runs_output/CSSP_stippling/INLAND_monthly_SAmerica.nc') #gcel
nc_sa
nc_close(ncsa) #close nc file
