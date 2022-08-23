#########################################################################
#######################Packages and Data#################################
#########################################################################
library(hexbin);library(RColorBrewer);library(gridExtra);library(grid);library(ggplot2);library(lattice)
str(all_pixels_compare_data_zero); summary(all_pixels_compare_data)

y <- all_pixels_compare_data[all_pixels_compare_data$region_name=='amazonia',]$gpp
x <- all_pixels_compare_data[all_pixels_compare_data$region_name=='amazonia',]$lai
z <- all_pixels_compare_data[all_pixels_compare_data$region_name=='amazonia',]$cica

amazonia_bin_gpplai<-hexbin(x,y, xbins = 40)
amazonia_bin_gppcica<-hexbin(z,y, xbins = 40)
amazonia_colors<-colorRampPalette(rev(brewer.pal(11,'Spectral')))
plot(amazonia_bin_gpplai, main="Amazonia GPP~LAI Bin",xlab='LAI',ylab='GPP', colramp=amazonia_colors) 
plot(amazonia_bin_gppcica, main="Amazonia GPP~CICA Bin",xlab='CiCa',ylab='GPP', colramp=amazonia_colors) 


hexbin_plot_fun <- function(comp_data) {
  regions <- unique(comp_data$region_name)
  colors <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
  par(mfrow=c(2,3))
  for (i in regions) {
    if(i=='amazon_nw'){
      y <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$gpp
      x <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$lai
      z <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$cica
      bin_gpplai<-hexbin(x,y, xbins = 40)
      bin_gppcica<-hexbin(z,y, xbins = 40)
      gpplai_nw <- plot(bin_gpplai, main="North West Amazonia GPP~LAI Bin",xlab='LAI',ylab='GPP', colramp=colors) 
      gppcica_nw <- plot(bin_gppcica, main="North West Amazonia GPP~CICA Bin",xlab='CiCa',ylab='GPP', colramp=colors) 
      }
    else if (i=='amazon_sw') {
      y <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$gpp
      x <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$lai
      z <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$cica
      bin_gpplai<-hexbin(x,y, xbins = 40)
      bin_gppcica<-hexbin(z,y, xbins = 40)
      gpplai_sw <- plot(bin_gpplai, main="South West Amazonia GPP~LAI Bin",xlab='LAI',ylab='GPP', colramp=colors) 
      gppcica_sw <- plot(bin_gppcica, main="South West Amazonia GPP~CICA Bin",xlab='CiCa',ylab='GPP', colramp=colors) 
    }
    else if (i=='amazon_ec') {
      y <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$gpp
      x <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$lai
      z <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$cica
      bin_gpplai<-hexbin(x,y, xbins = 40)
      bin_gppcica<-hexbin(z,y, xbins = 40)
      gpplai_ec <- plot(bin_gpplai, main="East central Amazonia GPP~LAI Bin",xlab='LAI',ylab='GPP', colramp=colors) 
      gppcica_ec <- plot(bin_gppcica, main="East central Amazonia GPP~CICA Bin",xlab='CiCa',ylab='GPP', colramp=colors) 
    }
    else if (i=='amazon_bs') {
      y <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$gpp
      x <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$lai
      z <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$cica
      bin_gpplai<-hexbin(x,y, xbins = 40)
      bin_gppcica<-hexbin(z,y, xbins = 40)
      gpplai_bs <- plot(bin_gpplai, main="Brazil Shield GPP~LAI Bin",xlab='LAI',ylab='GPP', colramp=colors) 
      gppcica_bs <- plot(bin_gppcica, main="Brazil Shield GPP~CICA Bin",xlab='CiCa',ylab='GPP', colramp=colors) 
    }
    else if (i=='amazon_gs') {
      y <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$gpp
      x <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$lai
      z <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$cica
      bin_gpplai<-hexbin(x,y, xbins = 40)
      bin_gppcica<-hexbin(z,y, xbins = 40)
      gpplai_gs <- plot(bin_gpplai, main="Guyana Shield GPP~LAI Bin",xlab='LAI',ylab='GPP', colramp=colors) 
      gppcica_gs <- plot(bin_gppcica, main="Guyana Shield GPP~CICA Bin",xlab='CiCa',ylab='GPP', colramp=colors) 
    }
    else{
      y <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$gpp
      x <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$lai
      z <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$cica
      bin_gpplai<-hexbin(x,y, xbins = 40)
      bin_gppcica<-hexbin(z,y, xbins = 40)
      gpplai_amazonia <- plot(bin_gpplai, main="Amazonia GPP~LAI Bin",xlab='LAI',ylab='GPP', colramp=colors) 
      gppcica_amazonia <- plot(bin_gppcica, main="Amazonia GPP~CICA Bin",xlab='CiCa',ylab='GPP', colramp=colors) 
    }
  }
  list(gpplai_nw,gpplai_sw,gpplai_ec,gpplai_bs,gpplai_gs,gpplai_amazonia,gppcica_nw,gppcica_sw,gppcica_ec,gppcica_bs,gppcica_gs,gppcica_amazonia)
}
#par(mfrow=c(2,3))
all_plot_list<-hexbin_plot_fun(all_pixels_compare_data_zero)
#hexbin_plot_fun(all_pixels_compare_data)

gpplai_plotlist<-list(gpplai_nw,gpplai_sw,gpplai_ec,gpplai_bs,gpplai_gs,gpplai_amazonia)
gppcica_plotlist<-list(gppcica_nw,gppcica_sw,gppcica_ec,gppcica_bs,gppcica_gs,gppcica_amazonia)
grid.arrange(gpplai_plotlist,ncol=3)
grid.arrange(gppcica_plotlist,ncol=3)


######extra unused#####
hexbin_plot_fun_wol <- function(comp_data) {
  regions <- unique(comp_data$region_name)
  colors <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
  par(mfrow=c(2,3))
  for (i in regions) {
    if(i=='amazon_nw'){
      y <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$gpp
      x <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$lai
      z <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$cica
      bin_gpplai<-hexbin(x,y, xbins = 40)
      bin_gppcica<-hexbin(z,y, xbins = 40)
      gpplai_nw <- plot(bin_gpplai, main="North West Amazonia GPP~LAI Bin",xlab='LAI',ylab='GPP', colramp=colors, legend=F ) 
      gppcica_nw <- plot(bin_gppcica, main="North West Amazonia GPP~CICA Bin",xlab='CiCa',ylab='GPP', colramp=colors, legend=F ) 
    }
    else if (i=='amazon_sw') {
      y <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$gpp
      x <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$lai
      z <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$cica
      bin_gpplai<-hexbin(x,y, xbins = 40)
      bin_gppcica<-hexbin(z,y, xbins = 40)
      gpplai_sw <- plot(bin_gpplai, main="South West Amazonia GPP~LAI Bin",xlab='LAI',ylab='GPP', colramp=colors, legend=F ) 
      gppcica_sw <- plot(bin_gppcica, main="South West Amazonia GPP~CICA Bin",xlab='CiCa',ylab='GPP', colramp=colors, legend=F ) 
    }
    else if (i=='amazon_ec') {
      y <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$gpp
      x <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$lai
      z <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$cica
      bin_gpplai<-hexbin(x,y, xbins = 40)
      bin_gppcica<-hexbin(z,y, xbins = 40)
      gpplai_ec <- plot(bin_gpplai, main="East central Amazonia GPP~LAI Bin",xlab='LAI',ylab='GPP', colramp=colors, legend=F ) 
      gppcica_ec <- plot(bin_gppcica, main="East central Amazonia GPP~CICA Bin",xlab='CiCa',ylab='GPP', colramp=colors, legend=F ) 
    }
    else if (i=='amazon_bs') {
      y <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$gpp
      x <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$lai
      z <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$cica
      bin_gpplai<-hexbin(x,y, xbins = 40)
      bin_gppcica<-hexbin(z,y, xbins = 40)
      gpplai_bs <- plot(bin_gpplai, main="Brazil Shield GPP~LAI Bin",xlab='LAI',ylab='GPP', colramp=colors, legend=F ) 
      gppcica_bs <- plot(bin_gppcica, main="Brazil Shield GPP~CICA Bin",xlab='CiCa',ylab='GPP', colramp=colors, legend=F ) 
    }
    else if (i=='amazon_gs') {
      y <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$gpp
      x <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$lai
      z <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$cica
      bin_gpplai<-hexbin(x,y, xbins = 40)
      bin_gppcica<-hexbin(z,y, xbins = 40)
      gpplai_gs <- plot(bin_gpplai, main="Guyana Shield GPP~LAI Bin",xlab='LAI',ylab='GPP', colramp=colors, legend=F ) 
      gppcica_gs <- plot(bin_gppcica, main="Guyana Shield GPP~CICA Bin",xlab='CiCa',ylab='GPP', colramp=colors, legend=F ) 
    }
    else{
      y <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$gpp
      x <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$lai
      z <- all_pixels_compare_data[all_pixels_compare_data$region_name==i,]$cica
      bin_gpplai<-hexbin(x,y, xbins = 40)
      bin_gppcica<-hexbin(z,y, xbins = 40)
      gpplai_amazonia <- plot(bin_gpplai, main="Amazonia GPP~LAI Bin",xlab='LAI',ylab='GPP', colramp=colors, legend=F ) 
      gppcica_amazonia <- plot(bin_gppcica, main="Amazonia GPP~CICA Bin",xlab='CiCa',ylab='GPP', colramp=colors, legend=F ) 
    }
  }
  layout(matrix(c(1, 2, 3, 4,5,6), 2, 3, byrow = TRUE))
  gpplai_nw
  gpplai_sw
  gpplai_ec
  gpplai_bs
  gpplai_gs
  gpplai_amazonia 
  gppcica_nw
  gppcica_sw 
  gppcica_ec 
  gppcica_bs
  gppcica_gs
  gppcica_amazonia
}