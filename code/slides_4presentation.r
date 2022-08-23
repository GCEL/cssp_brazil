
avitable_map <- brick('G://AGB/avitabile/Avitabile_AGB_Map/Avitabile_AGB_Map.tif')
plot(avitable_map)

amazonia_poly <-rasterToPolygons(biommort_10_16, dissolve=TRUE)
poly <- rasterToPolygons(masked2, dissolve=TRUE)
plot(biomass_amazon_mask_pol,add=T)

SA_extent <- extent(c(-83.75,-31.75,-56.75,13.75)) #use south america extent from INLAND SA
SA_r <- raster(SA_extent)
res(SA_r)<- 0.008333333 #use resolution to mirror other masks in ilamb/ could be 1x1
values(SA_r) <- 1
crs(SA_r) <- "+proj=longlat +datum=WGS84 +no_defs"

avitable_map_sa <- crop(avitable_map, SA_r)
plot(avitable_map_sa)
plot(biomass_amazon_mask_pol,add=T)
