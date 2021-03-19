get_ecomap<-function(){
  
  ### read in Ecoregion map
  ecomap1<-readOGR(paste(map.in,"na_cec_eco_l1",sep=""),
                   layer=c("NA_CEC_Eco_Level1"),verbose=F)
  
  proj4string(ecomap1)
  ecomap1 <- spTransform(ecomap1, CRS("+proj=longlat +datum=WGS84"))
  
  # Next the shapefile has to be converted to a dataframe for use in ggplot2
  ecomap1_df <- fortify(ecomap1,region="NA_L1NAME")
  ecomap1_df$code<-as.numeric(as.factor(ecomap1_df$id))
  
  # Now the shapefile can be plotted as either a geom_path or a geom_polygon.
  # Paths handle clipping better. Polygons can be filled.
  # You need the aesthetics long, lat, and group.
  ecomap_base <- ggplot() +
    geom_polygon(data = ecomap1_df, 
                 aes(x = long, y = lat, group = group, fill=id),
                 color = "darkgrey",
                 size = 0.5) +
    #scale_colour_grey() +
    theme_classic() +
    scale_fill_hue(c=20, l=90) +
    theme(legend.title=element_blank()) + 
    theme(legend.text=element_text(size=rel(0.5))) +
    xlab("Longitude") +
    ylab("Latitude") +
    coord_quickmap()
  #scale_colour_grey()
  
  ecomap_base2 <- ggplot() +
    geom_polygon(data = ecomap1_df, 
                 aes(x = long, y = lat, group = group, fill=id),
                 color = "darkgrey",
                 size = 0.5) +
    #scale_colour_grey() +
    theme_dark() +
    scale_fill_hue(c=20, l=90) +
    theme(legend.title=element_blank()) + 
    theme(legend.text=element_text(size=rel(0.5))) +
    xlab("Longitude") +
    ylab("Latitude") +
    coord_quickmap()
  #scale_colour_grey()
  
 return(list=ls(ecomap_base,ecomap_base2)) 
}
