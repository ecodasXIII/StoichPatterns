#download lake shapefiles locally
# they are downloaded into "raw-data" folder which is ignored in
# .gitignore files
lake_local_download = function(){

  
  # #This unzips shapefile & data for lake basins.
  lakes2012_landscape1= nars_data %>% #this is GIS file (Shapes, etc.)
    filter(Survey == "Lakes 2012",
           Indicator == "Landscape Data") %>%
    pull(filename) %>% '[[' (1) %>% #this is sample points and polys
    unzip(exdir = './raw-data')
  lakes2012_landscape2= nars_data %>% #this is GIS file (Shapes, etc.)
    filter(Survey == "Lakes 2012",
           Indicator == "Landscape Data") %>%
    pull(filename) %>% '[[' (2) %>% #this is lake basin poly
    unzip(exdir = "./raw-data")
}