
#attempt to speed it with lapply functions.

ten_lakes = list()
debugonce(get_nlcd_percents)
tic();ten_lakes = apply(basin_shapes_ex, 1, FUN = get_nlcd_percents);toc()

##need to work out why getting error "Error in st_sfc(x, crs = attr(x, "proj4string")) : 
#is.numeric(crs) || is.character(crs) || inherits(crs, "crs") is not TRUE"
as(basin_shapes_ex[1,], "Spatial")
as(list(basin_shapes_ex[1,]), "Spatial")
#the projection is NAD83 EPSG code (for crs = ) == 4269
#WSG 84 = 4326

lakes_loc <- siteinfo %>%
  st_as_sf(coords = c("LON_DD83","LAT_DD83"), crs = 4269)

#look at the first few sites
lakes_loc[1:5,] %>%
  st_transform(crs = 4269) %>%
  leaflet() %>%
  addTiles() %>%
  addMarkers()

#plot them all
plot(lakes_loc$geometry)

# Make buffer around first point

pt1_buffer <- lakes_loc %>%
  dplyr::slice(1) %>%
  st_buffer(0.01) #units in decimal degree. Should ~ 1000m

# plot buffer
leaflet(pt1_buffer) %>%
  addProviderTiles(providers$Esri.WorldImagery) %>%
  addPolygons()

# use buffer to get nlcd function

nlcd_nars  <- get_nlcd(as(pt1_buffer, "Spatial"), 
                       label = "first lake site",
                       dataset = "landcover",
                       year = "2011",
                       force.redo = TRUE)### had an error message pop up saying corrupt or incomplete .zip file

nlcd_nars %>% plot()##seems to be a pretty painting

leaflet(pt1_buffer) %>%
  addProviderTiles(providers$Esri.WorldImagery) %>%
  addPolygons() %>%
  addRasterImage(nlcd_nars)

#Mask raster using buffer object after transforming to same CRS

pt1_buffer_prj <- st_transform(pt1_buffer, crs = proj4string(nlcd_nars))
nlcd_nars_mask <- raster::mask(nlcd_nars, as(pt1_buffer_prj, "Spatial"))

plot(nlcd_nars_mask)

#tabulate number of cells in each type and merge with legend to see land cover types

lc_legend_df <- readr::read_csv("nlcd_legend_2011.csv")

raster::freq(nlcd_nars_mask) %>%
  as.data.frame() %>%
  dplyr::filter(!is.na(value)) %>%
  dplyr::left_join(lc_legend_df, by = c('value' = 'Class')) %>%
  dplyr::mutate(percent_cover = count/sum(count)) %>%
  dplyr::select(class_name, percent_cover)

##this seems to work well based on a circular buffer
##need to think about how to apply this to lakes and streams??

## Define function to do all of the above --get percent of land cover types
# within polygon. Returns polygon with columns class_name and percent cover

# this works for a single one now iterate across a subset of rows in lakes_loc
#automate the buffer getting for a subset of lakes

lakes_loc_ex <- siteinfo[1:10,]  %>%
  st_as_sf(coords = c("LON_DD83","LAT_DD83"), crs = 4269)

lakes_buffer <- lakes_loc_ex %>%
  st_buffer(0.01)

ten_lakes = list()
tic();for(i in seq_along(lakes_buffer[[1]]))ten_lakes[[i]] <- get_nlcd_percents(lakes_buffer[i,]);toc()
# for loop works, but likely slow. Prefer to use apply-group. Working on that
# got help from here: https://stackoverflow.com/questions/48505551/use-apply-with-a-simple-features-sf-function


debugonce(get_nlcd_percents)
ten_lakes = lapply(st_geometry(lakes_buffer), FUN = get_nlcd_percents)

##need to work out why getting error "Error in st_sfc(x, crs = attr(x, "proj4string")) : 
#is.numeric(crs) || is.character(crs) || inherits(crs, "crs") is not TRUE"

crs(lakes_buffer)#;crs(lakes_buffer)<-CRS("+proj=longlat +epsg=4269")
as(lakes_buffer, 'Spatial') 
st_geometry(lakes_buffer)
as(st_geometry(lakes_buffer),'Spatial')
crs(as(st_geometry(lakes_buffer),'Spatial')) <-CRS("+proj=longlat +epsg=4269")