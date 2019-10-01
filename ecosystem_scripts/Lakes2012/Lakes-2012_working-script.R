source("./functions/load_packages_function.R")

# run local_data_load() once to download all the data locally to './raw-data/' folder
# which is ignored for pushing purposes so csvs are not kept remotely
# just locally
# source("./functions/local_data_download_function.R")
# source("./ecosystem_functions/lakes_local_download_function.R")

####
 nars_data = read_csv("data/nars_data_table.csv");nars_data$filename <- basename(nars_data$data_link)
# #####
# #check the data surveys available
# nars_data %>%
#   pull(Survey) %>% unique()

### Lakes data
#check the indicator data available
nars_data %>%
  filter(Survey == "Lakes 2012") %>%
  pull(Indicator) %>% unique()

#check Water Chemistry available
lakes2012 <- nars_data %>%
  filter(Survey == "Lakes 2012", 
         Indicator == "Water Chemistry") %>%
  pull(filename) # 2 files

#first water chemistry file
lakes2012_wc1 <- nars_data %>%
  filter(Survey == "Lakes 2012", 
         Indicator == "Water Chemistry") %>%
  pull(filename) %>% '[['(1) %>% #This is site data
  paste0("raw-data/",.) %>% read_csv() 

# second water chemistry file
lakes2012_wc2 <- nars_data %>%
  filter(Survey == "Lakes 2012", 
         Indicator == "Water Chemistry") %>%
  pull(filename) %>% '[['(2) %>% #This is actual chemistry data
  paste0("raw-data/",.) %>% read_csv()

# site information file
lakes2012_site = nars_data %>%
  filter(Survey == "Lakes 2012",
         Indicator == "Site Information") %>%
  pull(filename) %>% paste0("raw-data/",.) %>% read_csv()

#### Working with nlcd data to pull lat-longs #####
## modified from Kelly's nlcd_feddata.Rmd ##
# read in lake shapefiles
basin_shapes = read_sf(dsn = "./raw-data", layer = "Lake_Basins") %>%
  st_transform(crs = 4269)#convert to NAD83 projection
basin_shapes_big = basin_shapes %>% filter(AreaSqKM >= 100000)

basin_shapes_small = basin_shapes %>% filter(AreaSqKM <= 100000)
#plot(basin_shapes$geometry)

#subset a few basins to try pulling land cover data
basin_shapes_ex<- basin_shapes_small[1:10,]

leaflet(basin_shapes_ex) %>% #view the subset of 
  addProviderTiles(providers$Esri.WorldImagery) %>%
  addPolygons()

#source the lake nlcd function to extract landcover of watersheds
source("./ecosystem_scripts/Lakes2012/functions/lake_nlcd_function.R")

#run on subset of lakes basin_shapes_ex <- basin_shapes[1:10,]
# works but takes a while ~9 secs per lake.
# ten_lakes = list()
# #debugonce(get_nlcd_percents)
# tic();for(i in seq_along(basin_shapes_ex[[1]]))ten_lakes[[i]] <- get_nlcd_percents(basin_shapes_ex[i,]);toc()
# names(ten_lakes) <- basin_shapes_ex$NLA12_ID
# 
# ten_lakes[[1]]#View the first lake example
# sum(ten_lakes[[1]]$percent_cover)#check adds to one

#full run
#need to break the list into parts because large vectors on some (~900MB!!)
#1
basin_shapes1 = basin_shapes[c(1:466),]
lakes_landuse1 = list()
tic();for(i in seq_along(basin_shapes[[1]]))lakes_landuse1[[i]] <- get_nlcd_percents(basin_shapes[i,]);toc()
names(lakes_landuse1) <- basin_shapes1$NLA12_ID
saveRDS(lakes_landuse1, file = "./ecosystem_scripts/Lakes2012/01_lakes_landuse_list.rds")
lakes_landuse1 = readRDS(file = "./ecosystem_scripts/Lakes2012/01_lakes_landuse_list.rds")
###
basin_shapes467 = basin_shapes[467,]
lakes_landuse467 = list()
tic();for(i in seq_along(basin_shapes467[[1]]))lakes_landuse467[[i]] <- get_nlcd_percents(basin_shapes467[i,]);toc()
names(lakes_landuse467) <- basin_shapes467$NLA12_ID
saveRDS(lakes_landuse467, file = "./ecosystem_scripts/Lakes2012/467_lakes_landuse_list.rds")
####
#2
basin_shapes2 = basin_shapes[c(468:700),]
lakes_landuse2 = list()
tic();for(i in seq_along(basin_shapes2[[1]]))lakes_landuse2[[i]] <- get_nlcd_percents(basin_shapes2[i,]);toc()
names(lakes_landuse2) <- basin_shapes2$NLA12_ID
saveRDS(lakes_landuse2, file = "./ecosystem_scripts/Lakes2012/02_lakes_landuse_list.rds")
#3
basin_shapes3 = basin_shapes[c(701:1000),]
lakes_landuse3 = list()
tic();for(i in seq_along(basin_shapes3[[1]]))lakes_landuse3[[i]] <- get_nlcd_percents(basin_shapes3[i,]);toc()
names(lakes_landuse3) <- basin_shapes3$NLA12_ID
saveRDS(lakes_landuse3, file = "./ecosystem_scripts/Lakes2012/03_lakes_landuse_list.rds")
#4
basin_shapes4 = basin_shapes[c(1001:1500),]
lakes_landuse4 = list()
tic();for(i in seq_along(basin_shapes4[[1]]))lakes_landuse4[[i]] <- get_nlcd_percents(basin_shapes4[i,]);toc()
names(lakes_landuse4) <- basin_shapes4$NLA12_ID
saveRDS(lakes_landuse4, file = "./ecosystem_scripts/Lakes2012/04_lakes_landuse_list.rds")
#5
basin_shapes5 = basin_shapes[c(1501:1714),]
lakes_landuse5 = list()
tic();for(i in seq_along(basin_shapes5[[1]]))lakes_landuse5[[i]] <- get_nlcd_percents(basin_shapes5[i,]);toc()
names(lakes_landuse5) <- basin_shapes5$NLA12_ID
saveRDS(lakes_landuse5, file = "./ecosystem_scripts/Lakes2012/05_lakes_landuse_list.rds")

#saveRDS(lakes_landuse, file = "./ecosystem_scripts/Lakes2012/lake_landuse_list.rds")
#attempt to speed it with lapply functions.
ten_lakes = list()
debugonce(get_nlcd_percents)
tic();ten_lakes = map_dfr(basin_shapes_ex, get_nlcd_percents);toc()

##need to work out why getting error "Error in st_sfc(x, crs = attr(x, "proj4string")) : 
#is.numeric(crs) || is.character(crs) || inherits(crs, "crs") is not TRUE"
as(basin_shapes_ex[1,], "Spatial")
as(list(basin_shapes_ex[1,]), "Spatial")
#the projection is NAD83 EPSG code (for crs = ) == 4269
#WSG 84 = 4326

lakes_loc <- lakes2012_site %>%
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

lakes_loc_ex <- lakes2012_site[1:10,]  %>%
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
####
  ######## Code for other ecosystems #######

nars_data %>% 
  filter(Survey == "Wetlands 2011") %>%
  pull(Indicator) %>% unique()


w2011_wc <- nars_data %>% 
  filter(Survey == "Wetlands 2011",
         Indicator == "Water Chemistry") %>%
  pull(filename) %>% paste0("data/", .) %>%
  read_csv()

w2011 <- nars_data %>%
  filter(Survey == "Wetlands 2011",
         Indicator == "All Indicators") %>%
  pull(filename) %>% paste0("data/", .) %>%
  read_csv()
names(w2011_wc)

w2011_carbon <- nars_data %>% 
  filter(Survey == "Wetlands 2011",
         Indicator == "Soil Chemistry") %>%
  pull(filename) %>% paste0("data/", .) %>%
  read_csv()

rivers2008_wc <- nars_data %>% 
  filter(Survey == "Rivers and Streams 2008-2009",
         Indicator == "Water Chemistry") %>%
  pull(filename)  %>% '[[' (2)# because there are 2 files.. 
  paste0("raw-data/", .) %>%
  read_csv()
  
nars_data %>%
  filter(Survey == "Rivers and Streams 2008-2009",
         Indicator == "Benthic Macroinvertebrates") %>%
  pull(filename) %>% unique()
  
rivers2008_benth.inv <- nars_data %>%
  filter(Survey == "Rivers and Streams 2008-2009",
         Indicator == "Benthic Macroinvertebrates") %>%
  pull(filename) %>% '[[' (5) %>%
  paste0("raw-data/",.) %>%
  read.csv
  
########### End other system code ##################