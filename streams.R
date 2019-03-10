##### Streams and Rivers 2008/9 ####
pacman::p_load(FedData, fs, furrr, leaflet, raster, tictoc)
library(sf)
library(tidyverse)

# Full set of NARS data files
nars_data <- read_csv("data/nars_data_table.csv") %>%
  mutate(filename = basename(.$data_link))

##### Water chemistry --> CNP ####
chem <- nars_data %>%
  filter(Survey %in% c("Rivers and Streams 2008-2009"), 
         Indicator %in% c("Water Chemistry"),
         str_detect(Data, 'Water Chemistry - Data')) %>% # Select for the raw data file
  pull(filename) %>% paste0("data/", .) %>%
  read_csv()
  
cnp <- chem %>% 
  select(YEAR, UID, SITE_ID, DOC, DOC_ALERT, NTL, NTL_ALERT, PTL, PTL_ALERT) %>% 
  mutate(CN = (DOC / (NTL/1000)) * (14.01/12.01), # Convert NTL to mg/L
         CP = (DOC / (PTL/1000)) * (30.97/12.01), # Convert PTL to mg/L
         NP = ((NTL/1000) / (PTL/1000)) * (30.97/14.01))  # Convert NTL, PTL to mg/L

# write_csv(chem, './data/streams/cnp_chem.csv')

#### Site info w/ PRISM daily/annual temp/precip ####
siteinfo <- read_csv('./data/streams/site_info_PRISM.csv')

#### Landscape data ####
landscape_all <- nars_data %>% 
  filter(Survey == "Rivers and Streams 2008-2009",
         Indicator == "Landscape Data",
         str_detect(Data, ' Landscape Metrics')) %>% # Other file is shapefiles
  pull(filename) %>% paste0("data/", .) %>%
  read_csv()

# Select focal landscape variables (temperatures, land cover)
landscape <- landscape_all %>%
  select(UID, SITE_ID, TMAX_ANN:TMIN_JUL, PCT_AG:PCT_SHRUB_GRASS, NHDWAT_PCT_IMPERV)

#### LANDSCAPE FROM SHAPEFILES ####
# unzip shapefile & data for lake basins.
sr_landscape2 <- unzip(zipfile = './data/nrsa0809_watersheds.zip', 
        exdir = "./stream_watersheds_data")

# read in lake shapefiles
basin_shapes = read_sf(dsn = "./stream_watersheds_data", layer = "NRSA0809_Watersheds") %>%
  st_transform(crs = 4269) # convert to NAD83 projection

#subset a few basins to try pulling land cover data
basin_shapes_ex <- basin_shapes[1:10,]

#leaflet(basin_shapes_ex) %>% #view the subset of 
#  addProviderTiles(providers$Esri.WorldImagery) %>%
#  addPolygons()

#source the stream nlcd function to extract landcover of watersheds
source("./stream_nlcd_function.R")

#run on subset of basin_shapes_ex <- basin_shapes[1:10,]

#ten_streams = list()

#tic();for(i in seq_along(basin_shapes_ex[[1]]))
#  ten_streams[[i]] <- get_nlcd_percents(basin_shapes_ex[i,]);toc()
#names(ten_streams) <- basin_shapes_ex$SITE_ID

#ten_streams[[1]] #View the first stream example

## Run all watersheds ####
all_streams = list()

for(i in seq_along(basin_shapes[[1]]))
  all_streams[[i]] <- get_nlcd_percents(basin_shapes[i,])
names(all_streams) <- basin_shapes$SITE_ID