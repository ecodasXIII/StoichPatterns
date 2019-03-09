##### Streams and Rivers 2008/9 ####
pacman::p_load(tidyverse)

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

#### PRISM daily/annual temp/precip ####
prism <- read_csv('./data/streams/site_info_PRISM.csv')

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
