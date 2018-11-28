library(tidyverse)
# run once to download all the data locally to './raw-data/' folder
# which is ignored for pushing purposes so csvs are not kept remotely
# just locally

#source("./data-download.Rmd")
####
nars_data = read_csv("data/nars_data_table.csv");nars_data$filename <- basename(nars_data$data_link)
#####
#check the data surveys available
nars_data %>%
  pull(Survey) %>% unique()

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

#This unzips shapefile & data for lake basins.
lakes2012_landscape= nars_data %>% #this is GIS file (Shapes, etc.)
  filter(Survey == "Lakes 2012",
         Indicator == "Landscape Data") %>%
  pull(filename) %>% '[[' (2) %>% #both 1 & 2 have sites and lake shapes
  paste0("raw-data/",.) %>% unzip() 

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
  pull(filename) %>% 
  head(1) %>% # because there are 2 files.. 
  paste0("data/", .) %>%
  read_csv()
########################################