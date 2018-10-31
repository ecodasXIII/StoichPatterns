library(tidyverse)
#source("./data-download.Rmd")

nars_data = read_csv("data/nars_data_table.csv")

nars_data %>% 
  filter(Survey == "Wetlands 2011") %>%
  pull(Indicator) %>% unique()
nars_data$filename <- basename(nars_data$data_link)


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

lakes
