#local data downloads 
local_data_download = function(){
  #from:
  #source("./data-download.Rmd")
  #create local folder
  ifelse(!dir.exists(file.path(getwd(), "/raw-data")), dir.create(file.path(getwd(), "/raw-data")), FALSE)
  
  #Read in the csv table created with `get_nars_links.R`. 
  
  nars_data <- read_csv("data/nars_data_table.csv")
  nars_data$filename <- basename(nars_data$data_link)
  
  #Download all the csv files in the data links column into a folder called "raw-data" (needs to exist first)
  map2(.x = nars_data$data_link,
       .y = file.path("raw-data", basename(nars_data$data_link)), 
       ~download.file(.x, .y))
}
local_data_download()