#PRISM example for daily temperature data, given the Streams 2008 dates and lat/long points
library(prism)
library(raster)
library(zoo)

#Get site locations and sample dates from a specified dataset
streams <- read.csv("data/siteinfo_0.csv", header = TRUE)

#may have to change date "format" depending on how dates are formatted in csv
streams$date <- as.Date(streams$DATE_COL, format = "%d-%b-%y") 
streams$lon <- streams$LON_DD83
streams$lat <- streams$LAT_DD83

#Set a location where all the PRISM files will be saved
options(prism.path = "prism/PRISM_temp_daily")

#This function downloads all the files. There is a single file per date for the entire
#continental US
#Change "type" depending on parameter you want (mean, max, min temp or precip)
#get_prism_monthlys is the option for monthly data
#get_prism_annual is the option for annual data

get_prism_dailys(type ="tmean", dates = unique(streams$date), keepZip = FALSE)

#Loop through all the date points we want to get temperature data for
for (i in 1:nrow(streams)){
  
  #These lines get the right file name for the sample date
  dt <- gsub("-", "", streams$date[i])
  file_nm <- grep(dt, ls_prism_data()[,1], value = T)
  file_nm <- grep("tmean", file_nm, value = T)
  
  #This pulls the data from the file for the lon/lat point you supply
  data <- prism_slice(c(streams$lon[i], streams$lat[i]), file_nm)
  
  #It returns a list so you have to specify the datapoint you want
  #Adds a column to the original dataset with the mean daily temperature
  streams$tmean_day[i] <- data$data$data[1]
}

