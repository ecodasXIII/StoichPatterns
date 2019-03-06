#PRISM example for daily temperature data, given the Streams 2008 dates and lat/long points
library(prism)
library(raster)
library(zoo)

#Get site locations and sample dates from a specified dataset
SiteData <- read.csv("data/assessed_ncca2010_siteinfo.revised.06212016.csv", header = TRUE)

#may have to change date "format" depending on how dates are formatted in csv
SiteData$date <- as.Date(SiteData$DATE_COL, format = "%d-%b-%y") 
SiteData$lon <- SiteData$ALON_DD
SiteData$lat <- SiteData$ALAT_DD

SiteData=SiteData[!is.na(SiteData$date), ]

#Set a location where all the PRISM files will be saved
options(prism.path = "prism/PRISM_temp_daily")

#This function downloads all the files. There is a single file per date for the entire
#continental US
#Change "type" depending on parameter you want (mean, max, min temp or precip)
#get_prism_monthlys is the option for monthly data
#get_prism_annual is the option for annual data

get_prism_dailys(type ="tmean", dates = unique(SiteData$date), keepZip = FALSE)

#Loop through all the date points we want to get temperature data for
for (i in 1:nrow(SiteData)){
  
  #These lines get the right file name for the sample date
  dt <- gsub("-", "", SiteData$date[i])
  file_nm <- grep(dt, ls_prism_data()[,1], value = T)
  file_nm <- grep("tmean", file_nm, value = T)
  
  #This pulls the data from the file for the lon/lat point you supply
  data <- prism_slice(c(SiteData$lon[i], SiteData$lat[i]), file_nm)
  
  #It returns a list so you have to specify the datapoint you want
  #Adds a column to the original dataset with the mean daily temperature
  SiteData$tmean_day[i] <- data$data$data[1]
}

