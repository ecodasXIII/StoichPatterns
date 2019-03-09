#PRISM example for daily temperature data, given the Streams 2008 dates and lat/long points
pacman::p_load(prism, raster, tidyverse, zoo)

#Get site locations and sample dates from a specified dataset
SiteData <- read_csv("data/siteinfo_0.csv") %>% 
  mutate(date = as.Date(DATE_COL, format = "%d-%b-%y"),
         lon = LON_DD83, 
         lat = LAT_DD83) %>% 
  filter(!is.na(date))

# Subset the file 
sites <- SiteData %>% 
  select(YEAR, date, lat, lon, UID, SITE_ID, STATE:EPA_REG, WSAREA_NARS)

# This function downloads all the files- one file per date for continental US
# Change "type" depending on parameter you want (mean, max, min temp or precip)
# get_prism_monthlys is the option for monthly data
# get_prism_annual is the option for annual data

# Daily mean temperature ####
options(prism.path = "prism/streams/PRISM_temp_daily")
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
  sites$tmean_day[i] <- data$data$data[1]
}

# Daily precip ####
options(prism.path = "prism/streams/PRISM_precip_daily")
get_prism_dailys(type ="ppt", dates = unique(SiteData$date), keepZip = FALSE)

#Loop through all the date points we want to get temperature data for
for (i in 1:nrow(SiteData)){
  
  #These lines get the right file name for the sample date
  dt <- gsub("-", "", SiteData$date[i])
  file_nm <- grep(dt, ls_prism_data()[,1], value = T)
  file_nm <- grep("ppt", file_nm, value = T)
  
  #This pulls the data from the file for the lon/lat point you supply
  data <- prism_slice(c(SiteData$lon[i], SiteData$lat[i]), file_nm)
  
  #It returns a list so you have to specify the datapoint you want
  #Adds a column to the original dataset with the mean daily temperature
  sites$ppt_day[i] <- data$data$data[1]
}

# Annual mean temperature ####
options(prism.path = "prism/streams/PRISM_temp_yearly")
get_prism_annual(type ="tmean", years = unique(SiteData$YEAR), keepZip = FALSE)

#Loop through all the date points we want to get temperature data for
for (i in 1:nrow(SiteData)){
  
  #These lines get the right file name for the sample date
  dt <- gsub("-", "", SiteData$YEAR[i])
  file_nm <- grep(dt, ls_prism_data()[,1], value = T)
  file_nm <- grep("tmean", file_nm, value = T)
  
  #This pulls the data from the file for the lon/lat point you supply
  data <- prism_slice(c(SiteData$lon[i], SiteData$lat[i]), file_nm)
  
  #It returns a list so you have to specify the datapoint you want
  #Adds a column to the original dataset with the mean daily temperature
  sites$tmean_year[i] <- data$data$data[1]
}

# Annual mean precip ####
options(prism.path = "prism/streams/PRISM_precip_yearly")
get_prism_annual(type ="ppt", years = unique(SiteData$YEAR), keepZip = FALSE)

#Loop through all the date points we want to get temperature data for
for (i in 1:nrow(SiteData)){
  
  #These lines get the right file name for the sample date
  dt <- gsub("-", "", SiteData$YEAR[i])
  file_nm <- grep(dt, ls_prism_data()[,1], value = T)
  file_nm <- grep("ppt", file_nm, value = T)
  
  #This pulls the data from the file for the lon/lat point you supply
  data <- prism_slice(c(SiteData$lon[i], SiteData$lat[i]), file_nm)
  
  #It returns a list so you have to specify the datapoint you want
  #Adds a column to the original dataset with the mean daily temperature
  sites$ppt_year[i] <- data$data$data[1]
}

# Save file! 
write_csv(sites, './data/streams/site_info_PRISM.csv', append=FALSE)
