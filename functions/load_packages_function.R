#load functions script
load_packages = function() {
if(!require("pacman")) install.packages("pacman")
library(pacman)
package.list = c("FedData","raster","sf","leaflet",
                 "readr","tidyverse","tictoc","fs",
                 "furrr")
p_load(char = package.list, install = T)
rm('package.list')
}
load_packages()