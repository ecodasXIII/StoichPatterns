#load functions script
load_packages = function() {
if(!require("pacman")) install.packages("pacman")
library(pacman)
package.list = c("FedData","raster","sf","leaflet",
                 "readr","tidyverse","tictoc","fs",
                 "furrr")
pacman::p_load(char = package.list, install = T)
pacman::p_load_gh('jimjunker1/junkR')
rm('package.list')
theme_mod <<- theme_bw() %+replace% theme(panel.grid = element_blank())
}
load_packages()
