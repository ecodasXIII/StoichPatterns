get_nlcd_percents <- function(aoi, nlcd_year = "2011"){
  aoi_sp <- as(aoi, "Spatial")
  #;crs(aoi_sp) <- CRS("+proj=longlat +epsg=4269")
  #print(crs(aoi_sp));print(crs(aoi))
  nlcd_nars <- FedData::get_nlcd(aoi_sp, label = "aoi",
                                 dataset = "landcover",
                                 year = nlcd_year,
                                 force.redo = TRUE,
                                 raw.dir = "./raw-data/RAW",
                                 extraction.dir = "./raw-data/EXTRACTIONS")
  aoi_prj <- st_transform(aoi, crs = proj4string(nlcd_nars))
  nlcd_nars_mask <- raster::mask(nlcd_nars, as(aoi_prj, "Spatial"))
  
  lc_legend_df <- readr::read_csv("nlcd_legend_2011.csv")
  
  cover <- raster::freq(nlcd_nars_mask) %>%
    as.data.frame() %>%
    dplyr::filter(!is.na(value)) %>%
    dplyr::left_join(lc_legend_df, by = c("value" = "Class")) %>%
    dplyr::mutate(percent_cover = count/sum(count)) %>%
    dplyr::select(class_name, percent_cover)
  
  return(cover)
}