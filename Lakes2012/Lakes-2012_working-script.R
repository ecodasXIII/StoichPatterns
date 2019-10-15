source("./functions/load_packages_function.R")

# run local_data_load() once to download all the data locally to './raw-data/' folder
# which is ignored for pushing purposes so csvs are not kept remotely
# just locally
# source("./functions/local_data_download_function.R")
# source("./ecosystem_functions/lakes_local_download_function.R")

####
 nars_data = read_csv("data/nars_data_table.csv");nars_data$filename <- basename(nars_data$data_link)
# #####
# #check the data surveys available
# nars_data %>%
#   pull(Survey) %>% unique()

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

#### Working with nlcd data to pull lat-longs #####
## modified from Kelly's nlcd_feddata.Rmd ##
# read in lake shapefiles
basin_shapes = read_sf(dsn = "./raw-data", layer = "Lake_Basins") %>%
  st_transform(crs = 4269)#convert to NAD83 projection
basin_shapes_big = basin_shapes %>% filter(AreaSqKM >= 100000)

basin_shapes_small = basin_shapes %>% filter(AreaSqKM <= 100000)
#plot(basin_shapes$geometry)

#subset a few basins to try pulling land cover data
basin_shapes_ex<- basin_shapes_small[1:10,]

leaflet(basin_shapes_ex) %>% #view the subset of 
  addProviderTiles(providers$Esri.WorldImagery) %>%
  addPolygons()

#source the lake nlcd function to extract landcover of watersheds
source("./ecosystem_scripts/Lakes2012/functions/lake_nlcd_function.R")

#run on subset of lakes basin_shapes_ex <- basin_shapes[1:10,]
# works but takes a while ~9 secs per lake.
# ten_lakes = list()
# #debugonce(get_nlcd_percents)
# tic();for(i in seq_along(basin_shapes_ex[[1]]))ten_lakes[[i]] <- get_nlcd_percents(basin_shapes_ex[i,]);toc()
# names(ten_lakes) <- basin_shapes_ex$NLA12_ID
# 
# ten_lakes[[1]]#View the first lake example
# sum(ten_lakes[[1]]$percent_cover)#check adds to one

#full run
#need to break the list into parts because large vectors on some (~900MB!!)
#1
basin_shapes1 = basin_shapes[c(1:466),]
lakes_landuse1 = list()
tic();for(i in seq_along(basin_shapes[[1]]))lakes_landuse1[[i]] <- get_nlcd_percents(basin_shapes[i,]);toc()
names(lakes_landuse1) <- basin_shapes1$NLA12_ID
saveRDS(lakes_landuse1, file = "./ecosystem_scripts/Lakes2012/01_lakes_landuse_list.rds")
lakes_landuse1 = readRDS(file = "./ecosystem_scripts/Lakes2012/01_lakes_landuse_list.rds")
###
basin_shapes467 = basin_shapes[467,]
lakes_landuse467 = list()
tic();for(i in seq_along(basin_shapes467[[1]]))lakes_landuse467[[i]] <- get_nlcd_percents(basin_shapes467[i,]);toc()
names(lakes_landuse467) <- basin_shapes467$NLA12_ID
saveRDS(lakes_landuse467, file = "./ecosystem_scripts/Lakes2012/467_lakes_landuse_list.rds")
####
#2
basin_shapes2 = basin_shapes[c(468:700),]
lakes_landuse2 = list()
tic();for(i in seq_along(basin_shapes2[[1]]))lakes_landuse2[[i]] <- get_nlcd_percents(basin_shapes2[i,]);toc()
names(lakes_landuse2) <- basin_shapes2$NLA12_ID
saveRDS(lakes_landuse2, file = "./ecosystem_scripts/Lakes2012/02_lakes_landuse_list.rds")
#3
basin_shapes3 = basin_shapes[c(701:1000),]
lakes_landuse3 = list()
tic();for(i in seq_along(basin_shapes3[[1]]))lakes_landuse3[[i]] <- get_nlcd_percents(basin_shapes3[i,]);toc()
names(lakes_landuse3) <- basin_shapes3$NLA12_ID
saveRDS(lakes_landuse3, file = "./ecosystem_scripts/Lakes2012/03_lakes_landuse_list.rds")
#4
basin_shapes4 = basin_shapes[c(1001:1500),]
lakes_landuse4 = list()
tic();for(i in seq_along(basin_shapes4[[1]]))lakes_landuse4[[i]] <- get_nlcd_percents(basin_shapes4[i,]);toc()
names(lakes_landuse4) <- basin_shapes4$NLA12_ID
saveRDS(lakes_landuse4, file = "./ecosystem_scripts/Lakes2012/04_lakes_landuse_list.rds")
#5
basin_shapes5 = basin_shapes[c(1501:1714),]
lakes_landuse5 = list()
tic();for(i in seq_along(basin_shapes5[[1]]))lakes_landuse5[[i]] <- get_nlcd_percents(basin_shapes5[i,]);toc()
names(lakes_landuse5) <- basin_shapes5$NLA12_ID
saveRDS(lakes_landuse5, file = "./ecosystem_scripts/Lakes2012/05_lakes_landuse_list.rds")


#######  CODE FOR RMARKDOWN OUTPUT #########

#create new CNP dataframe
# lakevars <- c("UID", "SITE_ID", "AREA_HA", "DATE_COL", "TEMPERATURE", "DOC_RESULT", "DOC_QA_FLAG", "NTL_RESULT", "NTL_QA_FLAG", "PTL_RESULT", "PTL_QA_FLAG")
lakevars <- c("UID", "SITE_ID", "AREA_HA", "DATE_COL", "DOC_RESULT", "DOC_QA_FLAG", "NTL_RESULT", "NTL_QA_FLAG", "PTL_RESULT", "PTL_QA_FLAG")


lakes2012_meta <- dplyr::select(lakes2012_wc1, one_of(lakevars)) %>%
  unique()
lakes2012_data <- dplyr::select(lakes2012_wc2, one_of(lakevars))
lakes2012_site <- dplyr::select(lakes2012_site, one_of(lakevars))

lakes2012_CNP <- dplyr::inner_join(lakes2012_meta, lakes2012_data) %>%
  inner_join(lakes2012_site %>% select(-DATE_COL)) %>%
  mutate(PTL_RESULT = PTL_RESULT/1000,#get to mg/L
         C_mol = (DOC_RESULT/12.011)*1000,#convert to umol
         N_mol = (NTL_RESULT/14.007)*1000,#convert to umol
         P_mol = (PTL_RESULT/30.974)*1000,#convert to umol
         CN_mol = (DOC_RESULT/NTL_RESULT)*(14.007/12.011),
         CP_mol = (DOC_RESULT/PTL_RESULT)*(30.974/12.011),
         NP_mol = (NTL_RESULT/PTL_RESULT)*(30.974/14.007))

# land cover work with landuse lists 1-5

#read in the lists of land use
land_use_temp = list.files(path = "./Lakes2012/", pattern = "*landuse_list.rds", full.names = TRUE)
land_use_files = lapply(land_use_temp, readRDS)#make list of landuse lists

land_use_dfs_list = map(land_use_files, function(x){
  x_names = names(x)
  x_df_list = map(x, ~..1 %>% spread(key = class_name, value = percent_cover))
  x_dfs = bind_rows(x_df_list)
  x_df_full = data.frame(SITE_ID = x_names, x_dfs, check.names = FALSE)
  colnames(x_df_full) = gsub(" ","_",colnames(x_df_full))
  return(x_df_full)
})

land_use_df = land_use_dfs_list %>% 
              bind_rows() %>%
              gather(key = class_name, value = percent_cover, `Deciduous_Forest`:`Perennial_Ice/Snow`)

land_use_df_wide = land_use_df %>% 
  group_by(SITE_ID) %>% summarise(total_use = sum(percent_cover, na.rm = TRUE))

#rename cover classes to aggregrated categories
#old names
old_cover_classes = list(c('Pasture/Hay','Cultivated_Crops'),
                         c('Developed_High_Intensity','Developed,_Low_Intensity',
                           'Developed,_Medium_Intensity','Developed,_Open_Space'),
                         c('Deciduous_Forest','Evergreen_Forest','Mixed_Forest'),
                         c('Open_Water','Woody_Wetlands','Emergent_Herbaceous_Wetlands'),
                         'Barren_Land_(Rock/Sand/Clay)',
                         'Grassland/Herbaceous',
                         'Shrub/Scrub',
                         'Perennial_Ice/Snow',
                         '<NA>')
#aggregated names
aggregate_cover_class = list('Agriculture',
                             'Developed',
                             'Forested',
                             'WetlandWater',
                             'Barren',
                             'Grassland',
                             'Shrubscrub',
                             'Ice_Snow',
                             "Unknown")

#keyval list of changes
keyval = setNames(rep(aggregate_cover_class, lengths(old_cover_classes)), unlist(old_cover_classes))

#make changes to class names
land_use_df = land_use_df %>%
              mutate(class_name = recode(class_name, !!!keyval)) %>%
              group_by(SITE_ID, class_name) %>%
              summarise(percent_cover = sum(percent_cover)) %>% ungroup()

# Quick check 
unique(levels(as.factor(land_use_df$class_name)))

lakes2012_full = inner_join(lakes2012_CNP, land_use_df %>% spread(key = class_name, value = percent_cover)) %>%
  mutate_at(vars(Agriculture:WetlandWater), replace_na, 0)


#Look at the distribution of watershed sizes
log_breaks = c(log10(0.5), log10(1), log10(5), log10(10), log10(50), log10(100), log10(500), log10(1000), log10(5000), log10(10000), log10(50000),log10(100000))
raw_breaks = round(10^(log_breaks),3)

Area_plot = ggplot(lakes2012_full, aes(x = log10(AREA_HA))) + geom_histogram() +
  scale_x_continuous(name = "Lake Area (Hectares)", breaks = log_breaks, labels = raw_breaks)+
  theme(panel.grid = element_blank(), axis.title.y = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1,margin = margin(t = 3)))

######## Create stoic isocline dataframe for plotting
NP_0.1 = data.frame(Stoic = "NP_0.1", Ppool = seq(0.0001,10000, by = 0.1), Npool = seq(0.0001,10000, by = 0.1)*0.1)
NP_0.2 = data.frame(Stoic = "NP_0.2", Ppool = seq(0.0001,10000, by = 0.1), Npool = seq(0.0001,10000, by = 0.1)*0.2)
NP_0.3 = data.frame(Stoic = "NP_0.3", Ppool = seq(0.0001,10000, by = 0.1), Npool = seq(0.0001,10000, by = 0.1)*0.3)
NP_0.4 = data.frame(Stoic = "NP_0.4", Ppool = seq(0.0001,10000, by = 0.1), Npool = seq(0.0001,10000, by = 0.1)*0.4)
NP_0.5 = data.frame(Stoic = "NP_0.5", Ppool = seq(0.0001,10000, by = 0.1), Npool = seq(0.0001,10000, by = 0.1)*0.5)
NP_0.6 = data.frame(Stoic = "NP_0.6", Ppool = seq(0.0001,10000, by = 0.1), Npool = seq(0.0001,10000, by = 0.1)*0.6)
NP_0.7 = data.frame(Stoic = "NP_0.7", Ppool = seq(0.0001,10000, by = 0.1), Npool = seq(0.0001,10000, by = 0.1)*0.7)
NP_0.8 = data.frame(Stoic = "NP_0.8", Ppool = seq(0.0001,10000, by = 0.1), Npool = seq(0.0001,10000, by = 0.1)*0.8)
NP_0.9 = data.frame(Stoic = "NP_0.9", Ppool = seq(0.0001,10000, by = 0.1), Npool = seq(0.0001,10000, by = 0.1)*0.9)
NP_1 = data.frame(Stoic = "NP_1", Ppool = seq(0.0001,10000, by = .1), Npool = seq(0.0001,10000, by = .1))
NP_2 = data.frame(Stoic = "NP_2", Ppool = seq(0.0001,10000, by = .1), Npool = seq(0.0001,10000, by = .1)*2)
NP_3 = data.frame(Stoic = "NP_3", Ppool = seq(0.0001,10000, by = .1), Npool = seq(0.0001,10000, by = .1)*3)
NP_4 = data.frame(Stoic = "NP_4", Ppool = seq(0.0001,10000, by = .1), Npool = seq(0.0001,10000, by = .1)*4)
NP_5 = data.frame(Stoic = "NP_5", Ppool = seq(0.0001,10000, by = .1), Npool = seq(0.0001,10000, by = .1)*5)
NP_6 = data.frame(Stoic = "NP_6", Ppool = seq(0.0001,10000, by = .1), Npool = seq(0.0001,10000, by = .1)*6)
NP_7 = data.frame(Stoic = "NP_7", Ppool = seq(0.0001,10000, by = .1), Npool = seq(0.0001,10000, by = .1)*7)
NP_8 = data.frame(Stoic = "NP_8", Ppool = seq(0.0001,10000, by = .1), Npool = seq(0.0001,10000, by = .1)*8)
NP_9 = data.frame(Stoic = "NP_9", Ppool = seq(0.0001,10000, by = .1), Npool = seq(0.0001,10000, by = .1)*9)
NP_10 = data.frame(Stoic = "NP_10", Ppool = seq(0.0001,10000, by = .1), Npool = seq(0.0001,10000, by = .1)*10)
NP_20 = data.frame(Stoic = "NP_20", Ppool = seq(0.0001,10000, by = .1), Npool = seq(0.0001,10000, by = .1)*20)
NP_30 = data.frame(Stoic = "NP_30", Ppool = seq(0.0001,10000, by = .1), Npool = seq(0.0001,10000, by = .1)*30)
NP_40 = data.frame(Stoic = "NP_40", Ppool = seq(0.0001,10000, by = .1), Npool = seq(0.0001,10000, by = .1)*40)
NP_50 = data.frame(Stoic = "NP_50", Ppool = seq(0.0001,10000, by = 0.1), Npool = seq(0.0001,10000, by = 0.1)*50)
NP_60 = data.frame(Stoic = "NP_60", Ppool = seq(0.0001,10000, by = 0.1), Npool = seq(0.0001,10000, by = 0.1)*60)
NP_70 = data.frame(Stoic = "NP_70", Ppool = seq(0.0001,10000, by = 0.1), Npool = seq(0.0001,10000, by = 0.1)*70)
NP_80 = data.frame(Stoic = "NP_80", Ppool = seq(0.0001,10000, by = 0.1), Npool = seq(0.0001,10000, by = 0.1)*80)
NP_90 = data.frame(Stoic = "NP_90", Ppool = seq(0.0001,10000, by = 0.1), Npool = seq(0.0001,10000, by = 0.1)*90)
NP_100 = data.frame(Stoic = "NP_100", Ppool = seq(0.0001,10000, by = 0.1), Npool = seq(0.0001,10000, by = 0.1)*100)
NP_200 = data.frame(Stoic = "NP_200", Ppool = seq(0.0001,10000, by = 0.1), Npool = seq(0.0001,10000, by = 0.1)*200)
NP_300 = data.frame(Stoic = "NP_300", Ppool = seq(0.0001,10000, by = 0.1), Npool = seq(0.0001,10000, by = 0.1)*300)
NP_400 = data.frame(Stoic = "NP_400", Ppool = seq(0.0001,10000, by = 0.1), Npool = seq(0.0001,10000, by = 0.1)*400)
NP_500 = data.frame(Stoic = "NP_500", Ppool = seq(0.0001,10000, by = 0.1), Npool = seq(0.0001,10000, by = 0.1)*500)
NP_1000 = data.frame(Stoic = "NP_1000", Ppool = seq(0.0001,10000, by = 0.1), Npool = seq(0.0001,10000, by = 0.1)*1000)
NP_2000 = data.frame(Stoic = "NP_2000", Ppool = seq(0.0001,10000, by = 0.1), Npool = seq(0.0001,10000, by = 0.1)*2000)
NP_3000 = data.frame(Stoic = "NP_3000", Ppool = seq(0.0001,10000, by = 0.1), Npool = seq(0.0001,10000, by = 0.1)*3000)
NP_4000 = data.frame(Stoic = "NP_4000", Ppool = seq(0.0001,10000, by = 0.1), Npool = seq(0.0001,10000, by = 0.1)*4000)
NP_5000 = data.frame(Stoic = "NP_5000", Ppool = seq(0.0001,10000, by = 0.1), Npool = seq(0.0001,10000, by = 0.1)*5000)

NP = data.frame(rbind(NP_0.1, NP_0.2,NP_0.3,NP_0.4, NP_0.5, NP_0.6, NP_0.7, NP_0.8, NP_0.9, NP_1, NP_2, NP_3, NP_4, NP_5, NP_6, NP_7, NP_8, NP_9, NP_10,
                      NP_20, NP_30, NP_40, NP_50, NP_60, NP_70, NP_80, NP_90, NP_100, NP_200, NP_300, NP_400, NP_500, NP_1000, NP_2000, NP_3000, NP_4000, 
                      NP_5000), stringsAsFactors = T)
NP$Stoic = factor(NP$Stoic)
 NP = NP[-which(NP$Npool > 3800 | NP$Npool < 0.9 | NP$Ppool > 120 | NP$Ppool < 0.1),]
# lines = c("solid", rep("dotted", 3), "dashed", rep("dotted", 4), "solid", "solid", "solid", rep("dotted", 6), "dashed", "dashed",rep("dotted", 8), rep("dotted", 4))
lines = c(rep(c('solid','dotted','dotted','dotted','dashed', 'dotted','dotted','dotted','dotted'),3),
          'solid','dotted','dotted','dotted','dashed','solid','dotted','dotted','dotted','dashed')
p1 = ggplot() +
  geom_line(data = NP, aes(x = Ppool, y = Npool, group = Stoic, linetype = Stoic), alpha = 0.7, size = 0.5, colour = "black") +
  scale_linetype_manual(values = lines) + theme_bw()+
  theme(panel.grid = element_blank(), legend.position = "none", axis.text = element_text(size = 20))#;p1

p1 = p1 + geom_point(data = lakes2012_full, aes(x = P_mol, y = N_mol), fill = 'red',  shape = 21, colour = "black", size = 2) + 
  xlab(expression("Phosphorus ("~mu*"mol)")) +
  ylab(expression("Nitrogen ("~mu*"mol)"))

p2 = p1 +   
   coord_cartesian(ylim = c(0.8,4000), 
                   xlim = c(0.08, max(lakes2012_full$P_mol))) +
  scale_y_log10() + scale_x_log10();p2

diss_NP.lm = lm(log10(N_mol)~log10(P_mol), data = lakes2012_full)

p3 = p2 + geom_text(aes(x = 0.09, y = 9.5, label = "100", size = 30), colour = "grey50", fontface = "bold", hjust = 1) +
  geom_text(aes(x = 0.09, y = 4.75, label = "50", size = 30), colour = "grey50", fontface = "bold", hjust = 1) +
  geom_text(aes(x = 0.09, y = 0.95, label = "10", size = 30), colour = "grey50", fontface = "bold", hjust = 1) +
  geom_text(aes(x = 0.09, y = 47.5, label = "500", size = 30), colour = "grey50", fontface = "bold", hjust = 1) +
  geom_text(aes(x = 0.09, y = 95, label = "1000", size = 30), colour = "grey50", fontface = "bold", hjust = 1) +
  geom_text(aes(x = 0.09, y = 450, label = "5000", size = 30), colour = "grey50", fontface = "bold",hjust = 1) +
  geom_text(aes(x = 120, y = 11, label = "0.1", size = 30), colour = "grey50", fontface = "bold",hjust = 0) +
  geom_text(aes(x = 120, y = 60, label = "0.5", size = 30), colour = "grey50", fontface = "bold",hjust = 0) +
  geom_text(aes(x = 120, y = 120, size = 30, label = "1"), colour = "grey50", fontface = "bold",hjust = 0) +
  geom_text(aes(x = 120, y = 600, label = "5", size = 30), colour = "grey50", fontface = "bold",hjust = 0) +
  geom_text(aes(x = 120, y = 1200, label = "10", size = 30), colour = "grey50", fontface = "bold",hjust = 0) +
  # geom_text(aes(x = 0.88, y = 10, label = "10", size = 30), colour = "grey50", fontface = "bold") +
  # geom_text(aes(x = 0.88, y = 50, label = "50", size = 30), colour = "grey50", fontface = "bold") +
  # geom_text(aes(x = 0.88, y = 100, label = "100", size = 30), colour = "grey50", fontface = "bold") +
  annotate("text",x = 0.08, y = 3000, label = paste("N:P scaling\nexponent =",round(diss_NP.lm$coefficients[2],2),
                                                   " ± ",round(summary(diss_NP.lm)$coefficients[4],2),sep=""), size = 3, hjust = 0);p3 


########## 
lakes2012_LU = lakes2012_full %>% select(Agriculture:WetlandWater) 

lakes2012_LU_sum = rowSums(lakes2012_LU)

landcover_mds = vegan::metaMDS(lakes2012_LU, distance = "bray", trymax =1000, autotransform = T)
plot()








######## END RMARKDOWN CODE ######


s
#saveRDS(lakes_landuse, file = "./ecosystem_scripts/Lakes2012/lake_landuse_list.rds")
#attempt to speed it with lapply functions.
ten_lakes = list()
debugonce(get_nlcd_percents)
tic();ten_lakes = map_dfr(basin_shapes_ex, get_nlcd_percents);toc()

##need to work out why getting error "Error in st_sfc(x, crs = attr(x, "proj4string")) : 
#is.numeric(crs) || is.character(crs) || inherits(crs, "crs") is not TRUE"
as(basin_shapes_ex[1,], "Spatial")
as(list(basin_shapes_ex[1,]), "Spatial")
#the projection is NAD83 EPSG code (for crs = ) == 4269
#WSG 84 = 4326

lakes_loc <- lakes2012_site %>%
  st_as_sf(coords = c("LON_DD83","LAT_DD83"), crs = 4269)

#look at the first few sites
lakes_loc[1:5,] %>%
  st_transform(crs = 4269) %>%
  leaflet() %>%
  addTiles() %>%
  addMarkers()

#plot them all
plot(lakes_loc$geometry)

# Make buffer around first point

pt1_buffer <- lakes_loc %>%
  dplyr::slice(1) %>%
  st_buffer(0.01) #units in decimal degree. Should ~ 1000m

# plot buffer
leaflet(pt1_buffer) %>%
  addProviderTiles(providers$Esri.WorldImagery) %>%
  addPolygons()

# use buffer to get nlcd function

nlcd_nars  <- get_nlcd(as(pt1_buffer, "Spatial"), 
                       label = "first lake site",
                       dataset = "landcover",
                       year = "2011",
                       force.redo = TRUE)### had an error message pop up saying corrupt or incomplete .zip file

nlcd_nars %>% plot()##seems to be a pretty painting

leaflet(pt1_buffer) %>%
  addProviderTiles(providers$Esri.WorldImagery) %>%
  addPolygons() %>%
  addRasterImage(nlcd_nars)

#Mask raster using buffer object after transforming to same CRS

pt1_buffer_prj <- st_transform(pt1_buffer, crs = proj4string(nlcd_nars))
nlcd_nars_mask <- raster::mask(nlcd_nars, as(pt1_buffer_prj, "Spatial"))

plot(nlcd_nars_mask)

#tabulate number of cells in each type and merge with legend to see land cover types

lc_legend_df <- readr::read_csv("nlcd_legend_2011.csv")

raster::freq(nlcd_nars_mask) %>%
  as.data.frame() %>%
  dplyr::filter(!is.na(value)) %>%
  dplyr::left_join(lc_legend_df, by = c('value' = 'Class')) %>%
  dplyr::mutate(percent_cover = count/sum(count)) %>%
  dplyr::select(class_name, percent_cover)

##this seems to work well based on a circular buffer
##need to think about how to apply this to lakes and streams??

## Define function to do all of the above --get percent of land cover types
# within polygon. Returns polygon with columns class_name and percent cover

# this works for a single one now iterate across a subset of rows in lakes_loc
#automate the buffer getting for a subset of lakes

lakes_loc_ex <- lakes2012_site[1:10,]  %>%
  st_as_sf(coords = c("LON_DD83","LAT_DD83"), crs = 4269)

lakes_buffer <- lakes_loc_ex %>%
  st_buffer(0.01)

ten_lakes = list()
tic();for(i in seq_along(lakes_buffer[[1]]))ten_lakes[[i]] <- get_nlcd_percents(lakes_buffer[i,]);toc()
# for loop works, but likely slow. Prefer to use apply-group. Working on that
# got help from here: https://stackoverflow.com/questions/48505551/use-apply-with-a-simple-features-sf-function


debugonce(get_nlcd_percents)
ten_lakes = lapply(st_geometry(lakes_buffer), FUN = get_nlcd_percents)

##need to work out why getting error "Error in st_sfc(x, crs = attr(x, "proj4string")) : 
#is.numeric(crs) || is.character(crs) || inherits(crs, "crs") is not TRUE"

crs(lakes_buffer)#;crs(lakes_buffer)<-CRS("+proj=longlat +epsg=4269")
as(lakes_buffer, 'Spatial') 
st_geometry(lakes_buffer)
as(st_geometry(lakes_buffer),'Spatial')
crs(as(st_geometry(lakes_buffer),'Spatial')) <-CRS("+proj=longlat +epsg=4269")
####
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
  pull(filename)  %>% '[[' (2)# because there are 2 files.. 
  paste0("raw-data/", .) %>%
  read_csv()
  
nars_data %>%
  filter(Survey == "Rivers and Streams 2008-2009",
         Indicator == "Benthic Macroinvertebrates") %>%
  pull(filename) %>% unique()
  
rivers2008_benth.inv <- nars_data %>%
  filter(Survey == "Rivers and Streams 2008-2009",
         Indicator == "Benthic Macroinvertebrates") %>%
  pull(filename) %>% '[[' (5) %>%
  paste0("raw-data/",.) %>%
  read.csv
  
########### End other system code ##################