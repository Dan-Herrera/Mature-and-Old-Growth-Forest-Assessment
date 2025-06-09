rm(list=ls())
library(sf) #vector processing
library(terra) #raster processing
library(raster) #raster processing
library(maxlike) #easily runs maxlike models (similar to maxent models, but probabilistic)
library(AICcmodavg) #easily constructs AIC table

#note that this script uses all variables at a 30x30-m scale. 
#to run these analyses at the native 1,590x1,590-m scale, use the coarser MOG dataset, and re-size the canopy and forest rasters accordingly.

USFS <- st_read("./data/Forest Service Boundaries/S_USA.AdministrativeForest.shp") #load forest service boundaries
USFS <- st_transform(x = USFS,
                     crs = 4326) #WGS84
USFS <- st_make_valid(USFS) #fix potential geometry issues within USFS shapefiles
USFS <- USFS[USFS$FORESTNAME %in% c("Gifford Pinchot National Forest", #select national forests within our study area
                                    "Columbia River Gorge National Scenic Area",
                                    "Mt. Hood National Forest",
                                    "Willamette National Forest"),]
PNW <- st_buffer(x = USFS,
                  dist = 10000)
PNW <- st_union(PNW)
PNW <- st_as_sf(PNW)
PNW <- vect(PNW)

#humans <- terra::rast("./human pop dens/usa_pd_2020_1km.tif") #load human population density data (UN OCHA: https://data.humdata.org/dataset/worldpop-population-density-for-united-states-of-america)
canopy <- terra::rast("./data/nlcd_tcc_conus_2021_v2021-4.tif") #load canopy cover dataset
NLCD <- terra::rast("./data/Annual_NLCD_LndCov_2023_CU_C1V0.tif") #load national landcover database (NLCD)

PNW <- project(x = PNW,
               y = crs(canopy))
#crop land covers to study area
canopy <- terra::crop(x = canopy,
                      y = PNW,
                      mask = TRUE, touches = TRUE)
canopy[canopy>100] <- 0 #replace values greater than 100 with zeros (ocean is listed as 250)

PNW <- project(x = PNW,
               y = crs(NLCD))
NLCD <- terra::crop(x = NLCD,
                    y = PNW,
                    mask = TRUE, touches = TRUE)

rules <- c(11,0, #create reclassification rules for landcovers
           12,0,
           21,0,
           22,0,
           23,0,
           24,0,
           31,0,
           41,1, #deciduous forest
           42,1, #evergreen forest
           43,1, #mixed forest
           51,0,
           52,0,
           71,0,
           72,0,
           73,0,
           74,0,
           81,0,
           82,0,
           90,1, #woody wetlands (aka swamp)
           95,0)
rules <- matrix(rules, ncol=2, byrow = TRUE) #convert rules to matrix (required)
forest <- classify(x = NLCD, #reclassify so all forest is 1 and non-forest is 0
                   rcl = rules)

forest <- terra::project(x = forest,
                  y = canopy)


#load MOG data
mog <- terra::rast("./data/USA_MOG_30x30.tif")

#crop mog to study area
mog <- terra::project(x = mog,
                      y = canopy)
PNW <- terra::project(x = PNW,
               y = crs(mog))
mog[is.na(mog)] <- 0 #convert NA's (non-mog) to zero (also non-mog)
mog <- terra::crop(x = mog,
                   y = PNW,
                   mask = TRUE, touches = TRUE)

beep(sound = "fanfare")


#load presence data for cascade torrent salamander
data.1 <- read.csv("./data/salamander data/RHCA_table_FY19report.csv")
data.2 <- st_read("./data/salamander data/rhca_merge_elev_asp_slp.shp")

#transform the two datasets to the same projections
data.2 <- st_transform(x = data.2,
                       crs = crs(NLCD))
data.1 <- st_as_sf(x = data.1,
                   coords = c("Lon..WGS84.","Lat..WGS84."),
                   crs = 4326) #WGS84 (original data source)
data.1 <- data.1[data.1$RHCA > 0, ] #thin to observations of focal species
data.1 <- st_transform(x = data.1,
                       crs = st_crs(data.2))
data.1 <- data.1[,"geometry"]
data.2 <- data.2[,"geometry"]
salamanders <- rbind(data.1,data.2) #combine salamander datasets
salamanders <- st_transform(x = salamanders,
                            crs = crs(canopy))
salamanders <- st_intersection(x = salamanders, #clip salamander observations to the dataset
                               y= st_as_sf(PNW))

#prepare data for model
forest <- raster(forest) #convert from terra back to raster
canopy <- raster(canopy)
mog <- raster(mog)


#place each variable in its own raster stack (required by model)
forest.stack <- stack(scale(forest))
names(forest.stack) <- "forest"
canopy.stack <- stack(scale(canopy))
names(canopy.stack) <- "canopy"
mog.stack <- stack(scale(mog))
names(mog.stack) <- "mog"

#prepare coordinate data for salamander detections (presence only)
xy.data <- as.data.frame(st_coordinates(salamanders))


#run MaxLike models
model.forest <- maxlike(formula = ~forest,
                        rasters = stack(forest.stack),
                        points = xy.data) 

model.canopy <- maxlike(formula = ~canopy,
                        rasters = stack(canopy.stack),
                        points = xy.data) 

model.mog <- maxlike(formula = ~mog,
                     rasters = stack(mog.stack),
                     points = xy.data) 

cand.models <- list()
cand.models[[1]] <- model.canopy
cand.models[[2]] <- model.forest
cand.models[[3]] <- model.mog
aictab(cand.models,sort=TRUE)
