rm(list=ls())
library(sf) #vector processing
library(terra) #raster processing
library(raster) #raster processing 
library(AICcmodavg) #easily constructs AIC table

#note that this script uses all variables at a 30x30-m scale. 
#to run these analyses at the native 1,590x1,590-m scale, use the coarser MOG dataset, and re-size the canopy and forest rasters accordingly.

#load non-MOG forest data
USA <- st_read("./data/utility_USA/tl_2020_us_state.shp")
USA <- st_transform(x = USA,
                    crs = 4326) #WGS84

canopy <- terra::rast("./data/nlcd_tcc_conus_2021_v2021-4.tif") #load canopy cover dataset
NLCD <- terra::rast("./data/Annual_NLCD_LndCov_2023_CU_C1V0.tif") #load national landcover database (NLCD)
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

mog <- terra::rast("./data/USA_MOG_30x30.tif")

NPS <- st_read("./data/NPS.parcels/National_Park_Service_Land_Resources_Division_Tract_and_Boundary_Service.shp")
region <- USA[USA$STUSPS %in% c("VA","WV","MD","DE","NJ","PA","NY"), ]
region <- st_transform(x = region,
                       crs = st_crs(NPS))
local.parks <- st_filter(x = NPS,
                         y = region)

#load observation data
bird.data <- read.csv("./data/midatlanticNPS.csv") #load NPS bird checklist data (https://irma.nps.gov/NPSpecies/Search/SpeciesList)
bird.data <- bird.data[!is.na(bird.data$nBirdSpp), ] #keep only parks with bird data
bird.data <- bird.data[!(bird.data$UNIT_NAME %in% "Appalachian National Scenic Trail"), ] #remove Appalachian trail since it goes well beyond our study region

for(i in 1:nrow(bird.data)){ #begin looping through bird data
  
  focal.park <- local.parks[local.parks$UNIT_CODE %in% bird.data[i,"UNIT_CODE"],]
  
  bird.data[i,"area.m2"] <- as.numeric(st_area(focal.park)) #save area of park
  
  temp.park <- st_transform(x = focal.park, #transform to crs of forest raster
                            crs = st_crs(forest))
  forest.area <- terra::extract(x = forest,
                                y = vect(temp.park))
  bird.data[i,"forest.area"] <- (sum(forest.area[,2], na.rm=TRUE)*(30*30))/ bird.data[i,"area.m2"] #proportion of park that is forested
  if(bird.data[i,"forest.area"] > 1){bird.data[i,"forest.area"] <- 1} #if rounding/imperfect alignment gets in the way, etc. truncate back down to 1
  
  temp.park <- st_transform(x = focal.park, #transform to crs of forest raster
                            crs = st_crs(canopy))
  canopy.cover <- terra::extract(x = canopy,
                                 y = vect(temp.park))
  bird.data[i,"canopy.cover"] <- mean(as.numeric(canopy.cover[,2]), na.rm=TRUE)
  
  
  temp.park <- st_transform(x = focal.park, #transform point to crs of MOG raster
                            crs = st_crs(mog))
  mog.park <- terra::extract(x = mog,
                             y = temp.park)
  bird.data[i,"MOG"] <- mean(mog.park[,2], na.rm = TRUE)
  
  message(i)
  
} #end looping through bird data


#run competing models
canopy.model <- glm(formula = bird.data$nBirdSpp ~ scale(bird.data$canopy.cover) + scale(bird.data$area.m2),
                    family = "poisson")
forest.model <- glm(formula = bird.data$nBirdSpp ~ scale(bird.data$forest.area)+ scale(bird.data$area.m2),
                    family = "poisson")
mog.model <- glm(formula = bird.data$nBirdSpp ~ scale(bird.data$MOG)+ scale(bird.data$area.m2),
                 family = "poisson")


cand.models <- list()
cand.models[[1]] <- canopy.model
cand.models[[2]] <- forest.model
cand.models[[3]] <- mog.model
aictab(cand.models,sort=TRUE)
