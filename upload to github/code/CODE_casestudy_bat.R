rm(list=ls())
library(sf) #vector processing
library(terra) #raster processing
library(raster) #raster processing
library(unmarked) #easily runs occupancy models
library(AICcmodavg) #easily constructs AIC table

#note that this script uses all variables at a 30x30-m scale. 
#to run these analyses at the native 1,590x1,590-m scale, use the coarser MOG dataset, and re-size the canopy and forest rasters accordingly.

species.of.interest <- "LASE" #Seminole bat


#load non-MOG forest data
ecoregions <- st_read("./data/utility_ecoregions/USA_ecoregions.shp")
ecoregions <- st_transform(x = ecoregions,
                           crs = 4326) #WGS84
USA <- st_read("./data/utility_USA/tl_2020_us_state.shp")
USA <- st_transform(x = USA,
                    crs = 4326) #WGS84

region <- ecoregions[ecoregions$NA_L2NAME %in% "MISSISSIPPI ALLUVIAL AND SOUTHEAST USA COASTAL PLAINS", ]

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

#load observation data
obs <- read.csv("./data/NAbat - bat data.csv")
obs.geom <- st_as_sf(st_as_sfc(structure(obs$grts_geometry, class = "WKB"), EWKB = TRUE)) #convert from wkb to something legible
obs.geom <- st_transform(obs.geom, #transform to match crs of spatial grid
                         crs = 4326) #WGS84
st_geometry(obs.geom) <- "geometry"
bat.obs <- cbind(obs.geom,obs) #combine bat observations with geometry
bat.obs <- st_filter(x = bat.obs, #subset to bat observations within the target region
                     y = region)

bat.obs <- bat.obs[substr(bat.obs$start_time,0,4) %in% "2022", ] #subset to 2023 data

sites <- as.vector(unique(bat.obs$geometry)) #create a vector of sites to loop through

data.df <- data.frame()
for(i in 1:length(sites)){ #begin looping through sites
  
  local.site <- bat.obs[bat.obs$geometry %in% sites[i], ]
  local.site <- local.site[1,"geometry"]
  data.df <- rbind(data.df, local.site)
  
  message(paste(i,"of",length(sites),sep=" "))
} #end looping through sites

sites <- data.df #rename as these are now your sites (individual sampling points will be related to their respecitive parent block)

df.detection <- data.frame()
df.effort <- data.frame()
df.forest <- data.frame()
for(i in 1:nrow(sites)){ #begin looping through sites (again)
  
  local.site <- sites[i,]
  
    
    local.obs <- st_filter(x = bat.obs,
                           y = local.site)
    
    sample.nights <- as.vector(unique(local.obs$night)) #create a vector of sample nights
    local.detection <- matrix(nrow = 1, ncol = 365, data = NA) #create an empty matrix to fill with values
    df.effort[i,1] <- length(sample.nights) #save the number of sampling nights
    
    for(n in 1:length(sample.nights)){ #begin looping through sample nights
      
      nightly.obs <- local.obs[local.obs$night %in% sample.nights[n], ] #subset to observations on night n
      df.effort[i,n] <- length(unique(nightly.obs$location_name)) #save the number of active sampling sites for the evening
      
      nightly.obs <- nightly.obs[nightly.obs$species_code %in% species.of.interest, ] #subset to the species of interest
      sum.nightly.obs <- sum(as.numeric(nightly.obs$count_auto_id))
      if(sum.nightly.obs > 0){df.detection[i,n] <- 1}
      if(sum.nightly.obs == 0){df.detection[i,n] <- 0}
      
    } #end looping through sample nights
    
    temp.site <- st_transform(x = local.site, #transform to crs of forest raster
                              crs = st_crs(forest))
    forest.area <- terra::extract(x = forest,
                                  y = vect(temp.site))
    df.forest[i,"forest.area"] <- (sum(forest.area[,2], na.rm=TRUE)*(30*30))/ as.numeric(st_area(temp.site)) #proportion of park that is forested
    if(df.forest[i,"forest.area"] > 1){df.forest[i,"forest.area"] <- 1} #if rounding/imperfect alignment gets in the way, etc. truncate back down to 1
    
    temp.site <- st_transform(x = local.site, #transform to crs of forest raster
                              crs = st_crs(canopy))
    canopy.cover <- terra::extract(x = canopy,
                                   y = vect(temp.site))
    df.forest[i,"canopy.cover"] <- mean(as.numeric(canopy.cover[,2]), na.rm=TRUE)
    
    temp.site <- st_transform(x = local.site, #transform point to crs of MOG raster
                              crs = st_crs(mog))
    mog.site <- terra::extract(x = mog,
                               y = temp.site)
    df.forest[i,"MOG"] <- mean(mog.site[,2], na.rm = TRUE)
    
  
  message(paste(i,"of",nrow(sites)))
} #end looping through sites

keep.index <- which(!is.na(df.detection$V1)) #identify sites that don't have missing data
df.detection <- as.data.frame(df.detection[keep.index, ])
df.effort <- as.data.frame(df.effort[keep.index, ])
df.forest <- as.data.frame(df.forest[keep.index, ])

#put data into unmarked frame for single-species single-season occupancy analysis
umf <- unmarkedFrameOccu(y = df.detection,
                         siteCovs = df.forest,
                         obsCovs = list(effort = df.effort)) #this needs to be a list

#run competing occupancy models
model.forest <- occu(formula = ~ scale(effort) ~scale(forest.area),
                     data = umf)
model.canopy <- occu(formula = ~scale(effort) ~scale(canopy.cover),
                     data = umf)
model.mog <- occu(formula = ~scale(effort) ~scale(MOG),
                  data = umf) 

cand.models <- list()
cand.models[[1]] <- model.canopy
cand.models[[2]] <- model.forest
cand.models[[3]] <- model.mog
aictab(cand.models,sort=TRUE)