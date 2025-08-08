# Assessing Mature and Old Growth Forest using FIA data
Mature and old growth forest support vast quantities of life, and are critical to the function of many ecosystems. While remotely sensed landcover datasets can indicate the spatial distribution of forests, they offer little insight to the age of forests. Fortunately, the United States Forest Service maintains an extensive forest inventory program across public and private lands called the [Forest Inventory Analysis](https://research.fs.usda.gov/programs/fia) (FIA). This function uses R (specifically the [RFIA package](https://doserlab.com/files/rfia/)) to access FIA data from the Forest Service data repository, then assesses the data to determine if it meets the definitions for [mature](https://www.sciencedirect.com/science/article/abs/pii/S0378112723005959) or [old growth](https://www.sciencedirect.com/science/article/abs/pii/S0378112723006710) forest. Definitions vary by region and forest-type, offering substantial flexibility even across ecosystems. Finally, this function offers the option to interpolate the mature/old growth (MOG) status of unsampled forests based on neighboring sampled forests. The following tutorial explains how to use the function.

## Identifying a study area
The `mapMOG()` function requires the user to identify a study area using a shapefile in `sf` format. The study area, called a `locale`, can be as large as a state, or as small as a forest. Please note, however, that the function cannot currently accommodate multi-state processing. If your study `locale` spans multiple states, the function must be run separately for each state contained in the `locale`. Additionally, small forests are not likely to contain enough FIA plots for adequate interpolation. Assessing data at the landscape-level (e.g., county, etc.) is likely to produce greater insights, even if you are only interested in a small portion of the landscape. In this tutorial we will use part of Olympic National Forest as a case study. A shapefile of Olympic National Forest is included in the function folder for your convenience.

## Downloading and Accessing the function and data
[Download the entire mapMOG folder from Dryad](https://drive.google.com/drive/folders/1EM1pKgUz62EmD2MrMg4RpMkV43ndUvaU?usp=sharing) (file sizes exceed Github's limits - email me if the link is broken), and unzip the folder if necessary. We recommend moving the folder from the Downloads folder of your computer to a more permanent location. **Do not remove files from the `mapMOG` folder, as these files are necessary for proper functionality.**

File size limitations barred us from including all necessary files in the ‘mapMOG’ folder available for download. Prior to running the function, download the following files and add them to the mapMOG folder:

*If you plan on interpolating MOG across unsampled landscapes*
1.	Download the most recent canopy cover layer from the National Land cover Database: https://www.mrlc.gov/data
2.	Unzip the downloaded folder into the ‘utility_canopy’ subfolder of ‘mapMOG’ such that the file structure is ‘[YOUR FILE PATH]/mapMOG/utility_canopy/[CANOPY RASTER]’
3.	No need to rename the raster
   
*If you plan on assessing MOG in the Pacific Northwest Region or Pacific Southwest Region*
1.	Download the PNV_Map_Layer (PNV_Map_zip_20220127.zip) from box: https://app.box.com/s/hfo4y7f9tg4vetm83zl4ahzhv0grqyqj/folder/166077218938
2.	Unzip the downloaded folder into your ‘mapMOG’ folder such that the file structure is ‘[YOUR FILE PATH]/mapMOG/[PNV RASTER.tif]’
3.	Rename the raster “utility_R5_PAZ” to ensure functionality


Once you have added the canopy and PNV files to your `mapMOG` folder (if necessary), open R Studio and load the function by sourcing the script from the `mapMOG` folder. Note that the entire filepath must be specified for your computer to find the file.
```
source("[YOUR FILE PATH]/mapMOG/FUNCTION_mapMOG.R") #loads function into R environment
```

To access the Olympic National Forest shapefile, use the `getOlympic()` function which is now saved in your environment. If you are using your own shapefile, ensure that the file is in `sf` format. The function will reproject the shapefile as needed. I use `ggplot2` to visualize shapefiles in this tutorial, but base R or other plotting packages can also be used.

```

my.folderpath <- "[YOUR FILE PATH]/mapMOG" #save file path to mapMog FOLDER (not function)

olympicNF <- getOlympic(source.path = my.folderpath) #reads in Olympic National Forest Shapefile

library(ggplot2) #load ggplot2

ggplot()+
   geom_sf(data = olympicNF, #plot the shapefile of Olympic National Park
           fill = "white", #make the inside of the shape white
           color = "black")+ #make the outline of the shape black
   coord_sf(crs = 4326) #reproject so it is oriented correctly (WGS84)
```
![olympicPoly](https://github.com/user-attachments/assets/723ef0cc-58ef-4d98-9540-83131b61eaf8)


## Using the mapMOG function
Two arguments are needed to run the function: the shapefile of the `locale`, and the file path of the `mapMOG` folder. The `recent` and `interpolate` arguments both default to `TRUE` but can be changed if desired. Only the most recent FIA plot data are used when `recent = TRUE`, and a raster of interpolated mature and old growth forests is produced when `interpolate = TRUE`. Finally, `api` defaults to the `rFIA` package. At this point, the `rFIA` package is the only supported API to access FIA data. Future versions of this function may incorporate alternative API's.

```
mapMOG(locale = olympicNF, #specify your study area
          source.path = my.folderpath) #specify the location of the folder that holds the function and auxiliary data
```

When you run the `mapMOG()` function, updates will print in your console ensuring that data has been successfully obtained. Error messages may also print if data are not available, etc. The `mapMOG()` function requires internet access to obtain data, but no error message will warn the user that the computer does not have access to internet. If an especially cryptic error message appears, check your internet connection.

## Visualizing and interpreting the output
After the function is finished running, it will produce two outputs: a shapefile of FIA points classified by their MOG score, and a raster of interpolated MOG scores for all forested land in the `locale`. Note that these objects are automatically named `MOG.points` and `MOG.raster` Let's take a look at the FIA points first:

```
ggplot()+
 geom_sf(data = olympicNF, #plot the shapefile of Olympic National Park
         fill = "white", #make the inside of the shape white
         color = "black")+ #make the outline of the shape black
 geom_sf(data = MOG.points, #plot FIA points
         aes(color = p.mog))+ #color FIA points by their MOG score
 coord_sf(crs = 4326) #reproject so it is oriented correctly (WGS84)
```
![olympicPoints](https://github.com/user-attachments/assets/223c8479-d5b2-489d-a03e-d001c4d416b3)

Each of the points on the map is a single FIA plot. The color of each point indicates the MOG score of that point. Scores range between zero and one, with zero indicating non-mature or old growth and one indicating old growth. A value of 0.5 indicates the forest is mature. In this example, all FIA plots in this section of Olympic National Forest are classified as being mature or old growth forest since all points have a score of 0.5 or greater. You may notice that some FIA points fall outside of the national forest. The function retrieves and processes data within 5 kilometers of the `locale` to minimize inaccurate interpolation at the locale edge due to sparse data.

You may notice that there are sections of the forest which do not have FIA data. We can use the interpolated MOG raster to offer insight about the MOG status of those areas based on the status of the three nearest FIA points. The fastest way to visualize this raster is to use base R: `plot(MOG.raster)`. However, for the sake of consistency, we will use `ggplot2()` in this example.

```
#library(terra) #if you just ran the mapMOG function, the terra package is already loaded. No need to load again

MOG.raster.transformed <- terra::project(x = MOG.raster, #reproject raster to ensure sensible plotting
                                  y = "epsg:4326") #WGS84
plot.df <- as.data.frame(x = MOG.raster.transformed, #convert raster into a dataframe
                         xy = TRUE) #save centroid coordinates of each cell
ggplot()+
  geom_tile(data = plot.df, #plot raster
              aes(x = x,
                  y = y,
                  fill = MOG))+ #color cells by MOG score
  geom_sf(data = olympicNF, #plot Olympic NF over raster to see administrative boundaries
          fill = "transparent", 
          color = "white")+ 
  coord_sf(crs = 4326) #reproject so it is oriented correctly (WGS84)
```
![olympicRaster](https://github.com/user-attachments/assets/abbac1cd-a446-4acc-b19e-f36eab974fd5)

While interpolation can be helpful, it is important to remember that this is simply an estimation. Unsampled forest areas are estimated using inverse-distance weighting at a 1600x1600 meter (~1 mile) resolution. FIA plot data have slightly obscured geographic coordinates to protect landowner safety, and the interpolation process only occurs in raster cells that meet a canopy cover threshold specific to the `locale`.

Dynamic data exploration can be achieved using the `mapview` package: `mapview(MOG.raster.transformed)`

If you encounter issues with the `mapMOG()` function, please email Dan Herrera at *herrerawildlife(at)gmail.com*.
