# Assessing Mature and Old Growth Forest using FIA data
Mature and old growth forest support vast quantities of life, and are critical to the function of many ecosystems. While remotely sensed landcover datasets can indicate the spatial distribution of forests, they offer little insight to the age of the forest. Fortunately, the United States Forest Service maintains an extensive forest inventory program across public and private lands called the [Forest Inventory Analysis](https://research.fs.usda.gov/programs/fia) (FIA). This function uses R to access FIA data from the [Forest Service data repository](https://research.fs.usda.gov/products/dataandtools/tools/fia-datamart), then assess the data to determine if it meets the definitions for [mature](https://www.sciencedirect.com/science/article/abs/pii/S0378112723005959) or [old growth](https://www.sciencedirect.com/science/article/abs/pii/S0378112723006710) forest. Definitions vary by region and forest-type, offering substantial flexibility even across ecosystems. Finally, this function offers the option to interpolate the mature/old growth status of unsampled forests based on neighboring sampled forests. The following tutorial explains how to use the function.

## Identifying a study area
The `assessMOG()` function requires the user to identify a study area using a shapefile in `sf` format. The study area, called a `locale`, can be as large as a state, or as small as a forest. Please note, however, that the function cannot currently accomodate multi-state processing. If your study `locale` spans multiple states, the function must be run seperately for each state contained in the `locale`. Additionally, small forests are not likely to contain enough FIA plots for adequate interpolation. Results are expected to be most insightful when the `locale` is greater than 160 square kilometers (approximately 10x10 miles). In this tutorial we will use Pisgah National Forest as a case study. A shapefile of Pisgah National Forest is included in the function folder for your convenience.

## Downloading and Accessing the function and data
Download the entire `assessMOG` folder from this Github repository, and unzip the folder if necessary. We recomend moving the folder from the Downloads folder of your computer to a more perminant location. **Do not add or remove files from the `assessMOG` folder, as these files are necessary for proper functionality.**

Once you have saved the folder, open R Studio and load the function by sourcing the script from the `assessMOG` folder. Note that the entire filepath must be specified for your computer to find the file.
```
source("[YOUR FILE PATH]/assessMOG/FUNCTION_assessMOG.R") #loads function into R environment
```

To access the Pisgah National Forest shapefile, use the `getPisgah()` function. If you are using your own shapefile, ensure that the file is in `sf` format. The function will reproject the shapefile as needed.

```
pisgah <- getPisgah() #reads in Pisgah National Forest Shapefile
plot(pisgah$geometry) #plot a map of the forest shapefile
```
DAN!!!! INSERT A MAP HERE!

## Using the assessMOG function
Two arguments are needed to run the function: the shapefile of the `locale`, and the file path of the `assessMOG` folder. The `recent` and `interpolate` arguments both default to `TRUE` but can be changed if desired. Only the most recent FIA plot data are used when `recent = TRUE`, and a raster of interpolated mature and old growth forests is produced when `interpolate = TRUE`. Finally, `api` defaults to the `rFIA` package. At this point, the `rFIA` package is the only supported API to access FIA data. Future versions of this function may incroporate alternative API's.

```
assessMOG(locale = pisgah, #sf shapefile of study area
          source.path = "[YOUR FILE PATH]/assessMOG") #file path to the assessMOG folder
```
