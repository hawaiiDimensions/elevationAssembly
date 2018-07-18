## IMPORT CLIMATIC RASTERS 
# Author: Jun Ying Lim
# 

## PACKAGES ============
library(geosphere) # calculating geographic distances
library(raster)
library(sp); library(rgdal)

## DEFINE DIRECTORIES ============
main.dir <- "~/Dropbox/Projects/2017/hawaiiCommunityAssembly/"
analysis.dir <- file.path(main.dir, "elevationAssemblyHawaii")
data.dir <- file.path(main.dir, "data")
fig.dir <- file.path(analysis.dir, "figures")
map.dir <- "/Users/junyinglim/Dropbox/Projects/2014/Hawaii/Maps"
source(file.path(analysis.dir, "ecoDataTools.R"))

## DEFINE SPATIAL EXTENT ============
siteData <- read.csv(file.path(data.dir, "clim.final.csv"), stringsAsFactors = FALSE)
siteData <- siteData[-grep(siteData$site.id, pattern = "BRG"),] # Exclude Rosie's samples
siteExtent <- extent(min(siteData$longitude)-0.1,
                     max(siteData$longitude)+0.1,
                     min(siteData$latitude)-0.1,
                     max(siteData$latitude)+0.1)

## IMPORT RASTERS ============
cloudFreqRaster <- raster(file.path(map.dir, "CloudFreq_month_ascii", "cl_frq_ann.txt"))
cloudFreqRaster_cropped <- crop(cloudFreqRaster, siteExtent)

annTempRaster <- raster(file.path(map.dir, "Tair_month_ascii", "tair_ann.txt"))
annTempRaster_cropped <- crop(annTempRaster, siteExtent)

annPrecipRaster <- raster(file.path(map.dir, "StateASCIIGrids_mm", "rfgrid_mm_state_ann.txt"))
annPrecipRaster_cropped <- crop(annPrecipRaster, siteExtent)

plot(annPrecipRaster)

## EXPORT CROPPED RASTERS ============
saveRDS(cloudFreqRaster_cropped, file.path(data.dir, "cloudFreqRaster_cropped.rds"))
saveRDS(annTempRaster_cropped, file.path(data.dir, "annTempRaster_cropped.rds"))
saveRDS(annPrecipRaster_cropped, file.path(data.dir, "annPrecipRaster_cropped.rds"))

test <- readOGR(dsn = file.path(map.dir, "Haw_geo", "Haw_St_geo_20070426_region.shp"), layer = "Haw_St_geo_20070426_region")

?crop
crop(test)
