#########################################
###### resistanceGA - J. Patiño 10/2021
#########################################

# load libraries
library("devtools")
library("tinytex")
library("ResistanceGA")
library("doParallel")
library("ape")
library("raster") 
library("rgdal")

# set WD and loading climatic variables
setwd("./")
annPrecip_stain <- readRDS("annPrecip_XXX.rds") # replace XXX by the right transect
annTemp_stain <- readRDS("annTemp_XXX.rds") # replace XXX by the right transect

# create my raster object by two climatic variables
r.stack <- stack(annPrecip_stain,annTemp_stain)

# upload and bymatting a beta diversity matrix
XXXOTUPlots.gd.df <- as.matrix(read.csv("betaBrayGeneral_OTUYYY_XXX.csv", row.names = 1)) # replace YYY and XXX by the right taxon and transect respectively
XXXOTUPlots.gen.dist <- as.vector(as.dist(XXXOTUPlots.gd.df, diag=FALSE, upper = FALSE)) # replace XXX by the right transect

# get locations by each sampling plot
XXXOTUPlots.points <- read.csv("coordinates_XXX.csv", sep=",", row.names = 1) # lon, lat; replace XXX by the right transect
XXXOTUPlots.points <- as.matrix(XXXOTUPlots.points) # replace XXX by the right transect
XXXOTUPlots.sample.locales <- SpatialPoints(XXXOTUPlots.points[ ,c(1,2)]) # replace XXX by the right transect
XXXOTUPlots.sample.locales # visualize the object; replace XXX by the right transect

# prepare all the files by lunching the analyses
# replace XXX by the right transect

dir.create(file.path("./","all.comb"))     
XXXOTUPlots.GA.inputs <- GA.prep(ASCII.dir = r.stack, Results.dir = "all.comb",
                         method = "LL",
                         max.cat = 500,
                         max.cont = 500,
                         seed = 555,
                         select.trans = list("A","A"),
                         parallel = 16,
                         maxiter = 1000)

XXXOTUPlots.gdist.inputs <- gdist.prep(length(XXXOTUPlots.sample.locales),
                               samples = XXXOTUPlots.sample.locales,
                               response = XXXOTUPlots.gen.dist,
                               method = 'commuteDistance') ## optimize using commute distance
                               
# run the final optimization
# replace XXX by the right transect

XXXOTUPlots.allComb_RESULTS.gdist <- all_comb(gdist.inputs = XXXOTUPlots.gdist.inputs,
                             GA.inputs = XXXOTUPlots.GA.inputs,
                             results.dir = "all.comb/",
                             max.combination = 2,
                             sample.prop = 0.75,
                             replicate = 2
                              )

# save all the outputs and replace replace XXX, YYY and ZZZ by the right taxon, transect and threshold respectively
save.image("XXX-YYYOTUZZZ.allComb.resistanceGA.RData")
