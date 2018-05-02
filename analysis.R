## ECOLOGICAL ASSEMBLY OF THE HAWAIIAN ARTHROPOD COMMUNITIES
# Author: Jun Ying Lim
# Rarefies OTU tables for analysis

## PACKAGES
library(stringr)
library(plyr)
library(reshape2)
library(picante)
library(betapart)
library(geosphere)

## IMPORT DATA
main.dir <- "~/Dropbox/Projects/2017/hawaiiCommunityAssembly/"
analysis.dir <- file.path(main.dir, "elevationAssemblyHawaii")
data.dir <- file.path(main.dir, "data")
fig.dir <- file.path(main.dir, "figures")
source(file.path(analysis.dir, "ecoDataTools.R"))

# Import otu data
arfSpeciesIDdata <- readRDS(file.path(data.dir, "arfSpeciesIDdata.rds"))
arfZOTUIDdata <- readRDS(file.path(data.dir, "arfZOTUIDdata.rds"))
arfOTUIDdata <- readRDS(file.path(data.dir, "arfOTUIDdata.rds"))

# Import site data
siteData <- read.csv(file.path(data.dir, "clim.final.csv"))

## CALCULATE CLIMATE DISTANCE BETWEEN SITES
rownames(siteData) <- siteData$site.id
climDist <- dist(siteData[c("PC1", "PC2", "PC3")])

## CALCULATE GEOGRAPHIC DISTANCE BETWEEN SITES
siteData$latitude, siteData$longitude
geogDist <- dist(siteData)

distVincentyEllipsoid()

## CALCULATE BETA DIVERSITY BETWEEN SITES
testData <- arfSpeciesIDdata[[1]]

testData_PA <- ifelse(testData > 0, 1, 0)
testData_beta <- beta.pair(testData_PA)
testData_betadist <- testData_beta$beta.sor

## MANTEL TESTS
climDist <- matchDist(testData_betadist, climDist)
geogDist <- matchDist(testData_betadist)
mantel(climDist, testData_betadist)





