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
library(vegan)

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

climPCA <- prcomp(siteData[c(paste0("BIO", 1:19))], center = TRUE, scale. = TRUE)
head(climPCA$x)
head(climPCA$rotation)
climPCA$sdev / sum(climPCA$sdev) * 100
which.max(as.data.frame(climPCA$rotation)$PC1)
as.data.frame(climPCA$rotation)$PC2

climPCAcoord <- as.data.frame(climPCA$x)
climPCAcoord$Site_ID <- rownames(climPCAcoord)
library(ggplot2)
ggplot(data = climPCAcoord) + geom_point(aes(y = PC1, x = PC2)) + geom_text_repel(aes(y = PC1, x = PC2, label = Site_ID))
library(ggrepel)
head(climPCAcoord)

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
#geogDist <- matchDist(testData_betadist)
mantel(climDist, testData_betadist)

## Mantel tests;
# * climate distance vs. turnover
# * geographic distance vs. turnover

# Plot abundance against PC1, PC2 (for groups that are found in more than 5 sites, in each transect)
# Plot abundance against PC1, PC2 for groups that are found in both sites, highlight by site
# Find the mean? Niche distance? Schoener's D
# Fit a truncated normal distribution? 

# Remove collembolla
# 




