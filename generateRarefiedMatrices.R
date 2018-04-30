## ECOLOGICAL ASSEMBLY OF THE HAWAIIAN ARTHROPOD COMMUNITIES
# Author: Jun Ying Lim
# Rarefies OTU tables for analysis

## PACKAGES
library(stringr)
library(plyr)
library(reshape2)
library(picante)

## IMPORT DATA
main.dir <- "~/Dropbox/Projects/2017/hawaiiCommunityAssembly/"
analysis.dir <- file.path(main.dir, "elevationAssemblyHawaii")
data.dir <- file.path(main.dir, "data")
fig.dir <- file.path(main.dir, "figures")
source(file.path(analysis.dir, "ecoDataTools.R"))

# Import site and otu tables
siteData <- read.csv(file.path(data.dir, "clim.final.csv"))
dataARF <- read.delim(file.path(data.dir, "otutabARF.txt"), header = TRUE, stringsAsFactors = FALSE)
dataMCO <- read.delim(file.path(data.dir, "otutabMCO.txt"), header = TRUE, stringsAsFactors = FALSE)

# Import specimen count table
specimenCounts <- read.delim(file.path(data.dir, "specimennumber.txt"), header = FALSE)
names(specimenCounts) <- c("Site_ID", "Transect", "Island", "Method", "Year", "SizeCategory", "Count", "RarefiedReadAbund")

## DATA PREPARATION
# Convert matrix into long-form data frame
dataARFmelt <- melt(dataARF, id.vars = "X.OTU.ID", value.name = "nReads", variable.name = "Site_SizeCategory")

# Label OTUs by site and size categories
temp <- str_split_fixed(dataARFmelt$Site_SizeCategory, pattern = "_", n = 3)
dataARFmelt$Site_ID <- gsub(temp[,1], pattern = "X", replacement = "")
dataARFmelt$SizeCategory <- temp[,2]

# Label reads by various species delimitation methods
temp <- str_split_fixed(dataARFmelt$X.OTU.ID, pattern = "_", n = 3)
dataARFmelt$Species_ID <- temp[,1]
dataARFmelt$zOTU_ID <- temp[,2]
dataARFmelt$OTU_ID <- temp[,3]

## RAREFY DATA BY SIZE CATEGORY AND SITE
# Checking rarefaction procedure
# 17 site by size category groupings have fewer reads than expected
sumReads <- function(df){
  totalReads <- sum(df$nReads)
  rarefiedReads <- subset(counts, Site_ID == unique(df$Site_ID) &
                            SizeCategory == unique(df$SizeCategory))$Count * 10
  return(data.frame(totalReads, rarefiedReads))
}
test <- ddply(.data = dataARFmelt, .var = .(SizeCategory, Site_ID), .fun = sumReads)
test2 <- test[test$totalReads < test$rarefiedReads,]
dim(test2)

# Rarefy datasets
# Samples where there are fewer reads than expected will be left unrarefied (only 17 site by size groupings)
nreps <- 100
arf_rarefied <- list()

for(i in 1:nreps){
  arf_rarefied[[i]] <- ddply(.data = dataARFmelt,
                             .var = .(SizeCategory,Site_ID),
                             .fun = rarefyOTUbySpecimenCount,
                             counts = specimenCounts,
                             readsPerIndividual = 10,
                             .progress = "text")
}

# Collapse data into OTUs
collapseOTUs <- function(df){
  # Collapses rarefied read abundance data frame by values in target column
  # Args:
  #   df = data.frame object
  #   collapse.by = string, name of column to collapse by
  # Returns:
  #   df = data.frame
  readAbund <- sum(df$rarefiedReadAbund)
  return(data.frame(readAbund))
}

collapseOTUsList <- function(x, collapse.by){
  # Wrapper function to iterate ddply function across each rarefied data frame using lapply
  ddply(.data = x, .fun = collapseOTUs, .variables = c("Site_ID", collapse.by), .progress = "text")
}

arfSpeciesID <- lapply(X = arf_rarefied, FUN = collapseOTUsList, collapse.by = "Species_ID")
arfZOTUID <- lapply(X = arf_rarefied, FUN = collapseOTUsList, collapse.by = "zOTU_ID")
arfOTUID <- lapply(X = arf_rarefied, FUN = collapseOTUsList, collapse.by = "OTU_ID")
head(arf_rarefied[[1]])

## EXPORT RAREFIED MATRICES
saveRDS(arfSpeciesID, file = file.path(data.dir,"arfSpeciesIDdata.rds"))
saveRDS(arfZOTUID, file = file.path(data.dir,"arfZOTUIDdata.rds"))
saveRDS(arfOTUID, file = file.path(data.dir,"arfOTUIDdata.rds"))