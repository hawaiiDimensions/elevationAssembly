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
dataMCOmelt <- melt(dataMCO, id.vars = "X.OTU.ID", value.name = "nReads", variable.name = "Site_SizeCategory")

# Clean up data
cleanOTUdata <- function(x){
  # Label OTUs by site and size categories
  temp <- str_split_fixed(x$Site_SizeCategory, pattern = "_", n = 3)
  x$Site_ID <- gsub(temp[,1], pattern = "X", replacement = "")
  x$SizeCategory <- temp[,2]
  temp <- str_split_fixed(x$X.OTU.ID, pattern = "_", n = 3)
  x$Species_ID <- temp[,1]
  x$zOTU_ID <- temp[,2]
  x$OTU_ID <- temp[,3]
  return(x)
}

arfOTU <- cleanOTUdata(dataARFmelt)
mcoOTU <- cleanOTUdata(dataMCOmelt)

## RAREFY DATA BY SIZE CATEGORY AND SITE
# Checking rarefaction procedure
# 17 site by size category groupings have fewer reads than expected
sumReads <- function(df){
  totalReads <- sum(df$nReads)
  rarefiedReads <- nrow(subset(df, nReads > 0)) * 20
    
  #subset(df, Site_ID == unique(df$Site_ID) &
  #                          SizeCategory == unique(df$SizeCategory))$Count * 20
  return(data.frame(totalReads, rarefiedReads))
}

# Determine the number of reads to rarefy from each site and size
test <- ddply(.data = megaData, .var = .(SizeCategory, Site_ID), .fun = sumReads)

test2 <- test[test$totalReads < test$rarefiedReads,]
dim(test2)

# Rarefy datasets
# Samples where there are fewer reads than expected will be left unrarefied (only 17 site by size groupings)
nreps <- 100
arf_rarefied <- list()
mco_rarefied <- list()

for(i in 1:nreps){
  arf_rarefied[[i]] <- ddply(.data = arfOTU,
                             .var = .(SizeCategory,Site_ID),
                             .fun = rarefyOTUbySpecimenCount,
                             counts = specimenCounts,
                             readsPerIndividual = 10,
                             .progress = "text")
  mco_rarefied[[i]] <- ddply(.data = mcoOTU,
                             .var = .(SizeCategory,Site_ID),
                             .fun = rarefyOTUbySpecimenCount,
                             counts = specimenCounts,
                             readsPerIndividual = 10,
                             .progress = "text")
}

# Combine data
megaData <- ddply(merge(arf_rarefied[[1]],
                        mco_rarefied[[1]], all.x=TRUE), 
                  .(X.OTU.ID, Site_SizeCategory), summarise, rarefiedReadAbund=sum(rarefiedReadAbund), .progress = "text")
megaData2 <- cleanOTUdata(megaData)
head(megaData2)
saveRDS(megaData2, file.path(data.dir, "combinedOTUdata.rds"))

# Calculate haplotype diversity per site
calcHaplotypeDiv <- function(df){
  # Calculate haplotype diversity for each site
  #df <- subset(combined_rarefied[[1]], Site_ID == 530 & Species_ID == "Species1")
  #df<- subset(megaData2, Species_ID == "Species100" & Site_ID == "541")
  temp <- subset(df, rarefiedReadAbund > 0)
  nHaplotype <- length(unique(temp$zOTU_ID))
  if(nHaplotype == 0 | nHaplotype == 1){
    haplotypeDiversity <- 0
  } else {
    readCount <- tapply(temp$rarefiedReadAbund, INDEX = temp$zOTU_ID, FUN = sum)
    propAbund <- readCount / sum(readCount)
    haplotypeDiversity <- (nHaplotype / (nHaplotype-1))* (1 - sum(propAbund^2))  
  }
  data.frame(haplotypeDiversity, nHaplotype)
}

test <- ddply(.data = megaData2, .var = .(Species_ID, Site_ID), .fun = calcHaplotypeDiv, .progress = "text")
test2 <- subset(test, nHaplotype > 0)

subset(megaData2, Species_ID == "Species100" & Site_ID == "541")
saveRDS(test2, file.path(data.dir, "haplotypeDiversity.rds"))


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

# Collapse OTUs by OTU delimiters
arfSpeciesID <- lapply(X = arf_rarefied, FUN = collapseOTUsList, collapse.by = "Species_ID")
arfZOTUID <- lapply(X = arf_rarefied, FUN = collapseOTUsList, collapse.by = "zOTU_ID")
arfOTUID <- lapply(X = arf_rarefied, FUN = collapseOTUsList, collapse.by = "OTU_ID")

# Convert all to matrices
arfSpeciesID <- lapply(arfSpeciesID, function(x){acast(x, Site_ID ~ Species_ID, value.var = "readAbund")})
arfZOTUID <- lapply(arfZOTUID, function(x){acast(x, Site_ID ~ zOTU_ID, value.var = "readAbund")})
arfOTUID <- lapply(arfOTUID, function(x){acast(x, Site_ID ~ OTU_ID, value.var = "readAbund")})

## EXPORT RAREFIED MATRICES
saveRDS(arfSpeciesID, file = file.path(data.dir,"arfSpeciesIDdata.rds"))
saveRDS(arfZOTUID, file = file.path(data.dir,"arfZOTUIDdata.rds"))
saveRDS(arfOTUID, file = file.path(data.dir,"arfOTUIDdata.rds"))