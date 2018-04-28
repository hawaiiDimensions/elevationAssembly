## ECOLOGICAL ASSEMBLY OF THE HAWAIIAN ARTHROPOD COMMUNITIES
# Author: Jun Ying Lim
# Created: 13th April 2018

## PACKAGES
library(stringr)
library(plyr)
library(reshape2)
library(picante)

## IMPORT DATA
main.dir <- "~/Dropbox/Projects/2017/hawaiiCommunityAssembly/"
data.dir <- file.path(main.dir, "data")
fig.dir <- file.path(main.dir, "figures")
source(file.path(main.dir, "ecoDataTools.R"))

# Import site and otu tables
siteData <- read.csv(file.path(data.dir, "clim.final.csv"))
speciesData <- otuData[,1:4]

dataARF <- read.delim(file.path(data.dir, "otutabARF.txt"), header = TRUE, stringsAsFactors = FALSE)
dataMCO <- read.delim(file.path(data.dir, "otutabMCO.txt"), header = TRUE, stringsAsFactors = FALSE)

# Collapse data by OTUs and zOTUs
temp <- str_split_fixed(dataARF$X.OTU.ID, n = 3, pattern = "_")
dataARF$SpeciesID <- temp[,1]
dataARF$ZOTU_ID <- temp[,2]
dataARF$OTU_ID <- temp[,3]

targetCol<- names(dataARF)[grep(names(dataARF), pattern = "X[1-9]")]

collapseData <- function(df){
  return(colSums(df[targetCol]))
}

otuARF <- ddply(.data = dataARF, .fun = collapseData, .variables = .(OTU_ID))
zotuARF <- ddply(.data = dataARF, .fun = collapseData, .variables = .(OTU_ID))

# Collapse the OTU and ZOTU tables by site, across size cateogries
zotuARF_melt <- melt(zotuARF, id.vars = "OTU_ID")
temp <- str_split_fixed(zotuARF_melt$variable, n = 3, pattern = "_")
zotuARF_melt$site_ID <- temp[,1]
zotuARF_melt$siteCategory <- temp[,2]

collapseSites <- function(df){
  return(data.frame("nreads" = sum(df$value)))
}

zotuARF_site <- ddply(.data = zotuARF_melt, .fun = collapseSites, .variables = .(site_ID, OTU_ID), progress = "text")
zotuARF_site_mat <- acast(zotuARF_site, site_ID ~ OTU_ID, value.var = "nreads")

rownames(zotuARF_site_mat) <- gsub(rownames(zotuARF_site_mat), pattern = "X", replacement = "")

# Calculate number of OTUs per sample
zotuARF_site_mat_pa <- ifelse(zotuARF_site_mat >0, 1, 0) # change to presence/absence
otuSite <- rowSums(zotuARF_site_mat_pa) # calculate no. of OTUs per site

readRandomize <- otuSite*100

sample(zotuARF_site_mat)

randomizeMatrix(zotuARF_site_mat, )
sample(zotuARF_site_mat)

x <- rrarefy(zotuARF_site_mat[1,], readRandomize[1])
x[1:5]
zotuARF_site_mat[1,1:5]
rowSums(zotuARF_site_mat_pa)
rowSums(x)


# Randomize otuARF


# Analyze OTUs
head(otuARF_melt)

