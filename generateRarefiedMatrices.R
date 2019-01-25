## ECOLOGICAL ASSEMBLY OF THE HAWAIIAN ARTHROPOD COMMUNITIES
# Author: Jun Ying Lim
# Rarefies OTU tables for analysis

## PACKAGES ===============
library(stringr)
library(plyr)
library(reshape2)
library(picante)

## IMPORT DATA ===============
main.dir <- "~/Dropbox/Projects/2017/hawaiiCommunityAssembly/"
analysis.dir <- file.path(main.dir, "elevationAssemblyHawaii")
data.dir <- file.path(main.dir, "data")
fig.dir <- file.path(main.dir, "figures")
source(file.path(analysis.dir, "metabarcodingTools.R"))

# Import site and otu tables
siteData <- read.csv(file.path(data.dir, "clim.final.csv"))
dataARF <- read.delim(file.path(data.dir, "otutabARF.txt"), header = TRUE, stringsAsFactors = FALSE)
dataMCO <- read.delim(file.path(data.dir, "otutabMCO.txt"), header = TRUE, stringsAsFactors = FALSE)

# Import specimen count table
specimenCounts <- read.delim(file.path(data.dir, "specimennumber.txt"), header = FALSE)
names(specimenCounts) <- c("Site_ID", "Transect", "Island", "Method", "Year", "SizeCategory", "Count", "RarefiedReadAbund")

## DATA PREPARATION ===============
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

## RAREFY DATA BY SIZE CATEGORY AND SITE ===============
# Checking rarefaction procedure
# 17 site by size category groupings have fewer reads than expected
# sumReads <- function(df){
#   totalReads <- sum(df$nReads)
#   rarefiedReads <- nrow(subset(df, nReads > 0)) * 20
#     
#   #subset(df, Site_ID == unique(df$Site_ID) &
#   #                          SizeCategory == unique(df$SizeCategory))$Count * 20
#   return(data.frame(totalReads, rarefiedReads))
# }
# 
# # Determine the number of reads to rarefy from each site and size
# test <- ddply(.data = megaData, .var = .(SizeCategory, Site_ID), .fun = sumReads)
# test2 <- test[test$totalReads < test$rarefiedReads,]

# Rarefy datasets
# Samples where there are fewer reads than expected will be left unrarefied (only 17 site by size groupings)
# ARF and MCO are rarefied separately and then combined
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
combData <- ddply(merge(arf_rarefied[[1]],
                        mco_rarefied[[1]], all.x=TRUE), 
                  .variables = .(X.OTU.ID, Site_SizeCategory),
                  summarise, rarefiedReadAbund=sum(rarefiedReadAbund), .progress = "text")
combData_clean <- cleanOTUdata(combData)
head(combData_clean)
saveRDS(combData_clean, file.path(data.dir, "combinedOTUdata.rds"))

## CREATE TAXONOMIC REFERENCE TABLES ===============
combData_clean <- readRDS(file.path(data.dir, "combinedOTUdata.rds"))

# Create reference table mapping species ID, zOTU and OTU
refTable <- combData_clean[c("Species_ID", "zOTU_ID", "OTU_ID")]
refTable <- unique(refTable)

# Import fasta file to get BLAST-derived taxonomies
seqData <- read.dna(file.path(data.dir, "ARFMCOAll_Species_ZOTU_OTUID042018.fasta"), format = "fasta")
taxonRaw <- rownames(seqData)
taxonTab <- str_split_fixed(taxonRaw, pattern = "_", n = 5)
taxonDF <- as.data.frame(taxonTab)

# Create reference table mapping species ID, zOTU, OTU and taxonomic order
refTableFinal <- merge(refTable,  taxonDF[,1:4], by.x = c("Species_ID", "zOTU_ID", "OTU_ID"), by.y = c("V1", "V2", "V3"))
saveRDS(refTableFinal, file.path(data.dir, "taxonData.rds"))

# NOTE: tables will only contain zOTUs that have above zero

# SUBSET THE DATAFRAME ====================
# Import taxonomic reference table
source(file.path(analysis.dir, "metabarcodingTools.R"))
taxonData <- readRDS(file.path(data.dir, "taxonData.rds"))

# Create master zOTU matrix
masterZOTU <- acast(zOTU_ID~Site_ID, data = combData_clean, fun.aggregate = sum, value.var = "rarefiedReadAbund")
saveRDS(masterZOTU, file.path(data.dir, "masterZOTU.rds"))
masterTable <- readRDS(file.path(data.dir, "masterZOTU.rds"))

# Collapse ZOTU matrix into OTU and Species matrices
OTUtab <- collapseOTU(x = masterTable, ref = taxonData, collapse.by = "OTU_ID", to.collapse = "zOTU_ID")#, subset.by = "V4")
SPPtab <- collapseOTU(x = masterTable, ref = taxonData, collapse.by = "Species_ID", to.collapse = "zOTU_ID")#, subset.by = "V4")

# Exclude the following orders
toExclude <- subset(taxonData, V4 %in% c("Collembola", "Chilopoda", "Diplopoda", "Blattodea", "Isopoda", "Mantodea", "Phasmatodea"))

masterTable_native <- masterTable[!rownames(masterTable) %in% toExclude$zOTU_ID, ]
saveRDS(masterTable_native, file.path(data.dir, "zOTU_native.rds"))

OTUtab_native <- OTUtab[!rownames(OTUtab) %in% toExclude$OTU_ID, ]
saveRDS(OTUtab_native, file.path(data.dir, "OTU_native.rds"))
SPPtab_native <- SPPtab[!rownames(SPPtab) %in% toExclude$Species_ID, ]
saveRDS(SPPtab_native, file.path(data.dir, "SPP_native.rds"))

orderList <- c("Araneae", "Hemiptera", "Lepidoptera", "Psocoptera", "Coleoptera", "Orthoptera")

for(i in orderList){
  targetTaxon <- subset(taxonData, V4 == i)
  masterTable_subset <- masterTable_native[rownames(masterTable_native) %in% targetTaxon$zOTU_ID, ]
  OTUtab_subset <- OTUtab_native[rownames(OTUtab_native) %in% targetTaxon$OTU_ID, ]
  SPPtab_subset <- SPPtab_native[rownames(SPPtab_native) %in% targetTaxon$Species_ID, ]
  saveRDS(masterTable_subset, file.path(data.dir, paste0("zOTU_", i, ".rds")))
  saveRDS(OTUtab_subset, file.path(data.dir, paste0("OTU_", i, ".rds")))
  saveRDS(SPPtab_subset, file.path(data.dir, paste0("SPP_", i, ".rds")))
}

# write.csv(masterTable, file.path(data.dir, "masterZOTU.csv"))
# write.csv(OTUtab, file.path(data.dir, "OTUtab.csv"))
# write.csv(SPPtab, file.path(data.dir, "SPPtab.csv"))
# write.csv(taxonData, file.path(data.dir, "taxonData.csv"))

# 


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

haplotDivBySiteByOTU <- ddply(.data = subset(combData_clean, (!zOTU_ID %in% toExclude$zOTU_ID)), .var = .(OTU_ID, Site_ID), .fun = calcHaplotypeDiv, .progress = "text")

saveRDS(haplotDivBySiteByOTU, file.path(data.dir, "haplotDivBySiteByOTU.rds"))

ReadAbundanceBySiteByOTU <- ddply(.data = subset(combData_clean, (!zOTU_ID %in% toExclude$zOTU_ID)), .var = .(OTU_ID, Site_ID), .fun = summarize, siteReadAbundance = sum(rarefiedReadAbund), .progress = "text")

saveRDS(ReadAbundanceBySiteByOTU, file.path(data.dir, "ReadAbundanceBySiteByOTU.rds"))

## PRELIMINARY PLOTS
ReadAbundanceBySiteByOTU_Final <- merge(ReadAbundanceBySiteByOTU, taxonData, by = "OTU_ID")

temp <- tapply(ReadAbundanceBySiteByOTU_Final$siteReadAbundance, ReadAbundanceBySiteByOTU_Final$OTU_ID, FUN = max)
temp2 <- data.frame("OTU_ID" = names(temp), "maxAbund" = as.vector(temp))
ReadAbundanceBySiteByOTU_Final2 <- merge(ReadAbundanceBySiteByOTU_Final, temp2, by = "OTU_ID")

z <- merge(ReadAbundanceBySiteByOTU_Final2, siteData[c("elevation", "site", "rf_ann", "t_ann", "site.id")], by.x = "Site_ID", by.y = "site.id")
z$relAbund <- z$siteReadAbundance/ z$maxAbund
z2 <- ggplot(aes(x = elevation, y = relAbund), data = subset(z, maxAbund > 1000)) + geom_point() + facet_wrap(vars(OTU_ID,site)) + geom_smooth(method = "loess")
ggsave("~/Desktop/readAbundElev.pdf", z2, width = 20, height = 20)



z3 <- merge(haplotDivBySiteByOTU, siteData[c("elevation", "site", "rf_ann", "t_ann", "site.id")], by.x = "Site_ID", by.y = "site.id")
z4 <- subset(z3, nHaplotype > 0)
z4 <- merge(z4, taxonData)
#z5 <- ddply(.data = z4, .variables = .(Site_ID), summarize, meanHaplo = mean(haplotypeDiversity))
#z6 <- merge(z4, siteData[c("elevation", "site", "rf_ann", "t_ann", "site.id")], by.x = "Site_ID", by.y = "site.id")

ggplot(aes(y = haplotypeDiversity, x = elevation, color = OTU_ID), data = z4) + facet_wrap(~site) + theme(legend.position = "none") + geom_smooth(method="glm", se = FALSE, alpha =0.5)
library(lme4)
mod <- glmer(haplotypeDiversity ~ elevation + (1|OTU_ID) + (elevation|OTU_ID), data = z4, family = "binomial")
summary(mod)
library(ggplot2)

install.packages("mcmcglmm")

library(nimble)
haploCode <- nimbleCode({
  for(i in 1:OTU_ID){
    
  }
})
table(subset(z4, site == "Laupahoehoe")$OTU_ID)

subset(z4, site == "Laupahoehoe" & OTU_ID == "OTU90")
# There are 25 laupahoehoe and 39 Stainback sites
