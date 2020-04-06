## ECOLOGICAL ASSEMBLY OF THE HAWAIIAN ARTHROPOD COMMUNITIES
# Author: Jun Ying Lim
# Rarefies OTU tables for analysis

## PACKAGES ===============
library(stringr)
library(plyr)
library(reshape2)
library(picante)

## IMPORT DATA ===============
main.dir <- "~/Dropbox/projects/2017/hawaiiCommunityAssembly/elevationAssemblyHawaii"
data.dir <- file.path(main.dir, "data")
fig.dir <- file.path(main.dir, "figures")
source(file.path(main.dir, "metabarcodingTools.R"))

# SUBSET THE DATAFRAME ====================
# Import taxonomic reference table
taxonData <- readRDS(file.path(data.dir, "taxonData.rds"))

# Define some output folders
masterOTU.dir <- file.path(data.dir, "masterOTU")
OTU_native.dir <- file.path(data.dir, "OTU_native")
zOTU_native.dir <- file.path(data.dir, "zOTU_native")
orderList <- c("Araneae", "Hemiptera", "Lepidoptera", "Psocoptera", "Coleoptera", "Orthoptera")
OTU_native_order.dir <- file.path(file.path(data.dir, "OTU_native_order"))


otufiles <- list.files(data.dir)[grep(list.files(data.dir), pattern = "combinedOTUdata_r")]
for(i in 1:100){
  # Create master zOTU table
  combData <- readRDS(file.path(data.dir, otufiles[i]))
  masterTable <- acast(zOTU_ID ~Site_ID, data = combData, fun.aggregate = sum, value.var = "rarefiedReadAbund")
  saveRDS(masterTable, file.path(masterOTU.dir, paste0("master_zOTU_r", i, ".rds")))
  
  # Collapse ZOTU matrix into OTU and Species matrices
  taxonData <- taxonData[taxonData$zOTU_ID %in% combData$zOTU_ID,]
  OTUtab <- collapseOTU(x = masterTable, ref = taxonData, collapse.by = "OTU_ID", to.collapse = "zOTU_ID")#, subset.by = "V4")
  
  # Exclude the following orders 
  toExclude <- subset(taxonData, V4 %in% c("Collembola", "Chilopoda", "Diplopoda", "Blattodea", "Isopoda", "Mantodea", "Phasmatodea"))
  masterTable_native <- masterTable[!rownames(masterTable) %in% toExclude$zOTU_ID, ]
  saveRDS(masterTable_native, file.path(zOTU_native.dir, paste0("zOTU_native_r", i, ".rds")))
  
  OTUtab_native <- OTUtab[!rownames(OTUtab) %in% toExclude$OTU_ID, ]
  saveRDS(OTUtab_native, file.path(OTU_native.dir, paste0("OTU_native_r", i, ".rds")))
  
  # Subset individual orders
  for(j in orderList){
    targetTaxon <- subset(taxonData, V4 == j)
    masterTable_subset <- masterTable_native[rownames(masterTable_native) %in% targetTaxon$zOTU_ID, ]
    OTUtab_subset <- OTUtab_native[rownames(OTUtab_native) %in% targetTaxon$OTU_ID, ]
    saveRDS(masterTable_subset, file.path(OTU_native_order.dir, paste0("zOTU_", j, "_r", i, ".rds")))
    saveRDS(OTUtab_subset, file.path(OTU_native_order.dir, paste0("OTU_", j, "_r", i,".rds")))
  }
}

# write.csv(masterTable, file.path(data.dir, "masterZOTU.csv"))
# write.csv(OTUtab, file.path(data.dir, "OTUtab.csv"))
# write.csv(SPPtab, file.path(data.dir, "SPPtab.csv"))
# write.csv(taxonData, file.path(data.dir, "taxonData.csv"))

# select random rarefaction
set.seed(12345)
otu_subsets <- sample(1:100, size = 5, replace = FALSE)
# 14, 51, 80, 90, 92

# split into two different sites for Jairo
resistance.dir <- file.path(data.dir, "resistance_analysis")
jairofname <- list.files(resistance.dir)
jairotab <- lapply(file.path(resistance.dir, jairofname), FUN = readRDS)

laupahoehoe_siteIDs <-subset(siteData, site == "Laupahoehoe")$site.id
stainback_siteIDs <-subset(siteData, site == "Stainback")$site.id

#laupahoehoe_tab <- list()
#stainback_tab <- list()
for(i in 1:length(jairotab)){
  #laupahoehoe_tab[[i]] <-
  saveRDS(jairotab[[i]][,laupahoehoe_siteIDs], file = file.path(resistance.dir, paste0("laupahoehoe_", jairofname[i])))
  #stainback_tab[[i]] <-
  saveRDS(jairotab[[i]][,stainback_siteIDs], file = file.path(resistance.dir, paste0("stainback_", jairofname[i])))
}