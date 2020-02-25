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
# Rarefy datasets
# Use the number of specimen counts in each size category to rarefy
# Samples where there are fewer reads than expected will be left unrarefied (only 17 site by size groupings)
# ARF and MCO are rarefied separately and then combined

nreps <- 100
arf_rarefied <- list()
mco_rarefied <- list()

#arfOTU <- subset(arfOTU, nReads >= 10) # remove OTU with less than 10 reads (possible contamination?)
#mcoOTU <- subset(mcoOTU, nReads >= 10)

arfOTU_summary <- ddply(.data =  arfOTU, .var = .(SizeCategory,Site_ID),.fun = summarize, arfTotalReads = sum(nReads))
mcoOTU_summary <- ddply(.data =  mcoOTU, .var = .(SizeCategory,Site_ID),.fun = summarize, mcoTotalReads = sum(nReads))

arfOTU_summary2 <- merge(arfOTU_summary, specimenCounts)
mcoOTU_summary2 <- merge(mcoOTU_summary, specimenCounts)

sum(mcoOTU_summary2$arfTotalReads < (mcoOTU_summary2$Count * 7)) # all samples have enough reads for rarefaction
sum(mcoOTU_summary2$arfTotalReads < (mcoOTU_summary2$Count * 8)) # all samples have enough reads for rarefaction

sum(arfOTU_summary2$arfTotalReads < (mcoOTU_summary2$Count * 7)) # all samples have enough reads for rarefaction
sum(arfOTU_summary2$arfTotalReads < (mcoOTU_summary2$Count * 8)) # 4 samples have too few reads for rarefaction at the same level as other samples

for(i in 1:nreps){
  pb = txtProgressBar(min = 0, max = nreps, style = 3, initial = 0) 
  setTxtProgressBar(pb, i)
  arf_rarefied[[i]] <- ddply(.data = arfOTU,
                             .var = .(SizeCategory, Site_ID),
                             .fun = rarefyOTUbySpecimenCount,
                             counts = specimenCounts,
                             readsPerIndividual = 7)
                             #.progress = "text",
                             #.parallel = T)
  mco_rarefied[[i]] <- ddply(.data = mcoOTU,
                             .var = .(SizeCategory,Site_ID),
                             .fun = rarefyOTUbySpecimenCount,
                             counts = specimenCounts,
                             readsPerIndividual = 7)
                             #.progress = "text",
                             #.parallel = T)
}
close(pb)

saveRDS(arf_rarefied, file = file.path(data.dir, "arf_rarefied.rds"))
saveRDS(mco_rarefied, file = file.path(data.dir, "mco_rarefied.rds"))

# arf_rarefied <- readRDS(file.path(data.dir, "arf_rarefied.rds"))
# mco_rarefied <- readRDS(file.path(data.dir, "arf_rarefied.rds"))

combinedOTU.dir <- file.path(data.dir, "combinedOTU")
for(i in 1:nreps){
  pb = txtProgressBar(min = 0, max = nreps, style = 3, initial = 0) 
  setTxtProgressBar(pb, i)
  temp <- merge(arf_rarefied[[i]],
                mco_rarefied[[i]], all.x = TRUE)
  combined <- ddply(temp, .variables = .(X.OTU.ID, Site_SizeCategory),
                                  summarize, rarefiedReadAbund = sum(rarefiedReadAbund))
  saveRDS(cleanOTUdata(combined), file.path(combinedOTU.dir, paste0("combinedOTUdata_r", i, ".rds")))
}
close(pb)


## CREATE TAXONOMIC REFERENCE TABLES ===============
# Create reference table mapping species ID, zOTU and OTU
combinedOTU <- rbind(arfOTU[c("Species_ID", "zOTU_ID", "OTU_ID")], mcoOTU[c("Species_ID", "zOTU_ID", "OTU_ID")])
refTable <- unique(combinedOTU)

# Import fasta file to get BLAST-derived taxonomies
seqData <- read.dna(file.path(data.dir, "ARFMCOAll_Species_ZOTU_OTUID042018.fasta"), format = "fasta")
taxonRaw <- rownames(seqData)
taxonTab <- str_split_fixed(taxonRaw, pattern = "_", n = 5)
taxonDF <- as.data.frame(taxonTab)

# Create reference table mapping species ID, zOTU, OTU and taxonomic order
refTableFinal <- merge(refTable,  taxonDF[,1:4], by.x = c("Species_ID", "zOTU_ID", "OTU_ID"), by.y = c("V1", "V2", "V3"))
saveRDS(refTableFinal, file.path(data.dir, "taxonData.rds"))
# NOTE: table will also contain zOTUs that have zero rarefied read abundance

