## ECOLOGICAL ASSEMBLY OF THE HAWAIIAN ARTHROPOD COMMUNITIES
# Author: Jun Ying Lim
# Rarefies OTU tables for analysis

## PACKAGES ============
library(stringr); library(plyr); library(reshape2) # data manipulation tools
library(vegan); library(betapart) # calculating beta diversity
library(geosphere) # calculating geographic distances
library(ggplot2); library(ggrepel)

## IMPORT DATA ============
main.dir <- "~/Dropbox/Projects/2017/hawaiiCommunityAssembly/"
analysis.dir <- file.path(main.dir, "elevationAssemblyHawaii")
data.dir <- file.path(main.dir, "data")
fig.dir <- file.path(analysis.dir, "figures")
source(file.path(analysis.dir, "ecoDataTools.R"))

# Import otu data
arfSpeciesIDdata <- readRDS(file.path(data.dir, "arfSpeciesIDdata.rds"))
arfZOTUIDdata <- readRDS(file.path(data.dir, "arfZOTUIDdata.rds"))
arfOTUIDdata <- readRDS(file.path(data.dir, "arfOTUIDdata.rds"))

# Import site data
siteData <- read.csv(file.path(data.dir, "clim.final.csv"), stringsAsFactors = FALSE)
siteData <- siteData[-grep(siteData$site.id, pattern = "BRG"),] # Exclude Rosie's samples

laupahoehoe_siteIDs <- subset(siteData, site.1 == "Laupahoehoe")$site.id
steinbeck_siteIDs <- subset(siteData, site.1 == "Stainback")$site.id

## CALCULATE HAPLOTYPE DIVERSITY
test <- arfZOTUIDdata[[1]]
test2 <- arfOTUIDdata[[1]]



## CALCULATE CLIMATE DISTANCE BETWEEN SITES ============
# Principal coordinate analysis of bioclim variables
rownames(siteData) <- siteData$site.id
climPCA <- prcomp(siteData[c(paste0("BIO", 1:19))], center = TRUE, scale. = TRUE)
climPCA$sdev / sum(climPCA$sdev) * 100 # variance explained

siteData$PC1 <- climPCA$x[,1]
siteData$PC2 <- climPCA$x[,2]
siteData$PC3 <- climPCA$x[,3]

# Plot sites in climate space
climPCAcoord <- as.data.frame(climPCA$x)
climPCAcoord$Site_ID <- rownames(climPCAcoord)
#ggplot(data = climPCAcoord) + geom_point(aes(y = PC1, x = PC2)) + geom_text_repel(aes(y = PC1, x = PC2, label = Site_ID))

# Calculate climatic distance between sites
climDist <- dist(siteData[c("PC1", "PC2", "PC3")])

## CALCULATE GEOGRAPHIC DISTANCE BETWEEN SITES
nSite <- length(siteData$site.id)
geogDist <- matrix(data = NA, nrow = nSite, ncol = nSite)

for(i in 1:nSite){
  for(j in 1:nSite){
    geogDist[i,j] <- distVincentyEllipsoid(p1 = c(siteData$longitude[i], siteData$latitude[i]),
                                           p2 = c(siteData$longitude[j], siteData$latitude[j]))
  }
}
rownames(geogDist) <- siteData$site.id
colnames(geogDist) <- siteData$site.id

## CALCULATE BETA DIVERSITY BETWEEN SITES
# Try this out with one site first
testData <- arfSpeciesIDdata[[1]]

testData_PA <- ifelse(testData > 0, 1, 0)
testData_beta <- beta.pair(testData_PA)
testData_betadist <- testData_beta$beta.sor

## MANTEL TESTS
climDist <- matchDist(testData_betadist, climDist)
geogDist <- matchDist(testData_betadist, geogDist)

climDist_steinbeck <- as.matrix(climDist)[steinbeck_siteIDs, steinbeck_siteIDs]
climDist_steinbeck_vector <- as.vector(as.dist(climDist_steinbeck))
geogDist_steinbeck <- as.matrix(geogDist)[steinbeck_siteIDs, steinbeck_siteIDs]
geogDist_steinbeck_vector <- as.vector(as.dist(geogDist_steinbeck))

climDist_laup <- as.matrix(climDist)[laupahoehoe_siteIDs, laupahoehoe_siteIDs]
climDist_laup_vector <- as.vector(as.dist(climDist_laup))
geogDist_laup <- as.matrix(climDist)[laupahoehoe_siteIDs, laupahoehoe_siteIDs]
geogDist_laup_vector <- as.vector(as.dist(geogDist_laup))

betaDist_steinbeck <- as.matrix(testData_betadist)[steinbeck_siteIDs, steinbeck_siteIDs]
betaDist_steinbeck_vector <- as.vector(as.dist(betaDist_steinbeck))
betaDist_laup <- as.matrix(testData_betadist)[laupahoehoe_siteIDs, laupahoehoe_siteIDs]
betaDist_laup_vector <- as.vector(as.dist(betaDist_laup))

plot(betaDist_steinbeck_vector ~ climDist_steinbeck_vector)
plot(betaDist_steinbeck_vector ~ geogDist_steinbeck_vector)

mantel(geogDist_steinbeck, betaDist_steinbeck)
mantel(geogDist_laup, betaDist_laup)

mantel(climDist_steinbeck, betaDist_steinbeck)
mantel(climDist_laup, betaDist_laup)

# Scaled to unit variance
geogDist_steinbeck_scaled <- geogDist_steinbeck_vector / sd(geogDist_steinbeck_vector)
climDist_steinbeck_scaled <- climDist_steinbeck_vector / sd(climDist_steinbeck_vector)

mod1 <- lm(betaDist_steinbeck_vector ~ geogDist_steinbeck_scaled + climDist_steinbeck_scaled)
summary(mod1)
# maybe geographic distance may not be the best; since all the climatic variation is collinear with distance (that's why it's an elevational transect!) What you would need is randomized sites that have varying geographic distance and climatic variation.

## SPECIES ABUNDANCE WITH CLIMATE
testData_melt <- melt(testData, value.name = "nReads", varnames = c("Site_ID", "Species_ID"))
testData_melt <- merge(testData_melt, siteData[c("site.id", "site.1", "PC1", "PC2")], by.x = "Site_ID", by.y= "site.id")

cleanData <- function(x, nSites){
  # Keep only species that are found in more than equal the specified number of sites
  if(sum(x$nReads > 0) >= nSites){
    return(x)
  } else {
    return(NULL)
  }
}

table(siteData$site.1) # 25 laupahoehoe sites; 39 steinbeck
rtest <- ddply(.data = subset(testData_melt, site.1 == "Laupahoehoe"),
               .fun = cleanData,
               .variables = .(Species_ID),
               nSites = 16)

speciesAbundancePlot <- ggplot(data = rtest) +
  geom_point(aes(y = log10(nReads+1), x = PC1, color = Species_ID)) +
  geom_smooth(aes(y = log10(nReads+1), x = PC1, color = Species_ID), se = FALSE) +
  facet_wrap(~Species_ID) +
  theme(legend.position = "none")

ggsave(speciesAbundancePlot,
       filename = file.path(fig.dir, "speciesAbundancePlot_laupahoehoe.pdf"),
       height = 10, width = 10)

rtest2 <- ddply(.data = subset(testData_melt, site.1 == "Stainback"),
               .fun = cleanData,
               .variables = .(Species_ID),
               nSites = 24)

speciesAbundancePlot_steinbeck <- ggplot(data = rtest2) +
  geom_point(aes(y = log10(nReads+1), x = PC1, color = Species_ID)) +
  geom_smooth(aes(y = log10(nReads+1), x = PC1, color = Species_ID), se = FALSE) +
  facet_wrap(~Species_ID) +
  theme(legend.position = "none")

ggsave(speciesAbundancePlot_steinbeck,
       filename = file.path(fig.dir, "speciesAbundancePlot_steinbeck.pdf"),
       height = 10, width = 10)


## Mantel tests;
# * climate distance vs. turnover
# * geographic distance vs. turnover

# Import fasta (still need genetic distance from henrik)
# Plot abundance against PC1, PC2 (for groups that are found in more than 5 sites, in each transect)
# Plot abundance against PC1, PC2 for groups that are found in both sites, highlight by site
# Find the mean? Niche distance? Schoener's D
# Fit a truncated normal distribution? 

# Remove collembolla
# 

# Calculate genetic distance
library(ape)
otuFasta <- read.dna(file.path(data.dir, "ARFMCOAll_Species_ZOTU_OTUID042018.fasta"), format = "fasta")
test <- rownames(otuFasta)
test <- gsub(pattern = "_size_.*", replacement = "", test)
temp <- str_split_fixed(test, pattern = "_", n = 6)

taxoRef <- as.data.frame(temp)

collapseOtuRefs <- function(x){
  #x = subset(taxoRef, V1 == "Species1008")
  Order <- names(which.max(table(x$V4)))
  Family <- names(which.max(table(x$V5)))
  Genus <- names(which.max(table(x$V6)))
  data.frame(Order, Family, Genus)
}

taxoRefSpecies <- ddply(.fun = collapseOtuRefs, .variables = .(V1), .data = taxoRef)
taxoRefZOTUs <- ddply(.fun = collapseOtuRefs, .variables = .(V2), .data = taxoRef)
taxoRefOTUs <- ddply(.fun = collapseOtuRefs, .variables = .(V3), .data = taxoRef)

write.csv(taxoRefSpecies, file = file.path(data.dir, "taxoRefSpecies.csv"))
write.csv(taxoRefZOTUs, file = file.path(data.dir, "taxoRefZOTUs.csv"))
write.csv(taxoRefOTUs, file = file.path(data.dir, "taxoRefOTUs.csv"))

## CALCULATE GENETIC DISTANCE BETWEEN OTUS ================

geneticDist <- dist.dna(otuFasta, model = "K80", as.matrix = TRUE) # Genetic distances based on kimura 2-rate parameter
geneticDist_df <- melt(geneticDist)

geneticDist_df$Species1 <- str_split_fixed(geneticDist_df$Var1, n = 5, pattern = "_")[,1]
geneticDist_df$Species2 <- str_split_fixed(geneticDist_df$Var2, n = 5, pattern = "_")[,1]

geneticDist_df$ZOTU1 <- str_split_fixed(geneticDist_df$Var1, n = 5, pattern = "_")[,2]
geneticDist_df$ZOTU2 <- str_split_fixed(geneticDist_df$Var2, n = 5, pattern = "_")[,2]

geneticDist_df$OTU1 <- str_split_fixed(geneticDist_df$Var1, n = 5, pattern = "_")[,3]
geneticDist_df$OTU2 <- str_split_fixed(geneticDist_df$Var2, n = 5, pattern = "_")[,3]


geneticDist_species <- ddply(.data = geneticDist_df, .var = .(Species1, Species2), .fun = function(x){ geneticDist <- mean(x$value); data.frame(geneticDist) }, .progress = "text")
geneticDist_ZOTU <- ddply(.data = geneticDist_df, .var = .(ZOTU1, ZOTU2), .fun = function(x){ geneticDist <- mean(x$value); data.frame(geneticDist) }, .progress = "text")
geneticDist_OTU <- ddply(.data = geneticDist_df, .var = .(OTU1, OTU2), .fun = function(x){ geneticDist <- mean(x$value); data.frame(geneticDist) }, .progress = "text")

saveRDS(geneticDist_species, file.path(data.dir, "geneticDist_species.rds"))
saveRDS(geneticDist_ZOTU, file.path(data.dir, "geneticDist_ZOTU.rds"))
saveRDS(geneticDist_OTU, file.path(data.dir, "geneticDist_OTU.rds"))

length(unique(geneticDist_df$Species1))
length(unique(geneticDist_df$ZOTU1))
length(unique(geneticDist_df$OTU1))

geneticDist_OTU <- readRDS(file.path(data.dir, "geneticDist_OTU.rds"))
geneticDist_species <- readRDS(file.path(data.dir, "geneticDist_species.rds"))
geneticDist_ZOTU <- readRDS(file.path(data.dir, "geneticDist_ZOTU.rds"))

geneticDist_OTU_mat <- acast(geneticDist_OTU, OTU1 ~ OTU2, value.var = "geneticDist")
geneticDist_species_mat <- acast(geneticDist_species, Species1 ~ Species2, value.var = "geneticDist")
geneticDist_ZOTU_mat <- acast(geneticDist_ZOTU, ZOTU1 ~ ZOTU2, value.var = "geneticDist")

saveRDS(geneticDist_species_mat, file.path(data.dir, "geneticDist_species_mat.rds"))
saveRDS(geneticDist_OTU_mat, file.path(data.dir, "geneticDist_OTU_mat.rds"))
saveRDS(geneticDist_ZOTU_mat, file.path(data.dir, "geneticDist_ZOTU_mat.rds"))

## CALCULATE PAIRWISE SITE GENETIC DIVERSITY ================
library(picante)

geneticDist_ZOTU_mat <- readRDS(file.path(data.dir, "geneticDist_ZOTU_mat.rds"))
arfZOTU_geneticDist_laup <- comdist(comm = arfZOTU_laup, dis = geneticDist_ZOTU_mat, abundance.weighted = FALSE)
arfZOTU_geneticDist_stein <- comdist(comm = arfZOTU_stein, dis = geneticDist_ZOTU_mat, abundance.weighted = FALSE)

saveRDS(arfZOTU_geneticDist_laup, file.path(data.dir, "arfZOTU_geneticDist_laup.rds"))
saveRDS(arfZOTU_geneticDist_stein, file.path(data.dir, "arfZOTU_geneticDist_stein.rds"))

geneticDist_OTU_mat <- readRDS(file.path(data.dir, "geneticDist_OTU_mat.rds"))
arfOTU_geneticDist_laup <- comdist(comm = arfOTU_laup, dis = geneticDist_OTU_mat, abundance.weighted = FALSE)
arfOTU_geneticDist_stein <- comdist(comm = arfOTU_stein, dis = geneticDist_OTU_mat, abundance.weighted = FALSE)

saveRDS(arfOTU_geneticDist_laup, file.path(data.dir, "arfOTU_geneticDist_laup.rds"))
saveRDS(arfOTU_geneticDist_stein, file.path(data.dir, "arfOTU_geneticDist_stein.rds"))

# Subset data for testing
write.csv(arfZOTUIDdata[[1]], file.path(data.dir, "arf_ZOTU_subset.csv"))
write.csv(arfOTUIDdata[[1]], file.path(data.dir, "arf_OTU_subset.csv"))

arfZOTU_laup <- arfZOTUIDdata[[1]][rownames(arfZOTUIDdata[[1]]) %in% laupahoehoe_siteIDs,]
arfZOTU_stein <- arfZOTUIDdata[[1]][rownames(arfZOTUIDdata[[1]]) %in% steinbeck_siteIDs,]

arfOTU_laup <- arfOTUIDdata[[1]][rownames(arfOTUIDdata[[1]]) %in% laupahoehoe_siteIDs,]
arfOTU_stein <- arfOTUIDdata[[1]][rownames(arfOTUIDdata[[1]]) %in% steinbeck_siteIDs,]

saveRDS(arfZOTU_laup, file.path(data.dir, "arfZOTU_laup.rds"))
saveRDS(arfZOTU_stein, file.path(data.dir, "arfZOTU_stein.rds"))

saveRDS(arfOTU_laup, file.path(data.dir, "arfOTU_laup.rds"))
saveRDS(arfOTU_stein, file.path(data.dir, "arfOTU_stein.rds"))