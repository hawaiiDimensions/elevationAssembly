## ECOLOGICAL ASSEMBLY OF THE HAWAIIAN ARTHROPOD COMMUNITIES
# Author: Jun Ying Lim
# Rarefies OTU tables for analysis

## PACKAGES ============
library(stringr); library(plyr); library(reshape2) # data manipulation tools
library(vegan) # calculating beta diversity
library(geosphere) # calculating geographic distances
library(ggplot2); library(ggrepel)

## IMPORT DATA ============
main.dir <- "~/Dropbox/Projects/2017/hawaiiCommunityAssembly/"
analysis.dir <- file.path(main.dir, "elevationAssemblyHawaii")
data.dir <- file.path(main.dir, "data")
fig.dir <- file.path(analysis.dir, "figures")
source(file.path(analysis.dir, "metabarcodingTools.R"))

# Import site data
siteData <- read.csv(file.path(data.dir, "clim.final.csv"), stringsAsFactors = FALSE)
siteData <- siteData[-grep(siteData$site.id, pattern = "BRG"),] # Exclude Rosie's samples

laupahoehoe_siteIDs <- subset(siteData, site.1 == "Laupahoehoe")$site.id
steinbeck_siteIDs <- subset(siteData, site.1 == "Stainback")$site.id

# Import otu data
zOTU_mat <- readRDS(file.path(data.dir, "zOTU_native.rds"))
OTU_mat <- readRDS(file.path(data.dir, "OTU_native.rds"))
SPP_mat <- readRDS(file.path(data.dir, "SPP_native.rds"))

# Import taxonmic reference
taxonData <- readRDS(file.path(data.dir, "taxonData.rds"))

## ABUNDANCE PATTERNS =========
aranaeae_ZOTU <- subset(taxonData, V4 == "Araneae")$zOTU_ID

OTU_df <- melt(OTU_mat, value.name = "rarefiedReadAbund", varnames = c("OTU_ID", "site.id"))
OTU_df2 <- merge(x = OTU_df, y = taxonData[c("OTU_ID", "V4")], by = "OTU_ID", all.x = TRUE)

# Remove duplicates
OTU_df2 <- OTU_df2[!duplicated(OTU_df2),]
OTU_df3 <- merge(x = OTU_df2, y = siteData[c("site.id", "rf_ann", "t_ann", "site")], by = "site.id")

library(hypervolume)
# Determine suitability (i.e., sample size) for hypervolume analysis
determineSuitability <- function(x){
  # Determine suitability for hypervolume analysis
  suitable <- ifelse(sum(x$rarefiedReadAbund > 0) >= 5, "yes", "no")
  return(data.frame(suitable))
}

OTU_df_suit <- ddply(.data = OTU_df3, .fun = determineSuitability, .variables = .(site,OTU_ID))
OTU_df_suit2 <- merge(OTU_df3, OTU_df_suit, by = c("site", "OTU_ID"))

# Exclude all unsuitable taxa
test <- subset(OTU_df_suit2, OTU_ID == "OTU100", site = "Laupahoehoe")

target.col <- c("rf_ann", "t_ann")
hypervolume_gaussian(data = test[target.col], weight = test$rarefiedReadAbund)

OTU_hypervols <- dlply(.data = subset(OTU_df_suit2, suitable == "yes"),
                       .variables = .(site, OTU_ID),
                       .fun = function(x){ list(hypervol = 
                                                  hypervolume_gaussian( x[c("rf_ann", "t_ann")],
                                                                 weight = x$rarefiedReadAbund),
                                                site = x$site,
                                                OTU_ID = x$OTU_ID,
                                                Order = x$V4) } )

plot(OTU_hypervols[[3]])


hypervolume_gaussian(data = OTU_df3[target.col], weight )


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
