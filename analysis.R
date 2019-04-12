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
res.dir <- file.path(analysis.dir, "results")
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

## Calculate hypervolumes ==== 
library(hypervolume)
# Determine suitability (i.e., sample size) for hypervolume analysis
determineSuitability <- function(x){
  suitable <- ifelse(sum(x$rarefiedReadAbund > 0) >= 5, "yes", "no")
  return(data.frame(suitable))
}

OTU_df_suit <- ddply(.data = OTU_df3, .fun = determineSuitability, .variables = .(site,OTU_ID))
OTU_df_suit2 <- merge(OTU_df3, OTU_df_suit, by = c("site", "OTU_ID"))

# Exclude all unsuitable taxa
target.col <- c("rf_ann", "t_ann")
hypervolume_gaussian(data = test[target.col], weight = test$rarefiedReadAbund)

OTU_hypervols <- dlply(.data = subset(OTU_df_suit2, suitable == "yes" & rarefiedReadAbund > 0 ),
                       .variables = .(site, OTU_ID),
                       .fun = function(x){ list(hypervol = 
                                                  hypervolume_gaussian( x[c("rf_ann", "t_ann")],
                                                                 weight = x$rarefiedReadAbund),
                                                site = x$site[1],
                                                OTU_ID = x$OTU_ID[1],
                                                Order = x$V4[1]) } )
saveRDS(OTU_hypervols, file.path(res.dir, "OTU_hypervol.rds"))

library(plyr)
OTU_hypervols <- readRDS(file.path(res.dir, "OTU_hypervol.rds"))
OTU_hypervol_df <- ldply(.data = OTU_hypervols, .fun = function(x){data.frame( hypervol = x$hypervol@Volume, site = x$site, OTU_ID = x$OTU_ID, Order = x$Order) } )

OTU_hypervol_df2 <- dcast(OTU_hypervol_df, formula = OTU_ID ~ site, value.var = "hypervol")
ggplot(OTU_hypervol_df2) + geom_point(aes(y = Laupahoehoe, x = Stainback)) + geom_abline(slope = 1, intercept = 0)


ggplot(data = OTU_hypervol_df) + geom_boxplot(aes(y = hypervol, x = site, fill = Order))
ggplot(data = OTU_hypervol_df) + geom_boxplot(aes(y = hypervol, x = site, fill = Order))

# Calculate quantiles ======
with(OTU_df3, tapply(rf_ann, INDEX = site, FUN = range))
with(OTU_df3, tapply(t_ann, INDEX = site, FUN = range))

findLimits <- function(x){
  rf_v <- rep(x$rf_ann, x$rarefiedReadAbund)
  temp_v <- rep(x$t_ann, x$rarefiedReadAbund)
  rf_q <- quantile(rf_v, probs = c(0.025, 0.5, 0.975))
  temp_q <- quantile(temp_v, probs = c(0.025, 0.5, 0.975))
  temp_lower <- temp_q[1] 
  temp_median <- temp_q[2]
  temp_upper <- temp_q[3]
  rf_lower <- rf_q[1]
  rf_median <- rf_q[2]
  rf_upper <- rf_q[3]
  data.frame(temp_lower, temp_upper, temp_median, rf_lower, rf_upper, rf_median, order = x$V4[1])
}
test <- ddply(.data = subset(OTU_df_suit2, suitable == "yes"), .variables = .(site, OTU_ID), .fun = findLimits)


subset(OTU_df_suit2, site == "Laupahoehoe" & OTU_ID == "OTU300" & rarefiedReadAbund >0)
subset(OTU_df_suit2, site == "Laupahoehoe" & OTU_ID == "OTU74" & rarefiedReadAbund >0)

test$rf_range <- test$rf_upper - test$rf_lower
test$rf_site_upperlimit <- ifelse(test$site == "Laupahoehoe", ifelse(test$rf_upper == 4402.763, 1, 0), ifelse(test$rf_upper == 6388.352, 1, 0))
test$rf_site_lowerlimit <- ifelse(test$site == "Laupahoehoe", ifelse(test$rf_lower == 2840.776, 1, 0), ifelse(test$rf_lower == 2480.524, 1, 0))

test$temp_range <- test$temp_upper - test$temp_lower
test$temp_site_upperlimit <- ifelse(test$site == "Laupahoehoe", ifelse(test$temp_upper == 17.89003, 1, 0), ifelse(test$temp_upper == 18.16707, 1, 0))
test$temp_site_lowerlimit <- ifelse(test$site == "Laupahoehoe", ifelse(test$temp_lower == 13.54576, 1, 0), ifelse(test$temp_lower == 12.28192, 1, 0))

ggplot(data = subset(test, !(temp_site_lowerlimit==1 | temp_site_upperlimit == 1))) + 
  geom_boxplot(aes(y = temp_range, x = site, fill = order))
library(lmerTest)
mod <- lmer(temp_range ~ site + (1|order), 
     data = subset(test, !(temp_site_lowerlimit==1 | temp_site_upperlimit == 1)))
MuMIn::r.squaredGLMM(mod)

summary(lm(temp_range ~ site, 
   data = subset(test, !(temp_site_lowerlimit==1 | temp_site_upperlimit == 1))))

y <- dcast( subset(test, !(temp_site_lowerlimit==1 | temp_site_upperlimit == 1)), formula = OTU_ID ~ site, value.var = "temp_median")
taxonDataOTU <- taxonData[c("OTU_ID", "V4")][!duplicated(taxonData[c("OTU_ID", "V4")]),]
y2 <- merge(x = y, y = taxonDataOTU, by = "OTU_ID", all = FALSE)
ggplot(aes(y = Laupahoehoe, x = Stainback),data = y2) + geom_point(aes(color = V4)) + geom_abline(slope = 1, intercept = 0) + geom_smooth(method = "lm")


# The weights do not deal with zero weights as expected...
test <- data.frame(y = rnorm(100, 0, 1), x = rnorm(100,0,1))
plot(y~ x, data = test)
x <- hypervolume_gaussian(data = test, weight = rep(1,100))
plot(x)
y <- hypervolume_gaussian(data = rbind(test, c(500, 500)), weight = c(rep(1,100), 0) )



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
