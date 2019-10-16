## ECOLOGICAL ASSEMBLY OF THE HAWAIIAN ARTHROPOD COMMUNITIES
# Author: Jun Ying Lim
# Rarefies OTU tables for analysis

## Packages ============
rm(list = ls())
library(stringr); library(plyr); library(reshape2) # data manipulation tools
library(vegan) # calculating beta diversity
library(geosphere) # calculating geographic distances
library(ggplot2); library(ggrepel)

## Import data ============
main.dir <- "~/Dropbox/Projects/2017/hawaiiCommunityAssembly/elevationAssemblyHawaii"
data.dir <- file.path(main.dir, "data")
fig.dir <- file.path(main.dir, "figures")
res.dir <- file.path(main.dir, "results")
source(file.path(main.dir, "metabarcodingTools.R"))
source(file.path(main.dir, "test.R"))

# Import site data
siteData <- read.csv(file.path(data.dir, "clim.final.csv"), stringsAsFactors = FALSE)
siteData <- siteData[-grep(siteData$site.id, pattern = "BRG"),] # Exclude Rosie's samples

laupahoehoe_siteIDs <- subset(siteData, site.1 == "Laupahoehoe")$site.id
stainback_siteIDs <- subset(siteData, site.1 == "Stainback")$site.id

# Import otu data
zOTU_mat <- readRDS(file.path(data.dir, "zOTU_native.rds"))
OTU_mat <- readRDS(file.path(data.dir, "OTU_native.rds"))
SPP_mat <- readRDS(file.path(data.dir, "SPP_native.rds"))

# Import taxonmic reference
taxonData <- readRDS(file.path(data.dir, "taxonData.rds"))
OTUtaxonData <- taxonData[c("OTU_ID", "V4")]
OTUtaxonData <- OTUtaxonData[!duplicated(OTUtaxonData),]

## Abundance patterns =========
# Only analyze OTUs that have been matched to certain orders
targetOrders <- c("Acari", "Araneae", "Coleoptera", "Lepidoptera", "Orthoptera", "Neuroptera", "Psocoptera", "Hemiptera")
otu_target <- which(taxonData$V4[match(rownames(OTU_mat), taxonData$OTU_ID)] %in% targetOrders)
OTU_df <- melt(data = OTU_mat[otu_target,], value.name = "nReads", varnames = c("OTU_ID", "site.id"))
OTU_clim_df <- merge(x = OTU_df, y = siteData[c("site","site.id", "rf_ann", "t_ann")], by = "site.id", all.x = TRUE)

# Count number of unique OTUs per site
site_nsp <- ddply(.data = subset(OTU_clim_df, nReads > 0), .variables = .(site.id), .fun = function(x) {data.frame(nsp = length(unique(x$OTU_ID)))}) 
write.csv(site_nsp, file.path(res.dir, "site_nsp.csv"), row.names = F)

# Remove OTUs with fewer than 10 sites
removeLowOcc <- function(data){
  if(sum(data$nReads > 0) >= 10){
    return(data)
  }
}
OTU_clim_df2 <- ddply(.data = OTU_clim_df, .variables = .(site, OTU_ID), .fun = removeLowOcc)

# Calculate niche and niche width for each OTU for each site
OTU_clim_site <- ddply(.data = OTU_clim_df2, .variables = .(site, OTU_ID), .fun = findLimits)

# Flag OTUs whose range boundaries are outside the sampling ranges for each site ( 1 = in bounds, 0 = out of bounds; since upper and lower bounds are 0.75 and 0.25; not more than a quarter of occurrences may be at the most extreme site)
OTU_clim_site$rf_bounds <- with(data = OTU_clim_site,
                                ifelse(site == "Laupahoehoe" & rf_upper >= 4402.763 |
                                       site == "Laupahoehoe" & rf_lower <= 2840.776 | 
                                       site == "Stainback" & rf_upper >= 6388.352 |
                                       site == "Stainback" & rf_lower <= 2480.524, 0, 1))

OTU_clim_site$temp_bounds <- with(data = OTU_clim_site, 
                               ifelse(site == "Laupahoehoe" & temp_upper >= 17.89003 |
                                      site == "Laupahoehoe" &  temp_lower <= 13.54576 |
                                      site == "Stainback" & temp_upper >= 18.16707 |
                                      site == "Stainback" & temp_lower <= 12.28192, 0, 1))

table(OTU_clim_site$rf_bounds)
table(OTU_clim_site$temp_bounds)

# Calculate niche breadth
OTU_clim_site$temp_breadth <- OTU_clim_site$temp_upper - OTU_clim_site$temp_lower
OTU_clim_site$rf_breadth <- OTU_clim_site$rf_upper - OTU_clim_site$rf_lower

# Only include species in both sites
bothSites <- function(x){
  if(sum(c("Laupahoehoe", "Stainback") %in% x$site) == 2){
    return(x)
  } else {
    return(NULL)
  }
}

t_res <- ddply(.data =  subset(OTU_clim_site, temp_bounds == 1), .variables = .(OTU_ID), .fun = bothSites)
rf_res <- ddply(.data = subset(OTU_clim_site, rf_bounds == 1), .variables = .(OTU_ID), .fun = bothSites)

OTU_taxonref <- taxonData[c("OTU_ID","V4")]
OTU_taxonref <- OTU_taxonref[!duplicated(OTU_taxonref),]


## Question 1: Are niches conserved across sites?
t_med_site <- dcast(data = t_res, formula = OTU_ID~site, value.var = "temp_median")
rf_med_site <- dcast(data = rf_res, formula = OTU_ID~site, value.var = "rf_median")

t_med_site <- merge(t_res_site, OTU_taxonref)
rf_med_site <- merge(rf_res_site, OTU_taxonref)

write.csv(t_med_site, file.path(res.dir, "t_median_site.csv"), row.names = FALSE)
write.csv(rf_med_site, file.path(res.dir, "rf_median_site.csv"), row.names = FALSE)

cor.test(x = t_med_site$Laupahoehoe, y = t_med_site$Stainback)
cor.test(x = rf_med_site$Laupahoehoe, y = rf_med_site$Stainback)

## Question 2: Is niche breadth correlated across sites?
t_breadth_site <- dcast(data = t_res, formula = OTU_ID~site, value.var = "temp_breadth")
rf_breadth_site <- dcast(data = rf_res, formula = OTU_ID~site, value.var = "rf_breadth")

t_breadth_site <- merge(t_breadth_site, OTU_taxonref)
rf_breadth_site <- merge(rf_breadth_site, OTU_taxonref)

cor.test(x = t_breadth_site$Laupahoehoe, y = t_breadth_site$Stainback)
cor.test(x = rf_breadth_site$Laupahoehoe, y = rf_breadth_site$Stainback)

write.csv(t_breadth_site, file.path(res.dir, "t_breadth_site.csv"), row.names = FALSE)
write.csv(rf_breadth_site, file.path(res.dir, "rf_breadth_site.csv"), row.names = FALSE)

## Question 3: Is niche breadth higher in one site
mod <- lm(temp_breadth ~ site, data = merge(t_res, OTU_taxonref))
summary(mod)

library(lmerTest)
library(lme4)
mod <- lmer(temp_breadth ~ site + (1|V4), data = merge(t_res, OTU_taxonref))
summary(mod)
mod2 <- lmer(temp_breadth ~ site + (1|V4), data = merge(t_res, OTU_taxonref))
summary(mod2)
anova(mod2,mod)
# Higher in Stainback

## Question 4: Is niche breadth higher across groups
mod <- lm(temp_breadth ~ V4, data = subset(merge(t_res, OTU_taxonref), site == "Laupahoehoe"))
summary(mod)
mod2 <- lm(temp_breadth ~ V4, data = subset(merge(t_res, OTU_taxonref), site == "Stainback"))
summary(mod2)

mod <- lm(rf_breadth ~ V4, data = subset(merge(rf_res, OTU_taxonref), site == "Laupahoehoe"))
summary(mod)
mod2 <- lm(rf_breadth ~ V4, data = subset(merge(rf_res, OTU_taxonref), site == "Stainback"))
summary(mod2)

