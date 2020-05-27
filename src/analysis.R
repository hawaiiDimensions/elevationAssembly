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
main.dir <- "~/Dropbox/projects/2017/hawaiiCommunityAssembly/elevationAssemblyHawaii/"
src.dir <- file.path(main.dir, "src")
data.dir <- file.path(main.dir, "data")
fig.dir <- file.path(main.dir, "figures")
res.dir <- file.path(main.dir, "results")
source(file.path(src.dir, "metabarcodingTools.R"))
source(file.path(src.dir, "test.R"))

# Import site data
siteData <- read.csv(file.path(data.dir, "clim.final.csv"), stringsAsFactors = FALSE)
siteData <- siteData[-grep(siteData$site.id, pattern = "BRG"),] # Exclude Rosie's samples

subset(siteData, site == "Stainback" & elevation < 800)

laupahoehoe_siteIDs <- subset(siteData, site.1 == "Laupahoehoe")$site.id
stainback_siteIDs <- subset(siteData, site.1 == "Stainback")$site.id

# Import otu data
otu_native_files <- list.files(file.path(data.dir, "OTU_native"))

OTU_mat_list <- lapply(otu_native_files, FUN = function(x){readRDS(file.path(data.dir, "OTU_native", x))})
#OTU_mat <- readRDS(file.path(data.dir, "OTU_native.rds"))
#OTU_mat <- OTU_mat_list[[1]]

# Import taxonmic reference
taxonData <- readRDS(file.path(data.dir, "taxonData.rds"))
OTUtaxonData <- taxonData[c("OTU_ID", "V4")]
OTUtaxonData <- OTUtaxonData[!duplicated(OTUtaxonData),]

## Calculate alpha diversity with elevation  =========
OTU_mat_df <- melt(OTU_mat_list, varnames = c("OTU_ID", "site_id"),
                   value.name = "nReads")

OTU_mat_nsp <- ddply(subset(OTU_mat_df, nReads > 0), .variables = .(L1, site_id),
                     .fun = function(x) {data.frame(nsp = length(as.vector(unique(x$OTU_ID))))},
                     .progress = "text")

OTU_mat_nsp_site <- merge(x = OTU_mat_nsp, y = siteData, by.x = "site_id", by.y= "site.id", all.x = T)
saveRDS(OTU_mat_nsp_site, file.path(res.dir, "OTU_nsp_site.rds"))

# # Look at pattern of unique OTUs vs. shared OTUs
# asd <- merge(x= subset(OTU_mat_df, nReads > 0), y = siteData[,c("site","site.id")], by.x = "site_id", by.y = "site.id")
# 
# tes <- ddply(asd, .variables = .(L1, OTU_ID), .fun = summarize,
#              n_site = length(unique(as.vector(site))))
# 
# tes2 <- merge(x = asd, y = tes, by.x = c("L1", "OTU_ID"), by.y = c("L1", "OTU_ID") )
# 
# tes3 <- ddply(.data = tes2, .variables = .(L1, site_id), .fun = summarize,
#       nSp = length(unique(as.vector(OTU_ID))),
#       nUniSp = length(unique(as.vector(OTU_ID[n_site == 1]))),
#       propEnd = (nUniSp / nSp)*100 )
# 
# tes4 <- merge(tes3, siteData[,c("site", "site.id", "elevation")], by.x = "site_id", by.y = "site.id")
# 
# tes5 <- merge(tes4, taxonData)
# head(tes5)
# ggplot(data = tes4) +
#   geom_point(aes(y = propEnd, x = elevation, colour = site, group = L1)) + 
#   geom_smooth(aes(y = propEnd, x= elevation, group = site),se = F,method = "lm")

## Calculate the average distance between sampling localities between sites
# library(geosphere)
# test <- as.matrix(siteData[c("longitude", "latitude")])
# rownames(test) <- siteData$site.id
# 
# pw_geogdist <- distm(test, test, fun = distGeo)
# mean(pw_geogdist[which(siteData$site == "Laupahoehoe"), which(siteData$site == "Stainback")])

## NMDS of beta diversity across elevation =========
# Identify sites lower than 800 m elevation
lowsiteids <- subset(siteData, elevation <= 800)$site.id

# Calculate the average rarefied read abundance for all OTUs
OTU_mat_avg_df <- ddply(OTU_mat_df, .variables = .(site_id, OTU_ID), summarize, avgRead = mean(nReads), .progress = "text")
OTU_mat_avg <- reshape2::acast(OTU_mat_avg_df, formula = site_id~OTU_ID, value.var = "avgRead")

OTU_mat_avg_subset <- OTU_mat_avg[!rownames(OTU_mat_avg) %in% (lowsiteids),]

# Run NMDS
site_nmds<- metaMDS(OTU_mat_avg_subset, k = 2)

site_nmds_df <- data.frame(site_nmds$points)
site_nmds_df$site.id <- rownames(site_nmds$points)
site_nmds_df2 <- merge(site_nmds_df, siteData)

species_nmds_df <- data.frame(site_nmds$species)
species_nmds_df$species <- rownames(site_nmds$species)

site_nmds_contour <- ordisurf(site_nmds, site_nmds_df2$elevation)
extract.xyz <- function(obj) {
  xy <- expand.grid(x = obj$grid$x, y = obj$grid$y)
  xyz <- cbind(xy, c(obj$grid$z))
  names(xyz) <- c("x", "y", "z")
  return(xyz)
}
site_nmds_contour2 <- extract.xyz(site_nmds_contour)

write.csv(site_nmds_df2, file.path(res.dir, "site_nmds.csv"), row.names = F)
write.csv(site_nmds_contour2, file.path(res.dir, "site_nmds_contour.csv"), row.names = F)

# Beta-diversity over elevation
x <- as.matrix(vegdist(OTU_mat_avg_subset, method = "bray", diag = T, upper = T))
xy <- t(combn(colnames(x), 2))
distDF <- data.frame(xy, dist= x[xy])

distDF2 <- merge(x = distDF, y = siteData[, c("site.id", "elevation", "site")], by.x = "X1", by.y = "site.id")
distDF3 <- merge(x = distDF2, y = siteData[, c("site.id", "elevation", "site")], by.x = "X2", by.y = "site.id", suffixes = c("_X1", "_X2"))

# Bin sites into elevational categories
elev_bins <- c("500", "800", "1100", "1400", "1700", "2000")
elev_bin_labels <- c("< 800", "800 - 1,100", "1,100 - 1,400", "1,400 - 1,700", "> 1,700")

distDF3$elevation_X1_class <- cut(distDF3$elevation_X1, breaks = elev_bins, labels = elev_bin_labels)
distDF3$elevation_X2_class <- cut(distDF3$elevation_X2, breaks = elev_bins, labels = elev_bin_labels)
write.csv(distDF3, file = file.path(res.dir, "elevation_dist.csv"), row.names = F)

# Only analyze between-transect comparisons
site_elevation_beta <- subset(distDF3, site_X1 != site_X2 & elevation_X1_class == elevation_X2_class)
table(site_elevation_beta$elevation_X2_class) # sampling is fairly uneven

site_elevation_beta_anova <- aov(dist ~ elevation_X1_class, site_elevation_beta)
summary(site_elevation_beta_anova)
site_elevation_beta_tukey <- TukeyHSD(x = site_elevation_beta_anova)
write.csv(site_elevation_beta_tukey$elevation_X1_class, file = file.path(res.dir, "site_elevation_beta_tukey.csv"))

## Beta-diversity within OTUs =========================
# Import zotu tables
zotu_native_files <- list.files(file.path(data.dir, "zOTU_native"))
zOTU_mat_list <- lapply(zotu_native_files, FUN = function(x){readRDS(file.path(data.dir, "zOTU_native", x))})

# Collapse rarefied zotu tables into a single dataframe
zOTU_mat_df <- melt(zOTU_mat_list, varnames = c("zOTU_ID", "site_id"),
                     value.name = "nReads")
zOTU_mat_df_subset <- subset(zOTU_mat_df, !site_id %in% lowsiteids)

# Take average no. of reads for each zotu tables across rarefied tables
zOTU_mat_avg_df <- ddply(zOTU_mat_df_subset,
                         .variables = .(site_id, zOTU_ID),
                         summarize, avgRead = mean(nReads), .progress = "text")
zOTU_taxon_df <- merge(zOTU_mat_avg_df, taxonData, by = "zOTU_ID")

# Only include OTUs with more than 1 OTU
zOTU_taxon_df_subset <- ddply(subset(zOTU_taxon_df, avgRead > 0),
      .variables = .(OTU_ID), .fun = function(x){ if(length(unique(x$zOTU_ID)) > 1){return(x)} else {return(NULL)}  })

# Calculate beta-diversity of zOTUs for each OTU
zOTU_beta <- ddply(.data = zOTU_taxon_df_subset,
      .variables = .(OTU_ID),
      .fun = function(x){
        #x = subset(zOTU_taxon_df, OTU_ID == "OTU31")
        OTU_mat <- reshape2::acast(x, formula = site_id~zOTU_ID, value.var = "avgRead", fill = 0)
        OTU_mat2 <- OTU_mat[rowSums(OTU_mat) > 0,]
        beta_mat <- as.matrix(vegdist(OTU_mat2, method = "bray", diag = T, upper = T))
        beta_dist <- t(combn(colnames(beta_mat), 2))
        distDF <- data.frame(beta_dist, dist= beta_mat[beta_dist])
        return(distDF)
      })

# Calculate average distance across OTUs
zOTU_beta_avg <- ddply(zOTU_beta, .variables = .(X2,X1),
      .fun = function(x){
        nSharedOTUs = length(x$OTU_ID) # number of OTUs compared between sites
        avgDist = mean(x$dist)
        data.frame(nSharedOTUs, avgDist)
      })

# Merge with site data
zOTU_beta_avg2 <- merge(x = zOTU_beta_avg, y = siteData[, c("site.id", "elevation", "site")], by.x = "X1", by.y = "site.id")
zOTU_beta_avg3 <- merge(x = zOTU_beta_avg2, y = siteData[, c("site.id", "elevation", "site")], by.x = "X2", by.y = "site.id", suffixes = c("_X1", "_X2"))

# Group into elevational bins
elev_bins <- c("500", "800", "1100", "1400", "1500", "2000")
elev_bin_labels <- c("< 800", "800 - 1,100", "1,100 - 1,400", "1,400 - 1,500", "> 1,500")

zOTU_beta_avg3$elevation_X1_class <- cut(zOTU_beta_avg3$elevation_X1, breaks = elev_bins, labels = elev_bin_labels)
zOTU_beta_avg3$elevation_X2_class <- cut(zOTU_beta_avg3$elevation_X2, breaks = elev_bins, labels = elev_bin_labels)

write.csv(zOTU_beta_avg3, file = file.path(res.dir, "elevation_ZOTU_dist.csv"), row.names = F)

# length(subset(siteData, elevation < 1600 & elevation> 800 & site == "Laupahoehoe")$site.id)
# length(subset(siteData, elevation < 1600 & elevation> 800 & site == "Stainback")$site.id)

# Analyze cross-site patterns
zOTU_beta_acrosssite <- subset(zOTU_beta_avg3, elevation_X1_class == elevation_X2_class & site_X1 != site_X2)
zOTU_beta_anova <- aov(avgDist ~ elevation_X1_class, data = zOTU_beta_acrosssite)
zOTU_beta_tukey <- TukeyHSD(zOTU_beta_anova)
write.csv(zOTU_beta_tukey$elevation_X1_class, file = file.path(res.dir, "site_elevation_zotu_beta_tukey.csv"))

ggplot(data = zOTU_beta_acrosssite) +
  geom_violin(aes(y = avgDist, x = elevation_X1_class, fill = elevation_X1_class)) +
  geom_point(aes(y = avgDist, x = elevation_X1_class)) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_blank()) +
  labs(y = "Pairwise Bray-Curtis distance", x = "Elevation bin (m)") +
  scale_fill_manual(values = wesanderson::wes_palette("IsleofDogs2", n = 4))

# Analyze within-site patterns (relative to elevation_X1_class)
zOTU_beta_laupa <- subset(zOTU_beta_avg3, elevation_X1_class == "800 - 1,100" &
                            elevation_X2_class %in% c("800 - 1,100", "1,100 - 1,400", "1,400 - 1,500", "> 1,500") &
                            site_X1 == "Laupahoehoe" &
                            site_X2 == "Laupahoehoe")
# SHOULD IT BE ELEV1 OR ELEV2 == "800 - 1,100"?
ggplot(data = zOTU_beta_laupa) +
  geom_violin(aes(y = avgDist, x = elevation_X2_class, fill = elevation_X2_class)) +
  geom_point(aes(y = avgDist, x = elevation_X2_class))

summary(aov(avgDist ~ elevation_X1_class, data = zOTU_beta_laupa))

zOTU_beta_stain <- subset(zOTU_beta_avg3, elevation_X1_class =="800 - 1,100" &
                            elevation_X2_class %in% c("800 - 1,100", "1,100 - 1,400", "1,400 - 1,700", "> 1,700") &
                            site_X1 == "Stainback" &
                            site_X2 == "Stainback")
zOTU_beta_stain$elevation_X2_class <- factor(as.vector(zOTU_beta_stain$elevation_X2_class), levels = elev_bin_labels)
ggplot(data = zOTU_beta_stain) +
  geom_violin(aes(y = avgDist, x = elevation_X2_class, fill = elevation_X2_class)) +
  geom_point(aes(y = avgDist, x = elevation_X2_class))

summary(aov(avgDist ~ elevation_X1_class, data = zOTU_beta_stain))



test <- subset(zOTU_beta_avg3, elevation_X1_class == elevation_X2_class & site_X1 != site_X2)
zotu_beta_anova<- aov(avgDist ~ elevation_X1_class, data = test)
zotu_beta_tukey <- TukeyHSD(zotu_beta_anova)

## ABUNDANCE PATTERNS =========
t_res_list_t10 <- list()
rf_res_list_t10 <- list()
t_res_list_t5 <- list()
rf_res_list_t5 <- list()

for(i in 1:length(OTU_mat_list)){
  # Counter
  pb = txtProgressBar(min = 0, max = length(OTU_mat_list), style = 3, initial = 0) 
  setTxtProgressBar(pb, i)

  OTU_mat <- OTU_mat_list[[i]]
  
  # Only analyze OTUs that have been matched to certain orders
  targetOrders <- c("Acari", "Araneae", "Coleoptera", "Lepidoptera", "Orthoptera", "Neuroptera", "Psocoptera", "Hemiptera")
  #<- c("Araneae", "Hemiptera", "Lepidoptera", "Psocoptera", "Coleoptera", "Orthoptera")
  otu_target <- which(taxonData$V4[match(rownames(OTU_mat), taxonData$OTU_ID)] %in% targetOrders)
  
  # Convert OTU matrix into long table format
  OTU_df <- melt(data = OTU_mat[otu_target,], value.name = "nReads", varnames = c("OTU_ID", "site.id"))
  
  # Combine OTU table with site climatic data
  OTU_clim_df <- merge(x = OTU_df, y = siteData[c("site","site.id", "rf_ann", "t_ann")], by = "site.id", all.x = TRUE)
  
  # Remove OTUs with fewer than 10 sites on each gradient
  removeLowOcc <- function(data, threshold){
    if(sum(data$nReads > 0) >= threshold){
      return(data)
    }
  }
  OTU_clim_df2_10 <- ddply(.data = OTU_clim_df, .variables = .(site, OTU_ID), .fun = removeLowOcc, threshold = 10)  
  OTU_clim_df2_5 <- ddply(.data = OTU_clim_df, .variables = .(site, OTU_ID), .fun = removeLowOcc, threshold = 5) 
  
  # Calculate niche and niche width for each OTU for each site
  calcNiche <- function(x){
    rf_v <- rep(x$rf_ann, x$nReads)
    temp_v <- rep(x$t_ann, x$nReads)
    rf_q <- quantile(rf_v, probs = c(0.25, 0.5, 0.75))
    temp_q <- quantile(temp_v, probs = c(0.25, 0.5, 0.75))
     
    # Niche position
    temp_median <- temp_q[2]
    rf_median <- rf_q[2]
    
    # Niche breadth
    temp_upper <- temp_q[3]
    temp_lower <- temp_q[1]
    temp_breadth <- temp_upper - temp_lower
    
    rf_upper <- rf_q[3]
    rf_lower <- rf_q[1]
    rf_breadth <- rf_upper - rf_lower
    
    data.frame(temp_lower, temp_upper, temp_median, temp_breadth, rf_lower, rf_upper, rf_median, rf_breadth)
  }
  OTU_clim_site_10 <- ddply(.data = OTU_clim_df2_10, .variables = .(site, OTU_ID), .fun = calcNiche)
  OTU_clim_site_5 <- ddply(.data = OTU_clim_df2_5, .variables = .(site, OTU_ID), .fun = calcNiche)
  
  # Flag OTUs whose range boundaries are outside the sampling ranges for each site ( 1 = in bounds, 0 = out of bounds; since upper and lower bounds are 0.75 and 0.25; not more than a quarter of occurrences may be at the most extreme site)
  flagLimit <- function(x){
    x$rf_bounds <- with(data = x,
                        ifelse(site == "Laupahoehoe" & rf_upper >= 4402.763 |
                            site == "Laupahoehoe" & rf_lower <= 2840.776 | 
                            site == "Stainback" & rf_upper >= 6388.352 |
                            site == "Stainback" & rf_lower <= 2480.524, 0, 1))
    x$temp_bounds <- with(data = x,
                          ifelse(site == "Laupahoehoe" & temp_upper >= 17.89003 |
                                   site == "Laupahoehoe" &  temp_lower <= 13.54576 |
                                   site == "Stainback" & temp_upper >= 18.16707 |
                                   site == "Stainback" & temp_lower <= 12.28192, 0, 1))
    return(x)
  }
  OTU_clim_site_10 <- flagLimit(x = OTU_clim_site_10)
  OTU_clim_site_5 <- flagLimit(x = OTU_clim_site_5)
  
  # table(OTU_clim_site_10$rf_bounds)
  # table(OTU_clim_site_10$temp_bounds)
  # table(OTU_clim_site_5$rf_bounds)
  # table(OTU_clim_site_5$temp_bounds)
  
  t_res_list_t10[[i]] <- subset(OTU_clim_site_10, temp_bounds == 1)
  rf_res_list_t10[[i]] <- subset(OTU_clim_site_10, rf_bounds == 1)
  t_res_list_t5[[i]] <- subset(OTU_clim_site_5, temp_bounds == 1)
  rf_res_list_t5[[i]] <- subset(OTU_clim_site_5, rf_bounds == 1)
  
}
close(pb)

saveRDS(t_res_list_t10, file.path(res.dir, "t_res_list_t10.rds"))
saveRDS(rf_res_list_t10, file.path(res.dir, "rf_res_list_t10.rds"))
saveRDS(t_res_list_t5, file.path(res.dir, "t_res_list_t5"))
saveRDS(rf_res_list_t5, file.path(res.dir, "rf_res_list_t5"))

# Only include species in both sites
bothSites <- function(x){
  ddply(.data = x, .variables = .(OTU_ID), .fun = function(x){
    if(sum(c("Laupahoehoe", "Stainback") %in% x$site) == 2){
      return(x)
    } else {
      return(NULL)
    }  
  })
}
t_res_list_t10_subset <- llply(.data = t_res_list_t10, .fun = bothSites)
rf_res_list_t10_subset <- llply(.data = rf_res_list_t10, .fun = bothSites)

t_res_list_t5_subset <- llply(.data = t_res_list_t5, .fun = bothSites)
rf_res_list_t5_subset <- llply(.data = rf_res_list_t10, .fun = bothSites)

saveRDS(t_res_list_t10_subset, file.path(res.dir, "t_res_list_t10_subset.rds"))
saveRDS(rf_res_list_t10_subset, file.path(res.dir, "rf_res_list_t10_subset.rds"))
saveRDS(t_res_list_t5_subset, file.path(res.dir, "t_res_list_t5_subset.rds"))
saveRDS(rf_res_list_t5_subset, file.path(res.dir, "rf_res_list_t5_subset.rds"))

## QUESTION 1: ARE NICHES CONSERVED ACROSS SITES? ========
nicheConsPrep <- function(x, value, otu_tab){
  site <- dcast(data = x, formula = OTU_ID~site, value.var = value)
  site_OTU <- merge(site, otu_tab)
  return(site_OTU)
}

t_med_res_prep_t10  <- llply(.data = t_res_list_t10_subset, .fun = nicheConsPrep, value = "temp_median", otu_tab = OTUtaxonData)
rf_med_res_prep_t10  <- llply(.data = rf_res_list_t10_subset, .fun = nicheConsPrep, value = "rf_median", otu_tab = OTUtaxonData)

saveRDS(t_med_res_prep_t10, file.path(res.dir, "t_med_res_prep_t10.rds"))
saveRDS(rf_med_res_prep_t10, file.path(res.dir, "rf_med_res_prep_t10.rds"))

nicheConsTest <- function(x){
  cor_test_res <- cor.test(x = x$Laupahoehoe, y = x$Stainback)
  t <- cor_test_res$statistic
  cor <- cor_test_res$estimate
  p <- cor_test_res$p.value
  n <- nrow(x)
  data.frame(t, cor, p, n)
}

t_med_res_final_t10 <- ldply(.data = t_med_res_prep_t10, .fun = nicheConsTest)
rf_med_res_final_t10 <- ldply(.data = rf_med_res_prep_t10, .fun = nicheConsTest)

write.csv(t_med_res_final_t10, file.path(res.dir, "t_med_res_final_t10.csv"), row.names = T)
write.csv(rf_med_res_final_t10, file.path(res.dir, "rf_med_res_final_t10.csv"), row.names = T)

round(mean(t_med_res_final_t10$cor), 2)
round(range(t_med_res_final_t10$cor), 2)
sum(t_med_res_final_t10$p < 0.05) / 100
round(range(t_med_res_final_t10$n),2 )

round(mean(rf_med_res_final_t10$cor), 2)
round(range(rf_med_res_final_t10$cor), 2)
sum(rf_med_res_final_t10$p < 0.05) / 100
round(range(rf_med_res_final_t10$n), 1)

# Re-run analysis with lower site treshold
t_med_res_prep_t5  <- llply(.data = t_res_list_t5_subset, .fun = nicheConsPrep, value = "temp_median", otu_tab = OTUtaxonData)
rf_med_res_prep_t5  <- llply(.data = rf_res_list_t5_subset, .fun = nicheConsPrep, value = "rf_median", otu_tab = OTUtaxonData)

saveRDS(t_med_res_prep_t5, file.path(res.dir, "t_med_res_prep_t5.rds"))
saveRDS(rf_med_res_prep_t5, file.path(res.dir, "rf_med_res_prep_t5.rds"))

t_med_res_final_t5 <- ldply(.data = t_med_res_prep_t5, .fun = nicheConsTest)
rf_med_res_final_t5 <- ldply(.data = rf_med_res_prep_t5, .fun = nicheConsTest)

write.csv(t_med_res_final_t5, file.path(res.dir, "t_med_res_final_t5.csv"), row.names = T)
write.csv(rf_med_res_final_t5, file.path(res.dir, "rf_med_res_final_t5.csv"), row.names = T)

round(mean(t_med_res_final_t5$cor), 2)
round(range(t_med_res_final_t5$cor), 2)
sum(t_med_res_final_t5$p < 0.05) / 100
round(range(t_med_res_final_t5$n),2 )

round(mean(rf_med_res_final_t5$cor), 2)
round(range(rf_med_res_final_t5$cor), 2)
sum(rf_med_res_final_t5$p < 0.05) / 100
round(range(rf_med_res_final_t5$n), 1)

# Re-run analysis with OTU that occur in ALL rarefied datasets
t_med_res_prep_all_t10 <- do.call("rbind",t_med_res_prep_t10)
rf_med_res_prep_all_t10 <- do.call("rbind",rf_med_res_prep_t10)

t_med_res_prep_all_t5 <- do.call("rbind",t_med_res_prep_t5)
rf_med_res_prep_all_t5 <- do.call("rbind",rf_med_res_prep_t5)

t_med_res_prep_all_subset_t10 <- ddply(.data = t_med_res_prep_all_t10,
                               .variables = .(OTU_ID),
                               .fun = function(x) {if(nrow(x) == 100){return(x)} })
rf_med_res_prep_all_subset_t10 <- ddply(.data = rf_med_res_prep_all_t10,
                                .variables = .(OTU_ID),
                                .fun = function(x) {if(nrow(x) == 100){return(x)} })

t_med_res_prep_all_subset_t5 <- ddply(.data = t_med_res_prep_all_t5,
                                       .variables = .(OTU_ID),
                                       .fun = function(x) {if(nrow(x) == 100){return(x)} })
rf_med_res_prep_all_subset_t5 <- ddply(.data = rf_med_res_prep_all_t5,
                                        .variables = .(OTU_ID),
                                        .fun = function(x) {if(nrow(x) == 100){return(x)} })


t_med_res_prep_avg_t10 <- ddply(.data = t_med_res_prep_all_subset_t10,
                        .variables = .(OTU_ID),
                        .fun = summarize, 
                        Stainback = mean(Stainback),
                        Laupahoehoe = mean(Laupahoehoe),
                        V4 = V4[1])
rf_med_res_prep_avg_t10  <- ddply(.data = rf_med_res_prep_all_subset_t10,
                          .variables = .(OTU_ID),
                          .fun = summarize, 
                          Stainback = mean(Stainback),
                          Laupahoehoe = mean(Laupahoehoe),
                          V4 = V4[1])

t_med_res_prep_avg_t5 <- ddply(.data = t_med_res_prep_all_subset_t5,
                                .variables = .(OTU_ID),
                                .fun = summarize, 
                                Stainback = mean(Stainback),
                                Laupahoehoe = mean(Laupahoehoe),
                                V4 = V4[1])
rf_med_res_prep_avg_t5  <- ddply(.data = rf_med_res_prep_all_subset_t5,
                                  .variables = .(OTU_ID),
                                  .fun = summarize, 
                                  Stainback = mean(Stainback),
                                  Laupahoehoe = mean(Laupahoehoe),
                                  V4 = V4[1])

write.csv(t_med_res_prep_avg_t10, file.path(res.dir, "t_med_res_prep_avg_t10.csv"), row.names = T)
write.csv(rf_med_res_prep_avg_t10, file.path(res.dir, "rf_med_res_prep_avg_t10.csv"), row.names = T)

write.csv(t_med_res_prep_avg_t5, file.path(res.dir, "t_med_res_prep_avg_t5.csv"), row.names = T)
write.csv(rf_med_res_prep_avg_t5, file.path(res.dir, "rf_med_res_prep_avg_t5.csv"), row.names = T)

nicheConsTest(t_med_res_prep_avg_t10)
nicheConsTest(rf_med_res_prep_avg_t10)
nicheConsTest(t_med_res_prep_avg_t5)
nicheConsTest(rf_med_res_prep_avg_t5)

## QUESTION 2: IS NICHE BREADTH CONSERVED ACROSS SITES? ========

t_bdt_res_prep_t10 <- llply(.data = t_res_list_t10_subset,
      .fun = nicheConsPrep, value = "temp_breadth", otu_tab = OTUtaxonData)
rf_bdt_res_prep_t10 <- llply(.data = rf_res_list_t10_subset,
                        .fun = nicheConsPrep, value = "rf_breadth", otu_tab = OTUtaxonData)
saveRDS(t_bdt_res_prep_t10, file.path(res.dir, "t_bdt_res_prep_t10.rds"))
saveRDS(rf_bdt_res_prep_t10, file.path(res.dir, "rf_bdt_res_prep_t10.rds"))

t_bdt_res_final_t10 <- ldply(.data = t_bdt_res_prep_t10, .fun = nicheConsTest)
rf_bdt_res_final_t10 <- ldply(.data = rf_bdt_res_prep_t10, .fun = nicheConsTest)

round(mean(t_bdt_res_final_t10$cor), 2)
round(range(t_bdt_res_final_t10$cor), 2)
sum(t_bdt_res_final_t10$p < 0.05) / 100
round(range(t_bdt_res_final_t10$n), 2)

round(mean(rf_bdt_res_final_t10$cor), 2)
round(range(rf_bdt_res_final_t10$cor), 2)
sum(rf_bdt_res_final_t10$p < 0.05) / 100
round(range(rf_bdt_res_final_t10$n), 2)

# Re-run analysis with a lower site threshold
t_bdt_res_prep_t5 <- llply(.data = t_res_list_t5_subset,
                            .fun = nicheConsPrep, value = "temp_breadth", otu_tab = OTUtaxonData)
rf_bdt_res_prep_t5 <- llply(.data = rf_res_list_t5_subset,
                             .fun = nicheConsPrep, value = "rf_breadth", otu_tab = OTUtaxonData)
saveRDS(t_bdt_res_prep_t5, file.path(res.dir, "t_bdt_res_prep_t5.rds"))
saveRDS(rf_bdt_res_prep_t5, file.path(res.dir, "rf_bdt_res_prep_t5.rds"))

t_bdt_res_final_t5 <- ldply(.data = t_bdt_res_prep_t5, .fun = nicheConsTest)
rf_bdt_res_final_t5 <- ldply(.data = rf_bdt_res_prep_t5, .fun = nicheConsTest)

round(mean(t_bdt_res_final_t5$cor), 2)
round(range(t_bdt_res_final_t5$cor), 2)
sum(t_bdt_res_final_t5$p < 0.05) / 100
round(range(t_bdt_res_final_t5$n), 2)

round(mean(rf_bdt_res_final_t5$cor), 2)
round(range(rf_bdt_res_final_t5$cor), 2)
sum(rf_bdt_res_final_t5$p < 0.05) / 100
round(range(rf_bdt_res_final_t5$n), 2)

# Re-run analysis with OTU that occur in ALL rarefied datasets
t_bdt_res_prep_all_t10 <- do.call("rbind",t_bdt_res_prep_t10)
rf_bdt_res_prep_all_t10 <- do.call("rbind", rf_bdt_res_prep_t10)

t_bdt_res_prep_all_t5 <- do.call("rbind",t_bdt_res_prep_t5)
rf_bdt_res_prep_all_t5 <- do.call("rbind", rf_bdt_res_prep_t5)

t_bdt_res_prep_all_subset_t10 <- ddply(.data = t_bdt_res_prep_all_t10,
                                   .variables = .(OTU_ID),
                                   .fun = function(x) {if(nrow(x) == 100){return(x)} })
rf_bdt_res_prep_all_subset_t10 <- ddply(.data = rf_bdt_res_prep_all_t10,
                                    .variables = .(OTU_ID),
                                    .fun = function(x) {if(nrow(x) == 100){return(x)} })

t_bdt_res_prep_all_subset_t5 <- ddply(.data = t_bdt_res_prep_all_t5,
                                       .variables = .(OTU_ID),
                                       .fun = function(x) {if(nrow(x) == 100){return(x)} })
rf_bdt_res_prep_all_subset_t5 <- ddply(.data = rf_bdt_res_prep_all_t5,
                                        .variables = .(OTU_ID),
                                        .fun = function(x) {if(nrow(x) == 100){return(x)} })

t_bdt_res_prep_avg_t10 <- ddply(.data = t_bdt_res_prep_all_subset_t10,
                            .variables = .(OTU_ID),
                            .fun = summarize, 
                            Stainback = mean(Stainback),
                            Laupahoehoe = mean(Laupahoehoe),
                            V4 = V4[1])
rf_bdt_res_prep_avg_t10  <- ddply(.data = rf_bdt_res_prep_all_subset_t10,
                              .variables = .(OTU_ID),
                              .fun = summarize, 
                              Stainback = mean(Stainback),
                              Laupahoehoe = mean(Laupahoehoe),
                              V4 = V4[1])

t_bdt_res_prep_avg_t5 <- ddply(.data = t_bdt_res_prep_all_subset_t5,
                                .variables = .(OTU_ID),
                                .fun = summarize, 
                                Stainback = mean(Stainback),
                                Laupahoehoe = mean(Laupahoehoe),
                                V4 = V4[1])
rf_bdt_res_prep_avg_t5  <- ddply(.data = rf_bdt_res_prep_all_subset_t5,
                                  .variables = .(OTU_ID),
                                  .fun = summarize, 
                                  Stainback = mean(Stainback),
                                  Laupahoehoe = mean(Laupahoehoe),
                                  V4 = V4[1])

write.csv(t_bdt_res_prep_avg_t10, file.path(res.dir, "t_bdt_res_prep_avg_t10.csv"), row.names = T)
write.csv(rf_bdt_res_prep_avg_t10, file.path(res.dir, "rf_bdt_res_prep_avg_t10.csv"), row.names = T)

write.csv(t_bdt_res_prep_avg_t5, file.path(res.dir, "t_bdt_res_prep_avg_t5.csv"), row.names = T)
write.csv(rf_bdt_res_prep_avg_t5, file.path(res.dir, "rf_bdt_res_prep_avg_t5.csv"), row.names = T)

nicheConsTest(t_bdt_res_prep_avg_t10)
nicheConsTest(rf_bdt_res_prep_avg_t10)

nicheConsTest(t_bdt_res_prep_avg_t5)
nicheConsTest(rf_bdt_res_prep_avg_t5)

## QUESTION 3: IS NICHE BREADTH HIGHER ACROSS TAXONOMIC GROUPS? =====================
# Previously, this was performed with only OTUs that were found in both sites, at least 5 or 10 samples per site, and OTUs that were found in all rarefied datasets

t_res_list_t10_prep <- llply(.data = t_res_list_t10,
                             .fun = nicheConsPrep, value = "temp_breadth", otu_tab = OTUtaxonData)
t_res_list_t5_prep <- llply(.data = t_res_list_t5,
                             .fun = nicheConsPrep, value = "temp_breadth", otu_tab = OTUtaxonData)
rf_res_list_t10_prep <- llply(.data = rf_res_list_t10,
                             .fun = nicheConsPrep, value = "rf_breadth", otu_tab = OTUtaxonData)
rf_res_list_t5_prep <- llply(.data = rf_res_list_t5,
                            .fun = nicheConsPrep, value = "rf_breadth", otu_tab = OTUtaxonData)

t_res_list_t10_final <- ddply(.data = do.call("rbind",t_res_list_t10_prep),
      .variables = .(OTU_ID),
      .fun = summarize,
      Stainback = mean(Stainback, na.rm = T),
      Laupahoehoe = mean(Laupahoehoe, na.rm = T),
      V4 = V4[1])

t_res_list_t10_final <- ddply(.data = do.call("rbind",t_res_list_t10_prep),
                              .variables = .(OTU_ID),
                              .fun = summarize,
                              Stainback = mean(Stainback, na.rm = T),
                              Laupahoehoe = mean(Laupahoehoe, na.rm = T),
                              V4 = V4[1])

rf_res_list_t10_final <- ddply(.data = do.call("rbind",rf_res_list_t10_prep),
                              .variables = .(OTU_ID),
                              .fun = summarize,
                              Stainback = mean(Stainback, na.rm = T),
                              Laupahoehoe = mean(Laupahoehoe, na.rm = T),
                              V4 = V4[1])

t_res_list_t5_final <- ddply(.data = do.call("rbind",t_res_list_t5_prep),
                              .variables = .(OTU_ID),
                              .fun = summarize,
                              Stainback = mean(Stainback, na.rm = T),
                              Laupahoehoe = mean(Laupahoehoe, na.rm = T),
                              V4 = V4[1])

rf_res_list_t5_final <- ddply(.data = do.call("rbind",rf_res_list_t5_prep),
                             .variables = .(OTU_ID),
                             .fun = summarize,
                             Stainback = mean(Stainback, na.rm = T),
                             Laupahoehoe = mean(Laupahoehoe, na.rm = T),
                             V4 = V4[1])


rf_bdt_melt_t10 <- melt(rf_res_list_t10_final, id.vars = c("OTU_ID", "V4"),
     variable.name = c("Site"),
     measure.vars = c("Stainback", "Laupahoehoe"),
     value.name = c("rf_breadth"))

rf_bdt_melt_t5 <- melt(rf_res_list_t5_final, id.vars = c("OTU_ID", "V4"),
                        variable.name = c("Site"),
                        measure.vars = c("Stainback", "Laupahoehoe"),
                        value.name = c("rf_breadth"))

t_bdt_melt_t10 <- melt(t_res_list_t10_final, id.vars = c("OTU_ID", "V4"),
                       variable.name = c("Site"),
                       measure.vars = c("Stainback", "Laupahoehoe"),
                       value.name = c("t_breadth"))

t_bdt_melt_t5 <- melt(t_res_list_t5_final, id.vars = c("OTU_ID", "V4"),
                      variable.name = c("Site"),
                      measure.vars = c("Stainback", "Laupahoehoe"),
                      value.name = c("t_breadth"))

write.csv(rf_bdt_melt_t10, file = file.path(res.dir, "rf_bdt_melt_t10.csv"), row.names = FALSE)
write.csv(rf_bdt_melt_t5, file = file.path(res.dir, "rf_bdt_melt_t5.csv"), row.names = FALSE)
write.csv(t_bdt_melt_t10, file = file.path(res.dir, "t_bdt_melt_t10.csv"), row.names = FALSE)
write.csv(t_bdt_melt_t5, file = file.path(res.dir, "t_bdt_melt_t5.csv"), row.names = FALSE)

table(as.vector(rf_bdt_melt_t5$V4), rf_bdt_melt_t5$Site)
table(as.vector(rf_bdt_melt_t10$V4), rf_bdt_melt_t10$Site)
table(as.vector(t_bdt_melt_t5$V4), t_bdt_melt_t5$Site)
table(as.vector(t_bdt_melt_t10$V4), t_bdt_melt_t10$Site)


summary(lm(rf_breadth ~ V4, data = subset(rf_bdt_melt_t10, Site == "Laupahoehoe")))
summary(lm(rf_breadth ~ V4, data = subset(rf_bdt_melt_t10, Site == "Stainback")))
summary(lm(rf_breadth ~ V4, data = subset(rf_bdt_melt_t5, Site == "Laupahoehoe")))
summary(lm(rf_breadth ~ V4, data = subset(rf_bdt_melt_t5, Site == "Stainback")))

summary(lm(t_breadth ~ V4, data = subset(t_bdt_melt_t10, Site == "Laupahoehoe")))
summary(lm(t_breadth ~ V4, data = subset(t_bdt_melt_t10, Site == "Stainback")))
summary(lm(t_breadth ~ V4, data = subset(t_bdt_melt_t5, Site == "Laupahoehoe")))
summary(lm(t_breadth ~ V4, data = subset(t_bdt_melt_t5, Site == "Stainback")))


pairwise.var.test <- function (x, g, p.adjust.method = "fdr", ...) {
  # Modified from pairwise.wilcox.test
  p.adjust.method <- match.arg(p.adjust.method)
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(g)))
  g <- factor(g)
  
  compare.levels <- function(i, j) {
    xi <- x[as.integer(g) == i]
    xj <- x[as.integer(g) == j]
    var.test(xi, xj, ...)$p.value
  }
  PVAL <- pairwise.table(compare.levels, levels(g), p.adjust.method)
  PVAL2 <- pairwise.table2(compare.levels, levels(g))
  ans <- list(data.name = DNAME, p.value.adj = PVAL, p.value = PVAL2,
              p.adjust.method = p.adjust.method)
  ans
}

pairwise.table2 <- function (compare.levels, level.names, p.adjust.method) {
  ix <- setNames(seq_along(level.names), level.names)
  pp <- outer(ix[-1L], ix[-length(ix)], function(ivec, jvec) sapply(seq_along(ivec), 
                                                                    function(k) {
                                                                      i <- ivec[k]
                                                                      j <- jvec[k]
                                                                      if (i > j) 
                                                                        compare.levels(i, j)
                                                                      else NA
                                                                    }))
  #pp[lower.tri(pp, TRUE)] <- p.adjust(pp[lower.tri(pp, FALSE)], 
  #                                    p.adjust.method)
  return(pp)
}

pairwise.var.test(subset(t_bdt_melt_t5, Site == "Laupahoehoe")$t_breadth, subset(t_bdt_melt_t5, Site == "Laupahoehoe")$V4)
pairwise.var.test(subset(t_bdt_melt_t5, Site == "Stainback")$t_breadth, subset(t_bdt_melt_t5, Site == "Stainback")$V4)

pairwise.var.test(subset(rf_bdt_melt_t5, Site == "Laupahoehoe")$rf_breadth, subset(rf_bdt_melt_t5, Site == "Laupahoehoe")$V4)
pairwise.var.test(subset(rf_bdt_melt_t5, Site == "Stainback")$rf_breadth, subset(rf_bdt_melt_t5, Site == "Stainback")$V4)

pairwise.var.test(subset(t_bdt_melt_t10, Site == "Laupahoehoe")$t_breadth, subset(t_bdt_melt_t10, Site == "Laupahoehoe")$V4)
pairwise.var.test(subset(t_bdt_melt_t10, Site == "Stainback")$t_breadth, subset(t_bdt_melt_t10, Site == "Stainback")$V4)

pairwise.var.test(subset(rf_bdt_melt_t10, Site == "Laupahoehoe")$rf_breadth, subset(rf_bdt_melt_t10, Site == "Laupahoehoe")$V4)
pairwise.var.test(subset(rf_bdt_melt_t10, Site == "Stainback")$rf_breadth, subset(rf_bdt_melt_t10, Site == "Stainback")$V4)


hist(subset(t_bdt_melt_t10, Site == "Laupahoehoe")$t_breadth)
hist(subset(rf_bdt_melt_t10, Site == "Laupahoehoe")$rf_breadth)
hist(subset(t_bdt_melt_t5, Site == "Laupahoehoe")$t_breadth)
hist(subset(rf_bdt_melt_t10, Site == "Laupahoehoe")$rf_breadth)
