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

## Site patterns =========
# ggplot(data = subset(siteData, island == "BigIsland")) +
#   geom_point(aes(y = rf_ann, x = elevation)) + facet_wrap(~ site)
# 
# ggplot(data = subset(siteData, island == "BigIsland")) +
#   geom_point(aes(y = t_ann, x = elevation)) + facet_wrap(~ site)
# 
# ggplot(data = subset(siteData, island == "BigIsland")) +
#   geom_point(aes(y = t_ann, x = rf_ann)) + facet_wrap(~ site)

## Abundance patterns =========
targetOrders <- c("Acari", "Araneae", "Coleoptera", "Lepidoptera", "Orthoptera", "Neuroptera", "Psocoptera", "Hemiptera")
otu_target <- which(taxonData$V4[match(rownames(OTU_mat), taxonData$OTU_ID)] %in% targetOrders)
OTU_df <- melt(data = OTU_mat[otu_target,], value.name = "nReads", varnames = c("OTU_ID", "site.id"))
OTU_clim_df <- merge(x = OTU_df, y = siteData[c("site","site.id", "rf_ann", "t_ann")], by = "site.id", all.x = TRUE)

# Remove OTUs with fewer than 10 sites
source(file.path(main.dir, "test.R"))
removeLowOcc <- function(data){
  if(sum(data$nReads > 0) >= 10){
    return(data)
  }
}
OTU_clim_df2 <- ddply(.data = OTU_clim_df, .variables = .(site, OTU_ID), .fun = removeLowOcc)

# Fit a gaussian distribution to each OTU in each site
OTU_gauss <- dlply(.data = OTU_clim_df2, .variables = .(site, OTU_ID), .fun = fitGaussianModel2, .progress = "text")
OTU_gauss_para <- do.call("rbind", lapply(X = OTU_gauss, FUN = function(x) x[[1]] ))

# Infer true ranges using gaussian parameters
OTU_gauss_q <- ddply(OTU_gauss_para, .variables = .(site,OTU_ID), .fun = calcQuantiles)
OTU_gauss_all <- merge(OTU_gauss_para, OTU_gauss_q)

# Flag OTUs whose range boundaries are outside the sampling ranges for each site
OTU_gauss_all$rf_bounds <- with(data = OTU_gauss_all,
                                ifelse(site == "Laupahoehoe" & rf_q3 > 4402.763 |
                                       site == "Laupahoehoe" & rf_q1 < 2840.776 | 
                                       site == "Stainback" & rf_q3 > 6388.352 |
                                       site == "Stainback" & rf_q1 < 2480.524, 0, 1))

OTU_gauss_all$t_bounds <- with(data = OTU_gauss_all, 
                               ifelse(site == "Laupahoehoe" & t_q3 > 17.89003 |
                                      site == "Laupahoehoe" &  t_q1 < 13.54576 |
                                      site == "Stainback" & t_q3 > 18.16707 |
                                      site == "Stainback" & t_q1 < 12.28192, 0, 1))

table(OTU_gauss_all$rf_bounds)
table(OTU_gauss_all$t_bounds)

tapply(X = siteData$rf_ann, INDEX = siteData$site, FUN = range)
tapply(X = siteData$t_ann, INDEX = siteData$site, FUN = range)

# Only include species whose 10-90 range remains within sampling bounds for each gradient
OTU_gauss_t_subset <- subset(OTU_gauss_all, t_bounds == 1)
OTU_gauss_rf_subset <- subset(OTU_gauss_all, rf_bounds == 1)

# Only include species that are found in both sites
foobar <- function(x){
  if("Laupahoehoe" %in% x$site & "Stainback" %in% x$site){
    return(x)
  }
}
OTU_gauss_t_subset <- ddply(.data = OTU_gauss_t_subset, .variables = .(OTU_ID), .fun = foobar)
OTU_gauss_rf_subset <- ddply(.data = OTU_gauss_rf_subset, .variables = .(OTU_ID), .fun = foobar)

# Merge back with read abundance data
OTU_gauss_t_subset_reads <- merge(OTU_gauss_t_subset, OTU_clim_df2[c("site.id", "OTU_ID", "site", "nReads")])
OTU_gauss_rf_subset_reads <- merge(OTU_gauss_rf_subset, OTU_clim_df2[c("site.id", "OTU_ID", "site", "nReads")])

OTU_gauss_subset_t_laup_mat <- acast(subset(OTU_gauss_t_subset_reads, site == "Laupahoehoe"),
                                     formula = site.id~OTU_ID, fill = 0, value.var = "nReads")
OTU_gauss_subset_t_stbk_mat <- acast(subset(OTU_gauss_t_subset_reads, site == "Stainback"),
                                     formula = site.id~OTU_ID, fill = 0, value.var = "nReads")
OTU_gauss_subset_rf_laup_mat <- acast(subset(OTU_gauss_rf_subset_reads, site == "Laupahoehoe"),
                                      formula = site.id~OTU_ID, fill = 0, value.var = "nReads")
OTU_gauss_subset_rf_stbk_mat <- acast(subset(OTU_gauss_rf_subset_reads, site == "Stainback"),
                                      formula = site.id~OTU_ID, fill = 0, value.var = "nReads")

# Calculate observed differences in median and range and test them against a null model
rf_test <- elevationNullModel(mat_x = OTU_gauss_subset_rf_laup_mat,
                           mat_y = OTU_gauss_subset_rf_stbk_mat,
                           site.data = siteData,
                           nreps = 999,
                           env.vars  = "rf_ann",
                           site.var = "site.id")

t_test <- elevationNullModel(mat_x = OTU_gauss_subset_t_laup_mat,
                              mat_y = OTU_gauss_subset_t_stbk_mat,
                              site.data = siteData,
                              nreps = 999,
                              env.vars  = "t_ann",
                              site.var = "site.id")

# p-value = the proportion of simulations that are HIGHER than the observed difference
# comparison is values in x minus y. Therefore z scores are higher if laupahoehoe is values in higher. 
# x = laupahoehoe
# y = stainback

# Plot results
meh <- names(rf_test)[grep(names(rf_test), pattern = "_x|_y")]
test_melt <- melt(rf_test,  measure.vars = meh)

ggplot(data = subset(test_melt, variable %in% c("rf_ann_range_x", "rf_ann_range_y"))) +
  geom_point(aes(y = value, x = variable)) +
  geom_path(aes(y = value, x= variable, group = ID, alpha = rf_ann_range_diff_p ))

ggplot(data = subset(test_melt, variable %in% c("rf_ann_median_x", "rf_ann_median_y"))) +
  geom_point(aes(y = value, x = variable)) +
  geom_path(aes(y = value, x= variable, group = ID, alpha = ifelse(rf_ann_median_diff_p < 0.05, 1, 0)))

# Plot temperature results
meh <- names(t_test)[grep(names(t_test), pattern = "_x|_y")]
t_test_melt <- melt(t_test,  measure.vars = meh)

ggplot(data = subset(t_test_melt, variable %in% c("t_ann_range_x", "t_ann_range_y"))) +
  geom_point(aes(y = value, x = variable)) +
  geom_path(aes(y = value,
                x= variable,
                group = ID,
                alpha = t_ann_range_diff_p))
library(lme4)
summary(lm(value ~ variable, data = subset(t_test_melt, variable %in% c("t_ann_range_x", "t_ann_range_y"))))
t.test(value~variable, data = subset(t_test_melt, variable %in% c("t_ann_range_x", "t_ann_range_y")), paired = T)
t.test(value~variable, data = subset(t_test_melt, variable %in% c("t_ann_range_x", "t_ann_range_y")), paired = F)

ggplot(data = subset(t_test_melt, variable %in% c("t_ann_median_x", "t_ann_median_y"))) +
  geom_point(aes(y = value, x = variable)) +
  geom_path(aes(y = value,
                x= variable,
                group = ID,
                alpha = ifelse(t_ann_median_diff_p < 0.05, 1,0)))

y <- subset(test_melt, variable %in% c("t_ann_range_x", "t_ann_range_y"))


#y$ID %in% 

head(test_melt)



colSums(OTU_gauss_subset_rf_stbk_mat)
str(randomizeMatrix(OTU_gauss_subset_rf_stbk_mat, null.model = "frequency"))

table(siteData$site)
dim(acast(data = OTU_gauss_subset_t_laup_mat, formula = OTU_ID~site.id, fill = 0))
head(OTU_gauss_subset_t_laup_mat)
subset(OTU_gauss_subset_t_laup_mat, OTU_ID == "OTU59")[c("site.id", "nReads")]
library(picante)
randomizeMatrix(samp = mat, null.model = "frequency", iterations = 100)
# site against species



# Identify OTUs with more than 5 occurrences at each locality
# otu_laup_nsites <- which(rowSums(ifelse(OTU_mat[,laupahoehoe_siteIDs] > 0, 1, 0)) > 5)
# otu_stbk_nsites <- which(rowSums(ifelse(OTU_mat[,stainback_siteIDs] > 0, 1, 0)) > 5)


# Identify OTUs whose ranges are close to the edges of sampled boundaries

  

targetOTU <- Reduce(list(otu_laup_nsites, otu_stbk_nsites, otu_target), f = intersect) 


#OTU_df <- melt(OTU_mat, value.name = "rarefiedReadAbund", varnames = c("OTU_ID", "site.id"))
OTU_df2 <- merge(x = OTU_df, y = taxonData[c("OTU_ID", "V4")], by = "OTU_ID", all.x = TRUE)

# Remove duplicates
OTU_df2 <- OTU_df2[!duplicated(OTU_df2),]
OTU_df3 <- merge(x = OTU_df2, y = siteData[c("site.id", "rf_ann", "t_ann", "site")], by = "site.id")

# Only include species that occur in more than 5 sites.
determineSuitability <- function(x){
  suitable <- ifelse(sum(x$rarefiedReadAbund > 0) >= 5, "yes", "no")
  return(data.frame(suitable))
}

OTU_df_suit <- ddply(.data = OTU_df3, .fun = determineSuitability, .variables = .(site,OTU_ID))
OTU_df_suit2 <- merge(OTU_df3, OTU_df_suit, by = c("site", "OTU_ID"))
OTU_df_suit2 <- subset(OTU_df_suit2, rarefiedReadAbund > 0)

## Species diversity =========
#OTU_df3 <- subset(OTU_df2, !V4 %in% c("Amphipoda", "Diptera", "Hymenoptera") & rarefiedReadAbund > 0)

# speciescount <- ddply(OTU_df3, .variables = .(site.id), .fun = function(x){
#   nsp = length(unique(as.vector(x$OTU_ID)))
#   simpson = sum((x$rarefiedReadAbund / sum(x$rarefiedReadAbund))^2)
#   shannon = sum((x$rarefiedReadAbund / sum(x$rarefiedReadAbund)) * log((x$rarefiedReadAbund / sum(x$rarefiedReadAbund))))
#   data.frame(nsp, simpson, shannon)})
# speciescountbysite <- merge(speciescount, siteData)
# ggplot(aes(y = nsp, x = elevation), data = speciescountbysite) + geom_point() + geom_smooth(method = "loess") + facet_wrap(~ site)
# ggplot(aes(y = simpson, x = elevation), data = speciescountbysite) + geom_point() + geom_smooth(method = "loess") + facet_wrap(~ site)
# ggplot(aes(y = -1*shannon, x = elevation), data = speciescountbysite) + geom_point() + geom_smooth(method = "loess") + facet_wrap(~ site)
# 
# OTU_mat <- acast(site.id~OTU_ID, data = OTU_df3, value.var = "rarefiedReadAbund", fill = 0)
# OTU_mds <- metaMDS(comm = OTU_mat, dist = "bray", k = 3)
# OTU_mds_coords <-  as.data.frame(OTU_mds$points)
# OTU_mds_coords$site.id <- rownames(OTU_mds_coords)
# OTU_mds_coords_final <- merge(OTU_mds_coords, siteData, by = "site.id")

# Calculate OTU climatic ranges and optima by site ======
# Define function to find climatic ranges for each OTU
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

# Calculate the OTU climatic ranges and optima
OTU_clim_site <- ddply(.data = subset(OTU_df_suit2, suitable == "yes"), .variables = .(site, OTU_ID), .fun = findLimits)
OTU_clim_glob  <- ddply(.data = subset(OTU_df_suit2, suitable == "yes"), .variables = .(OTU_ID), .fun = findLimits)

# Flag species that are most abundant at the sampling limits
site_rf_lim <- with(siteData, tapply(rf_ann, INDEX = site, FUN = range))
site_t_lim <- with(siteData, tapply(t_ann, INDEX = site, FUN = range))

OTU_clim_site$rf_range <- OTU_clim_site$rf_upper - OTU_clim_site$rf_lower
OTU_clim_site$rf_site_upperlimit <-
  ifelse(OTU_clim_site$site == "Laupahoehoe",
         ifelse(OTU_clim_site$rf_upper == site_rf_lim$Laupahoehoe[2], 1, 0),
         ifelse(OTU_clim_site$rf_upper == site_rf_lim$Stainback[2], 1, 0))
OTU_clim_site$rf_site_lowerlimit <-
  ifelse(OTU_clim_site$site == "Laupahoehoe",
         ifelse(OTU_clim_site$rf_lower == site_rf_lim$Laupahoehoe[1], 1, 0),
         ifelse(OTU_clim_site$rf_lower == site_rf_lim$Stainback[1], 1, 0))

OTU_clim_site$temp_range <- OTU_clim_site$temp_upper - OTU_clim_site$temp_lower
OTU_clim_site$temp_site_upperlimit <-
  ifelse(OTU_clim_site$site == "Laupahoehoe",
         ifelse(OTU_clim_site$temp_upper == site_t_lim$Laupahoehoe[2], 1, 0),
         ifelse(OTU_clim_site$temp_upper == site_t_lim$Stainback[2], 1, 0))
OTU_clim_site$temp_site_lowerlimit <-
  ifelse(OTU_clim_site$site == "Laupahoehoe",
         ifelse(OTU_clim_site$temp_lower == site_t_lim$Laupahoehoe[1], 1, 0),
         ifelse(OTU_clim_site$temp_lower == site_t_lim$Stainback[1], 1, 0))

OTU_clim_site_subset <- subset(OTU_clim_site, !order %in% c("Amphipoda", "Diptera", "Hymenoptera"))

## Analyze
laupahoehoueOTUlist <- unique(subset(OTU_clim_site_subset, site == "Laupahoehoe" &
                                       (rf_site_lowerlimit==0 & rf_site_upperlimit == 0)))$OTU_ID
stainbackOTUlist <- unique(subset(OTU_clim_site_subset, site == "Stainback" &
                                       (rf_site_lowerlimit==0 & rf_site_upperlimit == 0)))$OTU_ID

site_rf_OTU <- intersect(as.vector(stainbackOTUlist), as.vector(laupahoehoueOTUlist))

OTU_clim_site_rf_subset <- subset(OTU_clim_site_subset, OTU_ID %in% site_rf_OTU)

# Rules
# Have to 
# Quantile of species is neither in the upper or lower limit for both sites
# Species is in both sites
site_rf_OTU
(OTU_mat)



ggplot(aes(y = rf_range, x = site, group = OTU_ID), data = OTU_clim_site_rf_subset) + 
  geom_point(aes(color = order)) + 
  geom_path(aes(color = order)) +
  scale_color_brewer(palette = "Spectral") + 
  labs(y = "Mean annual rainfall range", x = "Site") +
  theme(legend.position = "bottom",legend.justification="center")








temp_r_site_plot <- ggplot(data = subset(OTU_clim_site_subset,!(temp_site_lowerlimit==1|temp_site_upperlimit == 1))) +
  geom_boxplot(aes(y = temp_range, x = site, fill = order)) + 
  scale_fill_brewer(palette = "Spectral") + 
  guides(fill = guide_legend(title = "Order")) +
  labs(y = "Mean annual temperature range", x = "Site")

rf_r_site_plot <- ggplot(data = subset(OTU_clim_site_subset,!(rf_site_lowerlimit==1|rf_site_upperlimit == 1))) +
  geom_boxplot(aes(y = rf_range, x = site, fill = order)) + 
  scale_fill_brewer(palette = "Spectral") + 
  guides(fill = guide_legend(title = "Order")) +
  labs(y = "Mean annual rainfall range", x = "Site") +
  theme(legend.position = "bottom",legend.justification="center")

order_leg <- get_legend(rf_r_site_plot)

# Plot climatic ranges for all OTUs by site and order
library(cowplot)
clim_plot_woleg <- plot_grid(temp_r_site_plot + theme(legend.position = "none"),
                             rf_r_site_plot + theme(legend.position = "none"), labels = "AUTO")
clim_plot_wleg <- plot_grid(clim_plot_woleg, order_leg, nrow = 2, rel_heights = c(4,1))
ggsave(file.path(fig.dir, "climrange.pdf"), clim_plot_wleg, width = 8, height = 5)

## Niche conservatism across sites
OTU_tempMcons <- dcast( subset(OTU_clim_site_subset, !(temp_site_lowerlimit==1 | temp_site_upperlimit == 1)), formula = OTU_ID ~ site, value.var = "temp_median")
OTU_tempMcons2 <- merge(OTU_tempMcons, OTUtaxonData)

OTU_rfMcons <- dcast( subset(OTU_clim_site_subset, !(temp_site_lowerlimit==1 | temp_site_upperlimit == 1)), formula = OTU_ID ~ site, value.var = "rf_median")
OTU_rfMcons2 <- merge(OTU_rfMcons, OTUtaxonData)

tempMcons_plot <- ggplot(data = na.omit(OTU_tempMcons2)) + geom_point(aes(y = Laupahoehoe, x= Stainback, fill = V4), pch = 21, size = 2) + geom_smooth(aes(y = Laupahoehoe, x = Stainback), method = "lm", col = "grey70") + 
  guides(fill = guide_legend(title = "Order", override.aes = list(size = 3))) + 
  labs(y = "Median annual temperature\n(Laupahoehoe)", x = "Median annual temperature\n(Stainback)") +
  scale_fill_brewer(palette = "Spectral") +
  theme(legend.position = "bottom")

rfMcons_plot <- ggplot(data = na.omit(OTU_rfMcons2)) + geom_point(aes(y = Laupahoehoe, x= Stainback, fill = V4), pch = 21, size = 2) + geom_smooth(aes(y = Laupahoehoe, x = Stainback), method = "lm", col = "grey70") + 
  guides(fill = guide_legend(title = "Order", override.aes = list(size = 3))) + 
  labs(y = "Median annual precipitation\n(Laupahoehoe)", x = "Median annual precipitation\n(Stainback)") +
  scale_fill_brewer(palette = "Spectral") +
  theme(legend.position = "bottom",
        legend.justification = "center")

mediancons_leg <- get_legend(rfMcons_plot)

mediancons_comb <- plot_grid(plot_grid(tempMcons_plot + theme(legend.position = "none"),
                                       rfMcons_plot + theme(legend.position= "none"), labels = "AUTO",
                                       nrow = 1, rel_widths = c(1,1), align = "v"),
                             mediancons_leg,
                             nrow = 2, rel_heights = c(1, 0.1))
ggsave(file.path(fig.dir, "medianconservatism.pdf"), mediancons_comb, width = 8, height = 5)

## Range conservatism ===============
OTU_tempRcons <- dcast( subset(OTU_clim_site_subset, !(temp_site_lowerlimit==1 | temp_site_upperlimit == 1)), formula = OTU_ID ~ site, value.var = "temp_range")
OTU_tempRcons2 <- merge(OTU_tempRcons, OTUtaxonData)

OTU_rfRcons <- dcast( subset(OTU_clim_site_subset, !(temp_site_lowerlimit==1 | temp_site_upperlimit == 1)), formula = OTU_ID ~ site, value.var = "rf_range")
OTU_rfRcons2 <- merge(OTU_rfRcons, OTUtaxonData)

tempRcons_plot <- ggplot(data = na.omit(OTU_tempRcons2)) + geom_point(aes(y = Laupahoehoe, x= Stainback, fill = V4), pch = 21, size = 2) + geom_smooth(aes(y = Laupahoehoe, x = Stainback), method = "lm", col = "grey70") + 
  guides(fill = guide_legend(title = "Order", override.aes = list(size = 3))) + 
  labs(y = "Range in median annual temperature\n(Laupahoehoe)",
       x = "Range in median annual temperature\n(Stainback)") +
  scale_fill_brewer(palette = "Spectral") +
  theme(legend.position = "bottom")

rfRcons_plot <- ggplot(data = na.omit(OTU_rfRcons2)) + geom_point(aes(y = Laupahoehoe, x= Stainback, fill = V4), pch = 21, size = 2) + geom_smooth(aes(y = Laupahoehoe, x = Stainback), method = "lm", col = "grey70") + 
  guides(fill = guide_legend(title = "Order", override.aes = list(size = 3))) + 
  labs(y = "Range in annual precipitation\n(Laupahoehoe)", x = "Range in annual precipitation\n(Stainback)") +
  scale_fill_brewer(palette = "Spectral") +
  theme(legend.position = "bottom",
        legend.justification = "center")

rangecons_leg <- get_legend(rfRcons_plot)

rangecons_comb <- plot_grid(plot_grid(tempRcons_plot + theme(legend.position = "none"),
                                       rfRcons_plot + theme(legend.position= "none"), labels = "AUTO",
                                       nrow = 1, rel_widths = c(1,1), align = "v"),
                             rangecons_leg,
                             nrow = 2, rel_heights = c(1, 0.1))
ggsave(file.path(fig.dir, "rangeconservatism.pdf"), rangecons_comb, width = 8, height = 5)

## Calculate climatic optima across both sites ========

#geneticDist <- read.csv(file.path(data.dir, "K2P_Distance_COI.csv"))
geneticDistOTU <- readRDS(file.path(data.dir, "geneticDist_OTU_mat.rds"))
taxonRefOTUs <- read.csv(file.path(data.dir, "taxoRefOTUs.csv"))

OTU_genDist <- melt(as.matrix(geneticDistOTU), value.name = "genDist")

rownames(OTU_clim_glob) <- OTU_clim_glob$OTU_ID

rf_upperlim <- OTU_clim_glob$OTU_ID[OTU_clim_glob$rf_upper %in% c(4402.763,6388.352)]
rf_lowerlim <- OTU_clim_glob$OTU_ID[OTU_clim_glob$rf_upper %in% c(2480.524,2840.776) ]
temp_upperlim <- OTU_clim_glob$OTU_ID[OTU_clim_glob$rf_upper %in% c(18.16707, 17.89003)]
temp_lowerlim <- OTU_clim_glob$OTU_ID[OTU_clim_glob$rf_upper %in% c(12.28192, 13.54576)]
OTU_lim_list <- Reduce(list(rf_upperlim, rf_lowerlim, temp_upperlim, temp_lowerlim), f = union)

OTU_temp_diff <- melt(as.matrix(dist(OTU_clim_glob["temp_median"])), value.name = "temp_median_diff")
OTU_rf_diff <- melt(as.matrix(dist(OTU_clim_glob["rf_median"])), value.name = "rf_median_diff")
(OTU_temp_diff)
OTU_diff <- Reduce(list(OTU_genDist, OTU_temp_diff, OTU_rf_diff), f = merge)

OTU_diff2 <- subset(OTU_diff, !Var1 == Var2)
OTU_diff3 <- merge(y = taxonRefOTUs[c("V3", "Genus", "Family", "Order")], x = OTU_diff2, by.y = "V3", by.x = "Var1", suffixes = "Var1")
OTU_diff4 <- merge(y = taxonRefOTUs[c("V3", "Genus", "Family", "Order")], x = OTU_diff3, by.y = "V3", by.x = "Var2", suffixes = "Var2")
names(OTU_diff4)

names(OTU_diff4)

nichecons <- ggplot(aes(x = genDist, y = temp_median_diff),
                    data= subset(OTU_diff4, FamilyNA == FamilyVar2 & (!OrderVar2 %in% c("Amphipoda", "Diptera", "Hymenoptera")))) + geom_point(aes(fill = OrderVar2), pch = 21, size = 2) + facet_wrap(~OrderVar2, nrow = 2) + geom_smooth(colour = "black", method = "lm") +   scale_fill_brewer(palette = "Spectral")
ggsave(file.path(fig.dir, "nicheevol_temp.pdf"), nichecons, width = 10, height = 5)
#GenusNA == GenusVar2 &
head(OTU_diff4)

nichecons2 <- ggplot(aes(x = genDist, y = rf_median_diff),
                     data= subset(OTU_diff4,
                                  GenusNA == GenusVar2 & 
                                    (!GenusNA == "") &
                                    (!OrderVar2 %in% c("Amphipoda", "Diptera", "Hymenoptera")) & 
                                    (! Var1 %in% OTU_lim_list) & 
                                    (! Var2 %in% OTU_lim_list))) + geom_point(aes(fill = OrderVar2), pch = 21, size = 2) + facet_wrap(~OrderVar2, nrow = 2) + geom_smooth(colour = "black", method = "lm") +   scale_fill_brewer(palette = "Spectral")
ggsave(file.path(fig.dir, "nicheevol_rf.pdf"), nichecons2, width = 10, height = 5)




## NOTES:
# Fig 1: Shared OTUs, Range L vs Range S, a) Temp, b) Precipitation
# Fig 2: All OTUs Box pt Range L vs Range S
# Fig 3: Shared OTUs, Mean t vs Mean t , Mean rf vs Mean rf
# Fig 4: Niche conservatism (using the global mean and )