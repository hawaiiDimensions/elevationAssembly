## Calculate haplotype diversity at each site =========
zotu_native_files <- list.files(file.path(data.dir, "zOTU_native"))
zOTU_mat_list <- lapply(zotu_native_files, FUN = function(x){readRDS(file.path(data.dir, "zOTU_native", x))})
zOTU_mat_df <- melt(zOTU_mat_list, varnames = c("zOTU_ID", "site.id"),
                    value.name = "nReads")

zOTU_mat_df_taxon <- merge(zOTU_mat_df, taxonData[c("OTU_ID", "zOTU_ID", "V4")], by = "zOTU_ID")

zOTU_mat_df2 <- subset(zOTU_mat_df_taxon, V4 %in% c("Araneae", "Coleoptera", "Hemiptera", "Lepidoptera", "Psocoptera", "Orthoptera"))

zOTU_mat_haplodiv <- ddply(subset(zOTU_mat_df2, nReads > 0), .variables = .(L1, OTU_ID, site.id),
                           .fun = function(x) {data.frame( nhaplotype = length(as.vector(unique(x$zOTU_ID))))},
                           .progress = "text")

zOTU_mat_haplodiv2 <- merge(x = zOTU_mat_haplodiv, y = siteSpecimenCounts, by.x = "site.id", by.y = "Site_ID", all.x = TRUE)
saveRDS(zOTU_mat_haplodiv2, file.path(res.dir,"zOTU_mat_haplodiv.rds"))

## Calculate average read count across rarefied datasets
zOTU_avg_df <- ddply(subset(zOTU_mat_df2, nReads > 0), .variables = .(site.id, OTU_ID, zOTU_ID),
                     summarize,
                     avgRead = mean(nReads),
                     .progress = "text")

zOTU_avg_df2 <- merge(zOTU_avg_df, siteData[c("site.id", "site")])
zOTU_avg_transect <- ddply(zOTU_avg_df2,
                           .variable = .(site, OTU_ID), summarize,
                           nHaplotype = length(unique(zOTU_ID)),
                           .progress = "text")

ggplot(data = zOTU_avg_transect) + 
  geom_boxplot(aes(x = site, y = log(nHaplotype)))
#geom_point(aes(x = site, y = log(nHaplotype), colour = OTU_ID)) + 
#geom_path(aes(y = log(nHaplotype), x = site, group = OTU_ID, colour = OTU_ID)) + theme(legend.position = "none")

zOTU_avg_transect_df <- dcast(zOTU_avg_transect, formula = OTU_ID~site)
zOTU_avg_transect_df <- zOTU_avg_transect_df[complete.cases(zOTU_avg_transect_df),]

t.test(x = zOTU_avg_transect_df$Laupahoehoe, y = zOTU_avg_transect_df$Stainback, paired= TRUE)

library(lme4)
summary(glmer(nHaplotype~ site + (1|OTU_ID), data = zOTU_avg_transect, family = "poisson"))

## Beta-diversity within OTUs =========================
# Import zotu tables
zotu_native_files <- list.files(file.path(data.dir, "zOTU_native"))
zOTU_mat_list <- lapply(zotu_native_files, FUN = function(x){readRDS(file.path(data.dir, "zOTU_native", x))})

# Collapse rarefied zotu tables into a single dataframe
zOTU_mat_df <- melt(zOTU_mat_list, varnames = c("zOTU_ID", "site.id"),
                    value.name = "nReads")
zOTU_mat_df_subset <- subset(zOTU_mat_df, !site.id %in% lowsiteids) # remove sites below 800 m elevation

# Merge with taxonomic infomration
zOTU_mat_df_taxon <- merge(zOTU_mat_df_subset, taxonData[c("OTU_ID", "zOTU_ID", "V4")], by = "zOTU_ID")

zOTU_mat_df2 <- subset(zOTU_mat_df_taxon, V4 %in% c("Araneae", "Coleoptera", "Hemiptera", "Lepidoptera", "Psocoptera", "Orthoptera"))

# Calculate alpha-diversity in zOTUs ====================
zOTU_mat_haplodiv <- ddply(subset(zOTU_mat_df2, nReads > 0), .variables = .(L1, OTU_ID, site.id),
                           .fun = function(x) {data.frame( nhaplotype = length(as.vector(unique(x$zOTU_ID))))},
                           .progress = "text")
saveRDS(zOTU_mat_haplodiv, file.path(res.dir,"zOTU_mat_haplodiv.rds"))

## Calculate beta-diversity in zOTUs ====================
# Average across rarefied tables
zOTU_mat_avg_df <- ddply(zOTU_mat_df2,
                         .variables = .(site.id, zOTU_ID),
                         summarize, 
                         avgRead = mean(nReads), 
                         V4 = V4[1], 
                         OTU_ID = OTU_ID[1],
                         progress = "text")

# Only include OTUs with more than 1 zOTU
zOTU_taxon_df_subset <- ddply(subset(zOTU_taxon_df, avgRead > 0),
                              .variables = .(OTU_ID), .fun = function(x){ if(length(unique(x$zOTU_ID)) > 1){ return(x) } else {return(NULL)}  })

# Only include OTUs with more than 1 site
zOTU_taxon_df_subset <- ddply(zOTU_taxon_df_subset, .variables= .(OTU_ID),
                              .fun = function(x){
                                if(length(unique(as.vector(x$site_id))) > 1){
                                  return(x)
                                } else {
                                  return(NULL)
                                }} )

# Calculate pairwise genetic distance
zOTU_gendist <- ddply(zOTU_taxon_df_subset, .variable = "OTU_ID", .fun = function(x) {
  x_mat <- reshape2::acast(x, formula = site_id ~ zOTU_ID, fill = 0, value.var = "avgRead")
  seqData_OTU <- seqData[names(seqData) %in% colnames(x_mat)]
  zotu_dist <- dist.dna(seqData_OTU, as.matrix = T)
  max_zotu_dist <- max(zotu_dist)
  site_site_avgdist <- matrix(nrow = nrow(x_mat), ncol = nrow(x_mat))
  rownames(site_site_avgdist) <- colnames(site_site_avgdist) <- rownames(x_mat)
  site_site_stddist <- matrix(nrow = nrow(x_mat), ncol = nrow(x_mat))
  rownames(site_site_stddist) <- colnames(site_site_stddist) <- rownames(x_mat)
  for(i in 1:nrow(x_mat)){
    for(j in 1:nrow(x_mat)){
      pw_comparisons <- x_mat[i,]  %*% t(x_mat[j,])
      total_dist <- sum(pw_comparisons * zotu_dist)
      avg_dist <- total_dist / sum(pw_comparisons)
      std_avg_dist <- avg_dist / max_zotu_dist
      site_site_avgdist[i,j] <- avg_dist
      site_site_stddist[i,j] <- std_avg_dist
    }
  }
  res <- merge(x= melt(site_site_avgdist, value.name = "avgdist"),
               y = melt(site_site_stddist, value.name = "stddist"),
               by = c("Var1","Var2"))
  return(res)
})

# Calculate beta-diversity of zOTUs for each OTU
zOTU_beta <- ddply(.data = zOTU_taxon_df_subset,
                   .variables = .(OTU_ID),
                   .fun = function(x){
                     #x = subset(zOTU_taxon_df_subset, OTU_ID == "OTU100")
                     OTU_mat <- reshape2::acast(x, formula = site_id~zOTU_ID, value.var = "avgRead", fill = 0)
                     OTU_mat2 <- OTU_mat[rowSums(OTU_mat) > 0,]
                     beta_mat <- as.matrix(vegdist(OTU_mat2, method = "bray", diag = T, upper = T))
                     distDF <- melt(beta_mat, value.name = "haplo_bray")
                     return(distDF)
                   })

# Combine average pairwise distance and beta diversity
zOTU_beta_combined <- merge(zOTU_gendist, zOTU_beta, by = c("OTU_ID", "Var1", "Var2"), all = T)

# Calculate elevation differences between sites
siteDiss <- siteData[c("elevation")]
rownames(siteDiss) <- siteData$site.id
siteElevDist <- melt(as.matrix(dist(siteDiss)), value.name = "elev_diff")

zOTU_beta_final <- 
  zOTU_beta_combined %>% 
  merge(y = siteElevDist, by = c("Var1", "Var2"), all.x = TRUE) %>%
  merge(y = siteData[,c("site.id", "elevation", "site")], by.x = "Var1", by.y = "site.id") %>%
  merge(y = siteData[,c("site.id", "elevation", "site")], by.x = "Var2", by.y = "site.id", suffixes = c("_X1", "_X2"))

# Perform mantel test for each OTU with elevation
laupahoehoe_zOTU_beta_mantel <- ddply(.data = subset(zOTU_beta_final, site_X1 == "Stainback" & site_X1 == site_X2),
                                      .variables = .(OTU_ID), 
                                      .fun = function(x){
                                        if(nrow(x) < 25){
                                          # i.e. OTU is detected in < 5 samples in the transect
                                          return(NULL)
                                        } else {
                                          stddist <- acast(x, formula = Var1 ~ Var2, value.var = "stddist")
                                          elevdist <- acast(x, formula = Var1 ~ Var2, value.var = "elev_diff")
                                          haplobray <- acast(x, formula = Var1 ~ Var2, value.var = "haplo_bray")
                                          mantel_gendist <- mantel(stddist, elevdist)
                                          gendist_r <- mantel_gendist$statistic
                                          gendist_p <- mantel_gendist$signif
                                          mantel_haplobray <- mantel(haplobray, elevdist)
                                          haplobray_r <- mantel_haplobray$statistic
                                          haplobray_p <- mantel_haplobray$signif
                                          data.frame(gendist_r, gendist_p, haplobray_r, haplobray_p)
                                        }
                                      })
sum(laupahoehoe_zOTU_beta_mantel$gendist_p < 0.05, na.rm = TRUE)
sum(laupahoehoe_zOTU_beta_mantel$gendist_p >= 0.05, na.rm = TRUE)

sum(laupahoehoe_zOTU_beta_mantel$haplobray_p < 0.05, na.rm = TRUE)
sum(laupahoehoe_zOTU_beta_mantel$haplobray_p >= 0.05, na.rm = TRUE)

ggplot(aes(x = elev_diff, y = stddist), data = test) + geom_point() + geom_smooth(method = "lm")

# Calculate average distance across OTUs
zOTU_beta_avg <- ddply(zOTU_beta_combined, .variables = .(Var2,Var1),
                       .fun = function(x){
                         nSharedOTUs = length(x$OTU_ID) # number of OTUs compared between sites
                         mean_gendist = mean(x$avgdist, na.rm = T)
                         mean_stddist = mean(x$stddist, na.rm = T)
                         mean_haplo_bray = mean(x$haplo_bray, na.rm = T)
                         data.frame(nSharedOTUs, mean_gendist, mean_stddist, mean_haplo_bray)
                       })

# Remove pairwise duplicates
zOTU_beta_avg_sort = t(apply(zOTU_beta_avg[,c("Var1", "Var2")], 1, sort))
zOTU_beta_avg_unique <- zOTU_beta_avg[!duplicated(zOTU_beta_avg_sort), ]
zOTU_beta_avg_unique <- subset(zOTU_beta_avg_unique, !Var1 == Var2)

# Merge with site data
zOTU_beta_avg_unique2 <- merge(x = zOTU_beta_avg_unique, y = siteData[, c("site.id", "elevation", "site")], by.x = "Var1", by.y = "site.id")
zOTU_beta_avg_unique3 <- merge(x = zOTU_beta_avg_unique2, y = siteData[, c("site.id", "elevation", "site")], by.x = "Var2", by.y = "site.id", suffixes = c("_X1", "_X2"))

zOTU_beta_avg_unique3$elevation_mean <- apply(zOTU_beta_avg_unique3[,c("elevation_X1", "elevation_X2")], 1, mean )
write.csv(zOTU_beta_avg_unique3, file = file.path(res.dir, "zOTU_beta_results.csv"))

## GRAPHICS
## Haplotype diversity ================
site_haplodiv <- readRDS(file.path(res.dir, "zOTU_mat_haplodiv.rds"))

## Genetic distance ================
zOTU_beta_distance <- read.csv(file.path(res.dir, "zOTU_beta_results.csv"))

zOTU_beta_distance_subset <- subset(zOTU_beta_distance, abs(elevation_X2 - elevation_X1) <= 100 & site_X1 != site_X2)

zOTU_betsite_bray <- ggplot(data = zOTU_beta_distance_subset, aes(y = mean_haplo_bray, x= elevation_mean)) + 
  geom_point(aes( size = nSharedOTUs), pch = 21, fill= "white") +
  labs(y = "Pairwise Bray-Curtis haplotype dissimilarity", x = "Mean elevation") +
  geom_smooth(method = "lm", colour = "grey20") +
  scale_radius(name = "Number of shared OTUs", range = c(1,6), breaks = c(40, 50, 60, 70, 80), labels = c(40, 50, 60, 70, 80), limits = c(40,80)) +
  theme(panel.background = element_blank())

zOTU_betsite_stddist <- ggplot(data = zOTU_beta_distance_subset, aes(y = mean_stddist, x= elevation_mean)) + 
  geom_point(aes( size = nSharedOTUs), pch = 21, fill= "white") +
  labs(y = "Pairwise genetic distance", x = "Mean elevation") +
  geom_smooth(method = "lm", colour = "grey20", linetype = "dashed") +
  scale_radius(name = "Number of shared OTUs", range = c(1,6), breaks = c(40, 50, 60, 70, 80), labels = c(40, 50, 60, 70, 80), limits = c(40,80)) +
  theme(panel.background = element_blank())

zOTU_betsite_beta <- plot_grid(zOTU_betsite_bray + theme(legend.position = "none"),
                               zOTU_betsite_stddist+ theme(legend.position = "none"), labels = "auto")
zOTU_betsite_beta_wleg <- plot_grid(zOTU_betsite_beta, get_legend(zOTU_betsite_bray+theme(legend.position = "bottom")),
                                    nrow = 2, rel_heights = c(1,0.2))
ggsave(zOTU_betsite_beta_wleg, filename = file.path(fig.dir, "zOTU_betsite_beta.pdf"), width = 9, height = 5.5)

zOTU_beta_distance_stainback <- subset(zOTU_beta_distance, ((elevation_X1 <= 1000 & elevation_X2 > 1000) | (elevation_X2 <= 1000 & elevation_X1 > 1000)) &
                                         site_X1 == "Stainback" & site_X1 == site_X2)
zOTU_beta_distance_stainback$elevationCompare <- apply(zOTU_beta_distance_stainback[,c("elevation_X1", "elevation_X2")], FUN = max, MARGIN = 1)

zOTU_beta_distance_laupahoehoe <- subset(zOTU_beta_distance, ((elevation_X1 <= 1000 & elevation_X2 > 1000) | (elevation_X2 <= 1000 & elevation_X1 > 1000)) &
                                           site_X1 == "Laupahoehoe" & site_X1 == site_X2)
zOTU_beta_distance_laupahoehoe$elevationCompare <- apply(zOTU_beta_distance_laupahoehoe[,c("elevation_X1", "elevation_X2")], FUN = max, MARGIN = 1)

stainback_zOTU_bray <- ggplot(aes(x = elevationCompare, y = mean_haplo_bray), 
                              data = zOTU_beta_distance_stainback) +
  geom_point(aes(size = nSharedOTUs), pch = 21, fill= "white") +
  labs(y = "Pairwise Bray-Curtis haplotype dissimilarity", x = "Elevation") +
  geom_smooth(method = "lm", colour = "grey20") +
  scale_radius(name = "Number of shared OTUs", range = c(1,6), breaks = c(25, 37.5, 50, 62.5, 75), labels = c(25, 37.5, 50, 62.5, 75), limits = c(25,75))

stainback_zOTU_dist <- ggplot(aes(x = elevationCompare, y = mean_stddist), 
                              data = zOTU_beta_distance_stainback) +
  geom_point(aes(size = nSharedOTUs), pch = 21, fill= "white") +
  labs(y = "Pairwise genetic distance", x = "Elevation") +
  geom_smooth(method = "lm", colour = "grey20") +
  scale_radius(name = "Number of shared OTUs", range = c(1,6), breaks = c(25, 37.5, 50, 62.5, 75), labels = c(25, 37.5, 50, 62.5, 75), limits = c(25,75))

laupahoehoe_zOTU_bray <- ggplot(aes(x = elevationCompare, y = mean_haplo_bray), 
                                data = zOTU_beta_distance_laupahoehoe) +
  geom_point(aes(size = nSharedOTUs), pch = 21, fill= "white") +
  labs(y = "Pairwise Bray-Curtis haplotype dissimilarity", x = "Elevation") +
  geom_smooth(method = "lm", colour = "grey20") +
  scale_radius(name = "Number of shared OTUs", range = c(1,6), breaks = c(45, 55, 65, 75, 85), labels = c(45, 55, 65, 75, 85), limits = c(45,85))

laupahoehoe_zOTU_dist <- ggplot(aes(x = elevationCompare, y = mean_stddist), 
                                data = zOTU_beta_distance_laupahoehoe) +
  geom_point(aes(size = nSharedOTUs), pch = 21, fill= "white") +
  labs(y = "Pairwise genetic distance", x = "Elevation") +
  geom_smooth(method = "lm", colour = "grey20") +
  scale_radius(name = "Number of shared OTUs", range = c(1,6), breaks = c(45, 55, 65, 75, 85), labels = c(45, 55, 65, 75, 85), limits = c(45,85))

