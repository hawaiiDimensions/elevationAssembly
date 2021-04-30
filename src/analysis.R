## ECOLOGICAL ASSEMBLY OF THE HAWAIIAN ARTHROPOD COMMUNITIES
# Author: Jun Ying Lim
# Rarefies OTU tables for analysis

## PACKAGES ============
#rm(list = ls())
library(stringr); library(plyr); library(reshape2) # data manipulation tools
library(vegan) # calculating beta diversity
library(geosphere) # calculating geographic distances
library(ggplot2); library(ggrepel)

## IMPORT DATA ============
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

laupahoehoe_siteIDs <- subset(siteData, site.1 == "Laupahoehoe")$site.id
stainback_siteIDs <- subset(siteData, site.1 == "Stainback")$site.id

# Import otu data
otu_native_files <- list.files(file.path(data.dir, "OTU_native"))
OTU_mat_list <- lapply(otu_native_files, FUN = function(x){readRDS(file.path(data.dir, "OTU_native", x))})

# Import taxonmic reference
taxonData <- readRDS(file.path(data.dir, "taxonData.rds"))
OTUtaxonData <- taxonData[c("OTU_ID", "V4")]
OTUtaxonData <- OTUtaxonData[!duplicated(OTUtaxonData),]
OTUtaxonData$V4 <- as.vector(OTUtaxonData$V4)

# Remove sites lower than 800 m elevation
lowsiteids <- subset(siteData, elevation <= 800)$site.id
inclsiteids <- subset(siteData, elevation > 800)$site.id

# Import information on the number of specimens 
# NOTE: Analyses are based on just six taxonomic orders whereas the number of individuals also includes specimens from other orders and hence will not be accurate
specimenCounts <- read.delim(file.path(data.dir, "specimennumber.txt"), header = FALSE)
names(specimenCounts) <- c("Site_ID", "Transect", "Island", "Method", "Year", "SizeCategory", "NoOfIndividuals", "RarefiedReadAbund")

siteTotalCounts <- ddply(.data = subset(specimenCounts, Site_ID %in% inclsiteids),
                            .variables = .(Transect, Site_ID), summarise,
                            TotalNoIndividualsPerSite = sum(NoOfIndividuals))
specimenCounts2 <- merge(y = subset(specimenCounts, Site_ID %in% inclsiteids),
                         x = siteTotalCounts) 

specimenCounts2$PercentTotalIndividuals <- specimenCounts2$NoOfIndividuals / specimenCounts2$TotalNoIndividualsPerSite

percentSize <- ddply(.data = specimenCounts2,
                     .variables = .(Transect, SizeCategory),
                     summarise,
                     avgPercent = round(mean(PercentTotalIndividuals) * 100, 2))

# Average no. of individuals for each transect
avg_siteTransectsCounts <- ddply(.data = siteTotalCounts,
                                     .variables = .(Transect), summarise,
                                     AvgTotalNoIndividualsPerSite = mean(TotalNoIndividualsPerSite))
# Total no. of invidiuals for eahc transect
totalTransectCounts <- ddply(.data = subset(specimenCounts, Site_ID %in% inclsiteids),
                             .variables = .(Transect), summarise,
                             TotalNoIndividuals = sum(NoOfIndividuals))

## ALPHA DIVERSITY WITH ELEVATION  =========
# Combine all rarefied datasets into a single dataframe
OTU_mat_df <- melt(OTU_mat_list, varnames = c("OTU_ID", "site.id"),
                   value.name = "nReads")

# Merge OTU data to get higher taxonomic orders
OTU_mat_df_taxon <- merge(OTU_mat_df, OTUtaxonData[c("OTU_ID", "V4")], all.x = TRUE)

# Include only the six taxonomic orders
OTU_mat_df2 <- subset(OTU_mat_df_taxon, V4 %in% c("Araneae", "Coleoptera", "Hemiptera", "Lepidoptera", "Psocoptera", "Orthoptera"))

# Merge to get site information
OTU_mat_df3 <- merge(OTU_mat_df2, siteData[c("site.id", "site", "elevation")], by = "site.id")

## Create list of OTUs (all sites including lowland sites) 
allOTUs <- unique(as.vector(subset(OTU_mat_df3, nReads > 0)$OTU_ID))
laupOTUs <- unique(as.vector(subset(OTU_mat_df3, site == "Laupahoehoe" & nReads > 0)$OTU_ID))
stnbkOTUs <- unique(as.vector(subset(OTU_mat_df3, site == "Stainback" & nReads > 0)$OTU_ID))

length(laupOTUs) # 584 otus
length(stnbkOTUs) # 584 otus
length(allOTUs) # 697 OTUs
length(union(laupOTUs, stnbkOTUs)) # 697 OTUs

laupEndemicOTUs <- laupOTUs[!laupOTUs %in% stnbkOTUs] #113
stnbkEndemicOTUs <- stnbkOTUs[!stnbkOTUs %in% laupOTUs] # 113

laupEndemicOTUs %in% stnbkEndemicOTUs
stnbkEndemicOTUs %in% laupEndemicOTUs

calculateRichness <- function(x){
  elevation <- x$elevation[1]
  # Calculate number of unique OTUs as a proxy for species richness
  OTU_list <- unique(as.vector(x$OTU_ID))
  S <- length(OTU_list)
  if(x$site[1] == "Laupahoehoe"){
    S_end <- length(OTU_list[OTU_list %in% laupEndemicOTUs])   
  } 
  if(x$site[1] == "Stainback"){
    S_end <- length(OTU_list[OTU_list %in% stnbkEndemicOTUs])
  }

  S_araneae <- length(as.vector(unique(x$OTU_ID[x$V4 == "Araneae"])))
  S_coleoptera <- length(as.vector(unique(x$OTU_ID[x$V4 == "Coleoptera"])))
  S_hemiptera <-  length(as.vector(unique(x$OTU_ID[x$V4 == "Hemiptera"])))
  S_lepidoptera <-  length(as.vector(unique(x$OTU_ID[x$V4 == "Lepidoptera"])))
  S_psocoptera <-  length(as.vector(unique(x$OTU_ID[x$V4 == "Psocoptera"])))
  S_orthoptera <-  length(as.vector(unique(x$OTU_ID[x$V4 == "Orthoptera"])))
  
  pi_araneae <- sum(x$nReads[x$V4 == "Araneae"]) / sum(x$nReads)
  pi_coleoptera <- sum(x$nReads[x$V4 == "Coleoptera"]) / sum(x$nReads)
  pi_hemiptera <- sum(x$nReads[x$V4 == "Hemiptera"]) / sum(x$nReads)
  pi_lepidoptera <- sum(x$nReads[x$V4 == "Lepidoptera"]) / sum(x$nReads)
  pi_psocoptera <- sum(x$nReads[x$V4 == "Psocoptera"]) / sum(x$nReads)
  pi_orthoptera <- sum(x$nReads[x$V4 == "Orthoptera"]) / sum(x$nReads)
  
  pi <- x$nReads / sum(x$nReads)
  H <- -1 * sum( pi * log(pi))
  data.frame(S, S_araneae, S_coleoptera, S_hemiptera, S_lepidoptera, S_psocoptera, S_orthoptera,
             H, pi_araneae, pi_coleoptera, pi_hemiptera, pi_lepidoptera, pi_psocoptera, pi_orthoptera,
             S_end, elevation)
}

OTU_mat_nsp_site <- ddply(subset(OTU_mat_df3, nReads >0 ), .variables = .(L1, site, site.id),
                     .fun = calculateRichness,
                     .progress = "text")
saveRDS(OTU_mat_nsp_site, file.path(res.dir, "OTU_nsp_site.rds"))


# t-test of OTU richness 
OTU_mat_nsp_site2 <- ddply(.data = subset(OTU_mat_nsp_site, !site.id %in% lowsiteids),
                           .variables = .(site,site.id), summarise,
                           avgS = mean(S),
                           avgS_araneae = mean(S_araneae),
                           avgS_coleoptera = mean(S_coleoptera),
                           avgS_hemiptera = mean(S_hemiptera),
                           avgS_lepidoptera = mean(S_lepidoptera),
                           avgS_psocoptera = mean(S_psocoptera),
                           avgS_orthoptera = mean(S_orthoptera), 
                           avgH = mean(H))
t.test(avgS~site, data = OTU_mat_nsp_site2) # Laup = 108, Stainback = 83
t.test(avgS_araneae~site, data = OTU_mat_nsp_site2) # opp sign
t.test(avgS_coleoptera~site, data = OTU_mat_nsp_site2)
t.test(avgS_hemiptera~site, data = OTU_mat_nsp_site2)
t.test(avgS_orthoptera~site, data = OTU_mat_nsp_site2) # not sig
t.test(avgS_lepidoptera~site, data = OTU_mat_nsp_site2)
t.test(avgS_psocoptera~site, data = OTU_mat_nsp_site2)

## BETA DIVERSITY PATTERNS =========
# Calculate the average rarefied read abundance for all OTUs
OTU_mat_avg_df <- ddply(OTU_mat_df3, .variables = .(site.id, OTU_ID),
                        summarise, avgRead = mean(nReads), .progress = "text")
OTU_mat_avg <- reshape2::acast(OTU_mat_avg_df, formula = site.id~OTU_ID, value.var = "avgRead")

OTU_mat_avg_subset <- OTU_mat_avg[!rownames(OTU_mat_avg) %in% (lowsiteids),]
OTU_mat_avg_subset <- OTU_mat_avg_subset[,colSums(OTU_mat_avg_subset) > 0]

# Run NMDS for all orders togehter
site_nmds_bray <- metaMDS(OTU_mat_avg_subset, dist = "bray", k = 2)
site_nmds_bray_df <- data.frame(site_nmds_bray$points)
site_nmds_bray_df$site.id <- rownames(site_nmds_bray$points)
site_nmds_bray_df2 <- merge(site_nmds_bray_df, siteData)
write.csv(site_nmds_bray_df2, file.path(res.dir, "site_bray_nmds.csv"), row.names = F)

# Run NMDS with presence-absence
site_nmds_bray_pa <- metaMDS(decostand(OTU_mat_avg_subset, method = "pa"), dist = "bray", k = 2)
site_nmds_bray_pa_df <- data.frame(site_nmds_bray_pa$points)
site_nmds_bray_pa_df$site.id <- rownames(site_nmds_bray_pa$points)
site_nmds_bray_pa_df2 <- merge(site_nmds_bray_pa_df, siteData)
write.csv(site_nmds_bray_pa_df2, file.path(res.dir, "site_bray_pa_nmds.csv"), row.names = F)

# Calculate beta-diversity over elevation
# RUn the commented lines if you want to exclude invasives
# invasive_class <- read.delim(file.path(data.dir, "AllInvasiveTaxa.txt"), header = FALSE)
# invasive_otus <- unique(invasive_class$V3)
# sum(colnames(OTU_mat_avg_subset) %in% invasive_otus) # 42 invasive OTUs
# OTU_mat_avg_subset <- OTU_mat_avg_subset[,!colnames(OTU_mat_avg_subset) %in% invasive_otus]

elev_bray_all <- as.matrix(vegdist(OTU_mat_avg_subset, method = "bray", diag = T, upper = T))
elev_braypa_all <- as.matrix(vegdist(OTU_mat_avg_subset, method = "bray", diag = T, upper = T, binary = T))

elev_bray_araneae <- as.matrix(vegdist(OTU_mat_avg_subset[,colnames(OTU_mat_avg_subset) %in% subset(taxonData, V4 == "Araneae")$OTU_ID], method = "bray", diag = T, upper = T))
elev_bray_hemiptera <- as.matrix(vegdist(OTU_mat_avg_subset[,colnames(OTU_mat_avg_subset) %in% subset(taxonData, V4 == "Hemiptera")$OTU_ID], method = "bray", diag = T, upper = T))
elev_bray_lepidoptera <- as.matrix(vegdist(OTU_mat_avg_subset[,colnames(OTU_mat_avg_subset) %in% subset(taxonData, V4 == "Lepidoptera")$OTU_ID], method = "bray", diag = T, upper = T))
elev_bray_coleoptera <- as.matrix(vegdist(OTU_mat_avg_subset[,colnames(OTU_mat_avg_subset) %in% subset(taxonData, V4 == "Coleoptera")$OTU_ID], method = "bray", diag = T, upper = T))
elev_bray_orthoptera <- as.matrix(vegdist(OTU_mat_avg_subset[,colnames(OTU_mat_avg_subset) %in% subset(taxonData, V4 == "Orthoptera")$OTU_ID], method = "bray", diag = T, upper = T))
elev_bray_psocoptera <- as.matrix(vegdist(OTU_mat_avg_subset[,colnames(OTU_mat_avg_subset) %in% subset(taxonData, V4 == "Psocoptera")$OTU_ID], method = "bray", diag = T, upper = T))

elev_braypa_araneae <- as.matrix(vegdist(OTU_mat_avg_subset[,colnames(OTU_mat_avg_subset) %in% subset(taxonData, V4 == "Araneae")$OTU_ID], method = "bray", diag = T, upper = T, binary = T))
elev_braypa_hemiptera <- as.matrix(vegdist(OTU_mat_avg_subset[,colnames(OTU_mat_avg_subset) %in% subset(taxonData, V4 == "Hemiptera")$OTU_ID], method = "bray", diag = T, upper = T, binary = T))
elev_braypa_lepidoptera <- as.matrix(vegdist(OTU_mat_avg_subset[,colnames(OTU_mat_avg_subset) %in% subset(taxonData, V4 == "Lepidoptera")$OTU_ID], method = "bray", diag = T, upper = T, binary = T))
elev_braypa_coleoptera <- as.matrix(vegdist(OTU_mat_avg_subset[,colnames(OTU_mat_avg_subset) %in% subset(taxonData, V4 == "Coleoptera")$OTU_ID], method = "bray", diag = T, upper = T, binary = T))
elev_braypa_orthoptera <- as.matrix(vegdist(OTU_mat_avg_subset[,colnames(OTU_mat_avg_subset) %in% subset(taxonData, V4 == "Orthoptera")$OTU_ID], method = "bray", diag = T, upper = T, binary = T))
elev_braypa_psocoptera <- as.matrix(vegdist(OTU_mat_avg_subset[,colnames(OTU_mat_avg_subset) %in% subset(taxonData, V4 == "Psocoptera")$OTU_ID], method = "bray", diag = T, upper = T, binary = T))

# library(betapart)
# elev_betapart <- beta.pair(as.data.frame(decostand(OTU_mat_avg_subset, method = "pa")))
# elev_sorenturnover <- as.matrix(elev_betapart$beta.sne)
# elev_sorennested <- as.matrix(elev_betapart$beta.sim)

xy <- t(combn(colnames(elev_bray_all), 2)) # this removes backward comparisons

distDF <- data.frame(xy,
                     #elev_sorenturnover = elev_sorenturnover[xy],
                     #elev_sorennested = elev_sorennested[xy],
                     bc_all= elev_bray_all[xy], 
                     bcpa_all = elev_braypa_all[xy],
                     bc_araneae = elev_bray_araneae[xy],
                     bc_hemiptera =  elev_bray_hemiptera[xy],
                     bc_lepidoptera = elev_bray_lepidoptera[xy],
                     bc_coleoptera = elev_bray_coleoptera[xy],
                     bc_orthoptera = elev_bray_orthoptera[xy],
                     bc_psocoptera = elev_bray_psocoptera[xy],
                     bcpa_araneae = elev_braypa_araneae[xy],
                     bcpa_hemiptera =  elev_braypa_hemiptera[xy],
                     bcpa_lepidoptera = elev_braypa_lepidoptera[xy],
                     bcpa_coleoptera = elev_braypa_coleoptera[xy],
                     bcpa_orthoptera = elev_braypa_orthoptera[xy],
                     bcpa_psocoptera = elev_braypa_psocoptera[xy])
distDF2 <- merge(x = distDF, y = siteData[, c("site.id", "elevation", "site")], by.x = "X1", by.y = "site.id")
distDF3 <- merge(x = distDF2, y = siteData[, c("site.id", "elevation", "site")], by.x = "X2", by.y = "site.id", suffixes = c("_X1", "_X2"))
write.csv(distDF3, file = file.path(res.dir, "elevation_dist.csv"), row.names = F)

# Average dissimilarity within sites and between sites
sum(!laupahoehoe_siteIDs %in% lowsiteids) # 23 sites
sum(!stainback_siteIDs %in% lowsiteids) # 36 sites
(59*58) / 2 == nrow(distDF3) # 828 total comparisons
23 * 22 / 2 # number of laupahoehoe comparisons = 253
36*35 /2 # number of stainback comparisons = 630
mean(subset(distDF3, site_X1 == "Laupahoehoe" & site_X2 == "Laupahoehoe")$bc_all) # 0.64
length(subset(distDF3, site_X1 == "Laupahoehoe" & site_X2 == "Laupahoehoe")$bc_all) # 253 pw comparisons
range(subset(distDF3, site_X1 == "Laupahoehoe" & site_X2 == "Laupahoehoe")$bc_all)
mean(subset(distDF3, site_X1 == "Stainback" & site_X2 == "Stainback")$bc_all) # 0.67
length(subset(distDF3, site_X1 == "Stainback" & site_X2 == "Stainback")$bc_all) # 630 comparisons
range(subset(distDF3, site_X1 == "Stainback" & site_X2 == "Stainback")$bc_all)
mean(subset(distDF3, site_X1 !=  site_X2)$bc_all) # 0.80
range(subset(distDF3, site_X1 !=  site_X2)$bc_all) # 0.80
nrow(subset(distDF3, site_X1 !=  site_X2)) # 828 comparisons

# Only compare samples that are within 100m of each other
distDF3_subset <- subset(distDF3, abs(distDF3$elevation_X1 - distDF3$elevation_X2) <= 100)
distDF3_subset$avgElevation <- rowMeans(distDF3_subset[,c("elevation_X1", "elevation_X2")])
distDF3_subset_compare <- subset(distDF3_subset, site_X1 != site_X2)

test.labels <- c("All Orders", "Araneae", "Hemiptera",
                 "Lepidoptera", "Coleoptera", "Orthoptera", 
                 "Psocoptera", "All Orders (PA)", "Araneae (PA)", "Hemiptera (PA)",
                 "Lepidoptera (PA)", "Coleoptera (PA)", "Orthoptera (PA)", 
                 "Psocoptera (PA)")
test.col <- c("bc_all", "bc_araneae", "bc_hemiptera", "bc_lepidoptera", "bc_coleoptera", "bc_orthoptera", "bc_psocoptera",
              "bcpa_all", "bcpa_araneae", "bcpa_hemiptera", "bcpa_lepidoptera", "bcpa_coleoptera", "bcpa_orthoptera", "bcpa_psocoptera")

res <- list()
for(i in 1:length(test.col)){
  cor.res <- cor.test(x = distDF3_subset_compare[[test.col[i]]],
                      y = distDF3_subset_compare[["avgElevation"]],
                      method = "spearman", exact = FALSE)
  
  res[[i]] <- data.frame("rho" = round(as.vector(cor.res$estimate), 2),
                         "p.value" = round(as.vector(cor.res$p.value), 2),
                         "model" = test.labels[i] )
  
  
}
bray_res <- do.call("rbind", res)

# How correlated are the abundance-weighted vs. presence-absence
cor.test(distDF3_subset_compare$bc_all, distDF3_subset_compare$bcpa_all)
cor.test(distDF3_subset_compare$bc_araneae, distDF3_subset_compare$bcpa_araneae) # marginally not sig.
cor.test(distDF3_subset_compare$bc_hemiptera, distDF3_subset_compare$bcpa_hemiptera)
cor.test(distDF3_subset_compare$bc_coleoptera, distDF3_subset_compare$bcpa_coleoptera)
cor.test(distDF3_subset_compare$bc_orthoptera, distDF3_subset_compare$bcpa_orthoptera)
cor.test(distDF3_subset_compare$bc_psocoptera, distDF3_subset_compare$bcpa_psocoptera) # not correlated

## ZOTU BETA DIVERSITY PATTERNS =============
zotu_native_files <- list.files(file.path(data.dir, "zOTU_native"))
zOTU_mat_list <- lapply(zotu_native_files, FUN = function(x){readRDS(file.path(data.dir, "zOTU_native", x))})

# Combine list into a single dataframe
zOTU_mat_df <- melt(zOTU_mat_list, varnames = c("zOTU_ID", "site.id"),
                    value.name = "nReads")
zOTU_mat_df2 <- zOTU_mat_df[!zOTU_mat_df$site.id %in% (lowsiteids),]

# Merge with taxonomic information
zOTU_mat_df_taxon <- merge(zOTU_mat_df2, taxonData[c("OTU_ID", "zOTU_ID", "V4")], by = "zOTU_ID")

# Only include the six focal orders
zOTU_mat_df3 <- subset(zOTU_mat_df_taxon, V4 %in% c("Araneae", "Coleoptera", "Hemiptera", "Lepidoptera", "Psocoptera", "Orthoptera"))

# Calculate mean read abundance across rarefied datasets
zOTU_avg_df <- ddply(subset(zOTU_mat_df3, nReads > 0), .variables = .(site.id, zOTU_ID),
                     summarise,
                     avgRead = mean(nReads),
                     V4 = V4[1], 
                     OTU_ID = OTU_ID[1],
                     .progress = "text")

# Only include OTUs with more than 1 zOTU (this does not change results but overall dissimilarity will be higher)
# zOTU_taxon_df_subset <- ddply(zOTU_avg_df,
#                               .variables = .(OTU_ID),
#                               .fun = function(x){ if(length(unique(x$zOTU_ID)) > 1){ return(x) } else {return(NULL)}})

# Only include OTUs with more than 1 site
zOTU_taxon_df_subset <- ddply(zOTU_avg_df, .variables= .(OTU_ID),
                              .fun = function(x){
                                if(length(unique(as.vector(x$site.id))) > 1){
                                  return(x)
                                } else {
                                  return(NULL)
                                }} )

# Calculate bray curtis
zOTU_beta <- ddply(.data = zOTU_taxon_df_subset,
                   .variables = .(OTU_ID),
                   .fun = function(x){
                     #x = subset(zOTU_taxon_df_subset, OTU_ID == "OTU100")
                     #x = subset(zOTU_taxon_df_subset, OTU_ID == "OTU1003") #this OTU only has one zOTU
                     OTU_mat <- reshape2::acast(x, formula = site.id~zOTU_ID, value.var = "avgRead", fill = 0)
                     OTU_mat2 <- OTU_mat[rowSums(OTU_mat) > 0,,drop = FALSE]
                     
                     distDF <- expand.grid(rownames(OTU_mat2), rownames(OTU_mat2), stringsAsFactors = FALSE)
                     beta_mat <- as.matrix(vegdist(OTU_mat2, method = "bray", diag = T, upper = T))
                     beta_pamat <- as.matrix(vegdist(OTU_mat2, method = "bray", diag = T, upper = T, binary = T))
                     distDF$bray <- as.vector(beta_mat)
                     distDF$braypa <- as.vector(beta_pamat)
                     #distDF <- melt(beta_mat, value.name = "haplo_bray")
                     return(distDF)
                   })

# Remove the reverse pairwise comparisons
zOTU_beta2 <- ddply(zOTU_beta, .variables = .(OTU_ID), .fun = function(x){
  x_sort = t(apply(x[,c("Var1", "Var2")], 1, sort))  
  x_unique <- x[!duplicated(x_sort), ]
  x_unique <- subset(x_unique, !Var1 == Var2)
  return(x_unique)
})
zOTU_dist <- merge(x = zOTU_beta2, y = siteData[, c("site.id", "elevation", "site")], by.x = "Var1", by.y = "site.id")
zOTU_dist2 <- merge(x = zOTU_dist, y = siteData[, c("site.id", "elevation", "site")], by.x = "Var2", by.y = "site.id", suffixes = c("_X1", "_X2"))

zOTU_dist_subset <- subset(zOTU_dist2, abs(elevation_X1 - elevation_X2) <= 100 & site_X1 != site_X2)
zOTU_dist_subset$avgElevation <- apply(zOTU_dist_subset[,c("elevation_X1", "elevation_X2")], MARGIN = 1, FUN = mean)

dim(zOTU_dist_subset) # 9021 pairwise OTU comparisons
length(unique(zOTU_dist_subset$OTU_ID)) # 318 unique OTUs
nrow(zOTU_dist_subset[!duplicated(zOTU_dist_subset[,c("Var1", "Var2")]),]) # 148 site comparisons

# Calculate average bray-curtis of OTUs shared between sites
zOTU_beta_avg <- ddply(zOTU_dist_subset, .variables = .(Var2,Var1),
                       .fun = function(x){
                         #x<- subset(zOTU_dist_subset, Var2 == 521 & Var1 == 693)
                         nSharedOTUs = length(x$OTU_ID) # number of OTUs compared between sites
                         #mean_gendist = mean(x$avgdist, na.rm = T)
                         #mean_stddist = mean(x$stddist, na.rm = T)
                         mean_bray = mean(x$bray, na.rm = T)
                         mean_bray_pa = mean(x$braypa, na.rm = T)
                         elevation_X1 = x$elevation_X1[1]
                         elevation_X2 = x$elevation_X2[1]
                         site_X1 = x$site_X1[1]
                         site_X2 = x$site_X2[1]
                         avgElevation = x$avgElevation[1]
                         data.frame(nSharedOTUs, mean_bray, mean_bray_pa,
                                    elevation_X1, elevation_X2, site_X1, site_X2,
                                    avgElevation)
                       })
saveRDS(zOTU_beta_avg, file.path(res.dir, "zOTU_beta_avg.rds"))

cor.test(x = zOTU_beta_avg$mean_bray, y =  zOTU_beta_avg$avgElevation, method = "spearman", exact = FALSE) # 0.548
cor.test(x = zOTU_beta_avg$mean_bray_pa, y =  zOTU_beta_avg$avgElevation, method = "spearman", exact = FALSE) # 0.759

# ggplot(data = subset(zOTU_dist_subset, OTU_ID %in% test_set)) + 
#   geom_point(aes(y = braypa, x = avgElevation, colour = OTU_ID)) +
#   geom_smooth(aes(y = braypa, x = avgElevation, colour = OTU_ID), se = FALSE, method = "lm") +
#   theme(legend.position = "none")
# 
# library(lme4)
# elev_mod1 <- lm(bray ~ avgElevation, data = zOTU_dist_subset)
# AIC(elev_mod1)
# elev_mod2 <- lmer(bray ~ avgElevation + (1|OTU_ID), data = zOTU_dist_subset)
# AIC(elev_mod2)
# summary(elev_mod2)
# 
# test_set <- names(table(zOTU_dist_subset$OTU_ID))[table(zOTU_dist_subset$OTU_ID) > 20] # 112 OTUs with more than 10 comparisons
# elev_mod3 <- lmer(bray ~ avgElevation + (1|OTU_ID) + (avgElevation|OTU_ID), data = subset(zOTU_dist_subset, OTU_ID %in% test_set))
# AIC(elev_mod3)

## INVASIVE SPECIES =============
#testdata <- subset(OTU_mat_df3, nReads >0 & V4 == "Araneae" & L1 == 1)
#spider_inv_otu <- ara_invasives_otu$V2[ara_invasives_otu$status == "invasive"]
#OTU_mat_avg_subset <- OTU_mat_avg_subset[,!colnames(OTU_mat_avg_subset) %in% spider_inv_otu]
ara_invasives <- read.table(file.path(data.dir,"InvasiveAraneae.txt"))
length(unique(ara_invasives$V2))

subset(invasive_class, V3 == "OTU765")
subset(taxonData, OTU_ID =="OTU765")

invasive_class[invasive_class$V9 == "Trombidiformes",]
OTU_mat_avg_df

length(unique(taxonData$OTU_ID))
head(invasive_class)
length(unique(invasive_class$V3))
head

invasive_otu_class <- ddply(.data = invasive_class,
                           .variables = .(V2),
                           summarise,
                           status =  names(which.max(table(V13))))

testdata2 <- merge(x= testdata, y = ara_invasives_otu, by.x ="OTU_ID", by.y = "V2", all.x = TRUE)
testdata3 <- ddply(testdata2, .variables = .(site.id), .fun =  function(x){
  temp <- tapply(INDEX = x$status, X = x$nReads, FUN = sum)
  temp2 <- temp/ sum(x$nReads)
  data.frame(status = names(temp), prop = temp2)
})

invasive_plot <- ggplot(data = testdata3) + geom_bar(aes(fill = factor(status), x = factor(site.id), weight = prop)) + theme(axis.text = element_text(angle = 90))
ggsave("~/Dropbox/projects/2017/hawaiiCommunityAssembly/elevationAssemblyHawaii/manuscript/invasive_plot.pdf", invasive_plot, width = 10, height = 4)





## NICHE CONSERVATISM =============
siteData2 <- subset(siteData, !site.id %in% lowsiteids)
site_temprange <- tapply(siteData2$t_ann, siteData2$site, range)
site_pptrange <- tapply(siteData2$rf_ann, siteData2$site, range)

t_res_list_t10 <- list()
rf_res_list_t10 <- list()
t_res_list_t5 <- list()
rf_res_list_t5 <- list()

# Calculate niche position for all OTUs for all rarefied datasets
for(i in 1:length(OTU_mat_list)){
  # Counter
  pb = txtProgressBar(min = 0, max = length(OTU_mat_list), style = 3, initial = 0) 
  setTxtProgressBar(pb, i)

  OTU_mat <- OTU_mat_list[[i]]
  
  # Remove sites below 800 metres
  OTU_mat <- OTU_mat[,!colnames(OTU_mat) %in% lowsiteids]
  
  # Only analyze OTUs that have been matched to certain orders
  targetOrders <- 
  #c("Acari", "Araneae", "Coleoptera", "Lepidoptera", "Orthoptera", "Neuroptera", "Psocoptera", "Hemiptera")
  otu_target <- subset(OTUtaxonData, V4 %in% c("Araneae", "Hemiptera", "Lepidoptera", "Psocoptera", "Coleoptera", "Orthoptera"))$OTU_ID
  
  # Convert OTU matrix into long table format
  OTU_df <- melt(data = OTU_mat[rownames(OTU_mat) %in% otu_target,], value.name = "nReads", varnames = c("OTU_ID", "site.id"))
  
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
  
  # Flag OTUs whose range boundaries are outside the sampling ranges for each site
  # ( 1 = in bounds, 0 = out of bounds; since upper and lower bounds are 0.75 and 0.25;
  # not more than a quarter of occurrences may be at the most extreme site)
  flagLimit <- function(x){
    x$rf_bounds <- with(data = x,
                        ifelse(site == "Laupahoehoe" & rf_upper >= site_pptrange$Laupahoehoe[2] |
                            site == "Laupahoehoe" & rf_lower <= site_pptrange$Laupahoehoe[1] | 
                            site == "Stainback" & rf_upper >= site_pptrange$Stainback[2] |
                            site == "Stainback" & rf_lower <= site_pptrange$Stainback[1], 0, 1))
    x$temp_bounds <- with(data = x,
                          ifelse(site == "Laupahoehoe" & temp_upper >= site_temprange$Laupahoehoe[2] |
                                   site == "Laupahoehoe" &  temp_lower <= site_temprange$Laupahoehoe[1] |
                                   site == "Stainback" & temp_upper >= site_temprange$Stainback[2] |
                                   site == "Stainback" & temp_lower <= site_temprange$Stainback[1], 0, 1))
    return(x)
  }
  OTU_clim_site_10 <- flagLimit(x = OTU_clim_site_10)
  OTU_clim_site_5 <- flagLimit(x = OTU_clim_site_5)
  
  t_res_list_t10[[i]] <- subset(OTU_clim_site_10, temp_bounds == 1)
  rf_res_list_t10[[i]] <- subset(OTU_clim_site_10, rf_bounds == 1)
  t_res_list_t5[[i]] <- subset(OTU_clim_site_5, temp_bounds == 1)
  rf_res_list_t5[[i]] <- subset(OTU_clim_site_5, rf_bounds == 1)
  
}
close(pb)

saveRDS(t_res_list_t10, file.path(res.dir, "t_res_list_t10.rds"))
saveRDS(rf_res_list_t10, file.path(res.dir, "rf_res_list_t10.rds"))
saveRDS(t_res_list_t5, file.path(res.dir, "t_res_list_t5.rds"))
saveRDS(rf_res_list_t5, file.path(res.dir, "rf_res_list_t5.rds"))

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

# Clean up data
nicheConsPrep <- function(x, value, otu_tab){
  site <- dcast(data = x, formula = OTU_ID~site, value.var = value)
  site_OTU <- merge(site, otu_tab)
  return(site_OTU)
}


nicheConsTest <- function(x){
  cor_test_res <- cor.test(x = x$Laupahoehoe,
                           y = x$Stainback,
                           method = "spearman",
                           exact = FALSE)
  t <- cor_test_res$statistic
  cor <- cor_test_res$estimate
  p <- cor_test_res$p.value
  n <- nrow(x)
  data.frame(t, cor, p, n)
}

t_med_t5  <- llply(.data = t_res_list_t5_subset, .fun = nicheConsPrep, value = "temp_median", otu_tab = OTUtaxonData)
rf_med_t5  <- llply(.data = rf_res_list_t5_subset, .fun = nicheConsPrep, value = "rf_median", otu_tab = OTUtaxonData)

saveRDS(t_med_t5, file.path(res.dir, "t_med_t5.rds"))
saveRDS(rf_med_t5, file.path(res.dir, "rf_med_t5.rds"))

t_med_corr_t5 <- ldply(.data = t_med_t5, .fun = nicheConsTest)
rf_med_corr_t5 <- ldply(.data = rf_med_t5, .fun = nicheConsTest)

write.csv(t_med_corr_t5, file.path(res.dir, "t_med_corr_t5.csv"), row.names = F)
write.csv(rf_med_corr_t5, file.path(res.dir, "rf_med_corr_t5.csv"), row.names = F)

round(mean(t_med_corr_t5$cor), 2)
round(range(t_med_corr_t5$cor), 2)
sum(t_med_corr_t5$p < 0.05) / 100
round(range(t_med_corr_t5$n),2 )

round(mean(rf_med_corr_t5$cor), 2)
round(range(rf_med_corr_t5$cor), 2)
sum(rf_med_corr_t5$p < 0.05) / 100
round(range(rf_med_corr_t5$n), 1)

# Re-run with a higher site threshold
t_med_t10  <- llply(.data = t_res_list_t10_subset, .fun = nicheConsPrep, value = "temp_median", otu_tab = OTUtaxonData)
rf_med_t10  <- llply(.data = rf_res_list_t10_subset, .fun = nicheConsPrep, value = "rf_median", otu_tab = OTUtaxonData)

saveRDS(t_med_t10, file.path(res.dir, "t_med_t10.rds"))
saveRDS(rf_med_t10, file.path(res.dir, "rf_med_t10.rds"))

t_med_corr_t10 <- ldply(.data = t_med_t10, .fun = nicheConsTest)
rf_med_corr_t10 <- ldply(.data = rf_med_t10, .fun = nicheConsTest)

write.csv(t_med_corr_t10, file.path(res.dir, "t_med_corr_t10.csv"), row.names = F)
write.csv(rf_med_corr_t10, file.path(res.dir, "rf_med_corr_t10.csv"), row.names = F)

round(mean(t_med_corr_t10$cor), 2)
round(range(t_med_corr_t10$cor), 2)
sum(t_med_corr_t10$p < 0.05) / 100
round(range(t_med_corr_t10$n),2 )

round(mean(rf_med_corr_t10$cor), 2)
round(range(rf_med_corr_t10$cor), 2)
sum(rf_med_corr_t10$p < 0.05) / 100
round(range(rf_med_corr_t10$n), 1)

# Re-run analysis with OTU that occur in ALL rarefied datasets
t_med_t5_df <- do.call("rbind",t_med_t5)
rf_med_t5_df <- do.call("rbind",rf_med_t5)
t_med_t5_counts <- table(t_med_t5_df$OTU_ID)
rf_med_t5_counts <- table(rf_med_t5_df$OTU_ID)
t_med_t5_subset <- lapply(t_med_t5, FUN = function(x) { subset(x, OTU_ID %in% names(t_med_t5_counts)[t_med_t5_counts ==100])})
rf_med_t5_subset <- lapply(rf_med_t5, FUN = function(x) { subset(x, OTU_ID %in% names(rf_med_t5_counts)[rf_med_t5_counts ==100]) }) 

t_med_t10_df <- do.call("rbind",t_med_t10)
rf_med_t10_df <- do.call("rbind",rf_med_t10)
t_med_t10_counts <- table(t_med_t10_df$OTU_ID)
rf_med_t10_counts <- table(rf_med_t10_df$OTU_ID)
t_med_t10_subset <- lapply(t_med_t10, FUN = function(x) { subset(x, OTU_ID %in% names(t_med_t10_counts)[t_med_t10_counts ==100])   })
rf_med_t10_subset <- lapply(rf_med_t10, FUN = function(x) { subset(x, OTU_ID %in% names(rf_med_t10_counts)[rf_med_t10_counts ==100])   })


t_med_corr_t5_subset <- ldply(.data = t_med_t5_subset, .fun = nicheConsTest)
rf_med_corr_t5_subset <- ldply(.data = rf_med_t5_subset, .fun = nicheConsTest)

t_med_corr_t10_subset <- ldply(.data = t_med_t10_subset, .fun = nicheConsTest)
rf_med_corr_t10_subset <- ldply(.data = rf_med_t10_subset, .fun = nicheConsTest)

write.csv(t_med_corr_t5_subset, file.path(res.dir, "t_med_corr_t5_subset.csv"), row.names = FALSE)
write.csv(rf_med_corr_t5_subset, file.path(res.dir, "rf_med_corr_t5_subset.csv"), row.names = FALSE)

write.csv(t_med_corr_t10_subset, file.path(res.dir, "t_med_corr_t10_subset.csv"), row.names = FALSE)
write.csv(rf_med_corr_t10_subset, file.path(res.dir, "rf_med_corr_t10_subset.csv"), row.names = FALSE)



round(mean(t_med_corr_t5_subset$cor), 2)
round(range(t_med_corr_t5_subset$cor), 2)
sum(t_med_corr_t5_subset$p < 0.05) / 100
round(range(t_med_corr_t5_subset$n),2 )

round(mean(rf_med_corr_t5_subset$cor), 2)
round(range(rf_med_corr_t5_subset$cor), 2)
sum(rf_med_corr_t5_subset$p < 0.05) / 100
round(range(rf_med_corr_t5_subset$n), 1)
