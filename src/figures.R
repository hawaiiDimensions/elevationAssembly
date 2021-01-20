## Generate figures
# Authors: Jun Ying Lim

## Packages ====================
#rm(list = ls())
library(readxl)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(viridis)
library(wesanderson)
library(ggnewscale)
library(plyr)
library(rgdal)
library(raster)
library(RStoolbox)
library(Hmisc) # capitalize function

berkcol3 = c("#003262", "#3B7EA1", "#9BBEA9", "#00B0DA", "#00A598", brewer.pal(9, "BuGn")[8], "#CFDD45",
             "#859438", "#FDB515", brewer.pal(9, "YlOrRd")[5], "#ED4E33", "#C4820E", "#D9661F",
             "#6C3302")

## Directories ====================
main.dir <- "~/Dropbox/Projects/2017/hawaiiCommunityAssembly/elevationAssemblyHawaii/"
res.dir <- file.path(main.dir, "results")
fig.dir <- file.path(main.dir, "figures")
data.dir <- file.path(main.dir, "data")
siteData <- read.csv(file.path(data.dir, "clim.final.csv"), stringsAsFactors = FALSE)
siteData <- subset(siteData, island == "BigIsland")

## CLIMATIC GRADIENTS =========
 rf_site_plot <- ggplot(data = subset(siteData,
                                     site %in%  c("Laupahoehoe", "Stainback") &
                                     elevation > 800 )) +
  geom_point(aes(y = rf_ann, x = elevation, fill = site ), pch = 21, size = 3, colour = "grey20" ) +
  scale_fill_manual(values = c("#003262", "#FDB515")) +
  labs(y = "Total annual rainfall (mm)", x = "Elevation (m)") +
  scale_x_continuous(breaks = c(800, 1000, 1200, 1400, 1600)) +
  guides(fill = guide_legend(title = NULL, override.aes = list(size = 5))) +
  theme(panel.background = element_blank(),
        legend.position = c(1, 1),
        #legend.background = element_rect(colour = "grey20"),
        legend.justification = c(1,1),
        legend.text = element_text(size = 10))

temp_site_plot <- ggplot(data = subset(siteData, 
                                       site %in%  c("Laupahoehoe", "Stainback") &
                                       elevation > 800)) +
  geom_point(aes(y = t_ann, x = elevation, fill = site ), pch = 21, size = 3, colour = "grey20" ) +
  scale_fill_manual(values = c("#003262", "#FDB515")) +
  labs(y = expression("Mean annual temperature ("~degree*C~")" ), x = "Elevation (m)") +
  scale_x_continuous(breaks = c(800, 1000, 1200, 1400, 1600)) +
  guides(fill = guide_legend(title = NULL, override.aes = list(size = 5))) +
  theme(panel.background = element_blank(),
        legend.position = c(1, 1),
        #legend.background = element_rect(colour = "grey20"),
        legend.justification = c(1,1),
        legend.text = element_text(size = 10))

transect_climate_plot <- plot_grid(rf_site_plot, temp_site_plot, labels = "auto")
ggsave(transect_climate_plot, filename = file.path(fig.dir, "transect_climate.pdf"), height = 4, width = 8)

## ALPHA DIVERSITY ================
OTU_mat_nsp_site <- readRDS(file.path(res.dir, "OTU_nsp_site.rds"))
OTU_mat_nsp_site$propE  <- OTU_mat_nsp_site$S_end / OTU_mat_nsp_site$S * 100

OTU_mat_nsp_site_avg <- ddply(OTU_mat_nsp_site, .variables = .(site.id), summarise,
                              avgH = mean(H),
                              avgS = mean(S),
                              avgE = mean(S_end),
                              avgPE = mean(propE),
                              elevation = elevation[1],
                              site = site[1],
                              minS = min(S),
                              maxS = max(S),
                              minH = min(H),
                              maxH = max(H),
                              avgS_araneae = mean(S_araneae),
                              avgS_coleoptera = mean(S_coleoptera),
                              avgS_hemiptera = mean(S_hemiptera),
                              avgS_lepidoptera = mean(S_lepidoptera),
                              avgS_psocoptera = mean(S_psocoptera),
                              avgS_orthoptera = mean(S_orthoptera))

# Plot OTU richness with elevation (Fig 2b) ==========
nsp_plot <- ggplot(data = subset(OTU_mat_nsp_site_avg, elevation > 800)) +
  geom_errorbar(aes(ymin = minS, ymax = maxS, x = elevation, colour = site), width = 0) +
  geom_point(aes(y = avgS, x = elevation, fill = site), data = subset(OTU_mat_nsp_site_avg, elevation > 800), colour = "black",
             pch = 21, size = 3) +
  geom_smooth(aes(y = avgS, x = elevation, group = site, colour = site),
              se = TRUE, method = "loess", alpha = 0.2, size = 1) +
  guides(colour = FALSE, fill = guide_legend(title = NULL, override.aes = list(size = 5))) +
  scale_fill_manual(values = c("#003262", "#FDB515")) +
  scale_colour_manual(values = c("#003262", "#FDB515")) +
  scale_x_continuous(breaks = c(800,1000,1200,1400, 1600)) +
  labs(x = "Elevation (m)",  y= "No. of OTUs") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.background = element_blank(),
        legend.position = c(1,1),
        legend.key = element_blank(),
        legend.justification = c(1,1))
ggsave(nsp_plot, filename = file.path(fig.dir, "nsp.pdf"), height = 5, width = 5)

# Plot OTU richness with elevation (Fig S2) ==========
H_plot <- ggplot(data = subset(OTU_mat_nsp_site_avg, elevation > 800)) +
  geom_errorbar(aes(ymin = minH, ymax = maxH, x = elevation, colour = site), width = 0) +
  geom_point(aes(y = avgH, x = elevation, fill = site), data = subset(OTU_mat_nsp_site_avg, elevation > 800), colour = "black",
             pch = 21, size = 3) +
  geom_smooth(aes(y = avgH, x = elevation, group = site, colour = site),
              se = TRUE, method = "loess", alpha = 0.2, size = 1) +
  guides(colour = FALSE, fill = guide_legend(title = NULL, override.aes = list(size = 5))) +
  scale_fill_manual(values = c("#003262", "#FDB515")) +
  scale_colour_manual(values = c("#003262", "#FDB515")) +
  scale_x_continuous(breaks = c(800,1000,1200,1400, 1600)) +
  labs(x = "Elevation (m)",  y= "Shannon's Index") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.background = element_blank(),
        legend.position = c(1,1),
        legend.key = element_blank(),
        legend.justification = c(1,1))
ggsave(H_plot, filename = file.path(fig.dir, "H.pdf"), height = 5, width = 5)

## ORDER-LEVEL ALPHA DIVERSITY  ================
OTU_nsp_order_melt <- melt(OTU_mat_nsp_site_avg, 
                            id.vars = c("site.id", "site", "elevation"),
                            measure.vars = c("avgS_araneae",
                                             "avgS_coleoptera",
                                             "avgS_hemiptera",
                                             "avgS_lepidoptera",
                                             "avgS_orthoptera",
                                             "avgS_psocoptera"),
                           value.name = "avgS")

levels(OTU_nsp_order_melt$variable) <- capitalize(gsub(levels(OTU_nsp_order_melt$variable), pattern = "avgS_", replacement = ""))

nsp_order_plot <- ggplot(data = subset(OTU_nsp_order_melt, elevation > 800)) +
  #geom_errorbar(aes(ymin = minS, ymax = maxS, x = elevation, colour = site), width = 0) +
  geom_point(aes(y = avgS, x = elevation, fill = site), colour = "black",
             pch = 21, size = 3) +
  geom_smooth(aes(y = avgS, x = elevation, group = site, colour = site),
              se = TRUE, method = "loess", alpha = 0.2, size = 1) +
  facet_wrap(~variable, nrow =2, ncol = 3, scales = "free_y") +
  guides(colour = FALSE, fill = guide_legend(title = NULL, override.aes = list(size = 5))) +
  scale_fill_manual(values = c("#003262", "#FDB515")) +
  scale_colour_manual(values = c("#003262", "#FDB515")) +
  scale_x_continuous(breaks = c(800,1000,1200,1400, 1600)) +
  labs(x = "Elevation (m)",  y= "No. of OTUs") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = "bottom")
ggsave(nsp_order_plot, filename = file.path(fig.dir, "nsp_order_plot.pdf"), height = 5, width = 8)
ggsave(nsp_order_plot, filename = file.path(fig.dir, "nsp_order_plot.png"), height = 12.5, width = 20, units = cm)


OTU_order<- ddply(.data = OTU_mat_nsp_site, .variables = .(site.id), summarize, 
                  avgS_araneae = mean(S_araneae / S),
                  avgS_coleoptera = mean(S_coleoptera / S),
                  avgS_hemiptera = mean(S_hemiptera / S),
                  avgS_lepidoptera = mean(S_lepidoptera/ S),
                  avgS_psocoptera = mean(S_psocoptera / S),
                  avgS_orthoptera = mean(S_orthoptera / S),
                  avgpi_araneae = mean(pi_araneae),
                  avgpi_coleoptera = mean(pi_coleoptera),
                  avgpi_hemiptera = mean(pi_hemiptera),
                  avgpi_lepidoptera = mean(pi_lepidoptera),
                  avgpi_psocoptera = mean(pi_psocoptera),
                  avgpi_orthoptera = mean(pi_orthoptera),
                  elevation = elevation[1],
                  site = site[1],
                  .progress = "text")

OTU_order <- subset(OTU_order, elevation >= 800)

# Plot proportion of species with elevation (Fig S4a) ) ========
OTU_order_S <- melt(OTU_order, id.vars = c("site.id", "site", "elevation"),
                    measure.vars = c("avgS_araneae", "avgS_coleoptera", "avgS_hemiptera", "avgS_lepidoptera", "avgS_psocoptera", "avgS_orthoptera"))
OTU_order_S$variable <- factor(OTU_order_S$variable, 
                               levels = c("avgS_araneae", "avgS_coleoptera", "avgS_hemiptera", "avgS_lepidoptera", "avgS_orthoptera", "avgS_psocoptera"),
                               labels = c("Araneae", "Coleoptera", "Hemiptera", "Lepidoptera", "Orthoptera", "Psocoptera"))

OTU_order_S$site.id <- factor(OTU_order_S$site.id, 
                              levels = unique(OTU_order_S[order(OTU_order_S$elevation),]$site.id)) # Order sites by elevation
OTU_order_S_plot <- ggplot(data = OTU_order_S) +
  geom_bar(aes(x = factor(site.id), fill = variable, weight = value)) + 
  facet_wrap(~site, scales = "free_x") +
  labs(y = "Proportion of species", x = "Site ID") +
  #guides(fill = guide_legend(nrow = 1)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        legend.position = "right") +
  scale_fill_brewer(palette = "Spectral")

# Plot proportion of reads with elevation (Fig S4b) ========
OTU_order_pi <- melt(OTU_order, id.vars = c("site.id", "site", "elevation"),
                    measure.vars = c("avgpi_araneae", "avgpi_coleoptera", "avgpi_hemiptera", "avgpi_lepidoptera", "avgpi_psocoptera", "avgpi_orthoptera"))
OTU_order_pi$variable <- factor(OTU_order_pi$variable, 
                               levels = c("avgpi_araneae", "avgpi_coleoptera", "avgpi_hemiptera", "avgpi_lepidoptera", "avgpi_orthoptera", "avgpi_psocoptera"),
                               labels = c("Araneae", "Coleoptera", "Hemiptera", "Lepidoptera", "Orthoptera", "Psocoptera"))
OTU_order_pi$site.id <- factor(OTU_order_pi$site.id, 
                              levels = unique(OTU_order_pi[order(OTU_order_pi$elevation),]$site.id)) # Order sites by elevation

OTU_order_pi_plot <- ggplot(data = OTU_order_pi) +
  geom_bar(aes(x = factor(site.id), fill = variable, weight = value)) + 
  facet_wrap(~site, scales = "free_x") +
  labs(y = "Proportion of reads", x = "Site ID") +
  #guides(fill = guide_legend(nrow = 1)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        legend.position = "right") +
  scale_fill_brewer(palette = "Spectral")

OTU_order_plot <- 
  plot_grid(
  plot_grid(OTU_order_S_plot + theme(legend.position = "none"),
                            OTU_order_pi_plot + theme(legend.position = "none"), nrow = 2, labels = "auto"),
  get_legend(OTU_order_S_plot), nrow = 1, rel_widths = c(1, 0.2))
ggsave(OTU_order_plot, filename = file.path(fig.dir, "OTU_order_plot.pdf"), width = 12, height = 6)

## ELEVATIONAL TURNOVER  ================
# Plot weighted NMDS (Fig 3a) =========
site_bray_nmds <- read.csv(file.path(res.dir, "site_bray_nmds.csv"))
nmds_bray_plot <- 
  ggplot() + 
  geom_point(aes(x = MDS1, y = MDS2, shape = site, colour = elevation), data = site_bray_nmds, size = 3) +
  #ggrepel::geom_text_repel(aes(x = MDS1, y = MDS2, label = elevation), data = site_bray_nmds) +
  guides(colour = guide_colorbar(title = "Elevation (m)"), shape = guide_legend(title = "Site")) +
  scale_colour_viridis() +
  theme(#legend.title = element_blank(),
        panel.background = element_blank(), 
        legend.position = c(0,1),
        legend.justification = c(0,1))
ggsave(nmds_bray_plot, filename = file.path(fig.dir, "nmds_bray_plot.pdf"), width = 6, height = 6)

# Plot unweighted NMDS (Fig 3a) =========
site_nmds_bray_pa <- read.csv(file.path(res.dir, "site_bray_pa_nmds.csv"))
nmds_bray_pa_plot  <-
  ggplot() +
  geom_point(aes(x = MDS2, y = MDS1, shape = site, colour = elevation), data = site_nmds_bray_pa, size = 3) +
  #ggrepel::geom_text_repel(aes(x = MDS1, y = MDS2, label = elevation), data = site_nmds_bray_pa) +
  guides(colour = guide_colorbar(title = "Elevation (m)"), shape = guide_legend(title = "Site")) +
  scale_colour_viridis() +
  theme(panel.background = element_blank(),
    legend.position = c(0,1),
    legend.justification = c(0,1))
ggsave(nmds_bray_pa_plot, filename = file.path(fig.dir, "nmds_bray_pa_plot.pdf"), width = 6, height = 6)

# Prepare bray-curtis data
elevation_dist <- read.csv(file.path(res.dir, "elevation_dist.csv"))
elevation_dist_subset <- subset(elevation_dist, abs(elevation_dist$elevation_X1 - elevation_dist$elevation_X2) <= 100)
elevation_dist_subset$avgElevation <- rowMeans(elevation_dist_subset[,c("elevation_X1", "elevation_X2")])
elevation_dist_compare <- subset(elevation_dist_subset, site_X1 != site_X2)

# Plot bray-curtis (all orders) with elevation (Fig. 3b)  =========
all_cor <- cor.test(elevation_dist_compare$bc_all,
                    elevation_dist_compare$avgElevation,
                    method = "spearman", exact = FALSE)
plot_label <- as.character(round(all_cor$estimate,2))

elevation_dist_plot <-
  ggplot(data = elevation_dist_compare, aes(y = bc_all, x = avgElevation)) + 
  geom_point(fill = "white", colour = "black", pch = 21, size = 3) + 
  annotate("text", x = -Inf, y = -Inf,
           label = as.expression(bquote(rho==.(plot_label)~"*")), hjust = -0.3, vjust= -0.5, parse = T, size = 4) +
  labs(y = "Bray-Curtis dissimilarity\n(All orders)", x = "Mean elevation (m)") +
  xlab("Mean elevation (m)") + 
  geom_smooth(method = "lm", colour = "grey20", size = 0.5) +
  theme(panel.background = element_blank())

# Plot bray-curtis (all orders + presence-absence) with elevation (Fig. S5b)
all_cor_pa <- cor.test(elevation_dist_compare$bcpa_all,
                       elevation_dist_compare$avgElevation,
                       method = "spearman", exact = FALSE)
plot_pa_label <- as.expression(bquote(rho==.(as.character(round(all_cor_pa$estimate,2)))~"*"))
elevation_distpa_plot <-
  ggplot(data = elevation_dist_compare, aes(y = bcpa_all, x = avgElevation)) + 
  geom_point(fill = "white", colour = "black", pch = 21, size = 3) + 
  annotate("text", x = -Inf, y = -Inf,
             label = plot_pa_label, hjust = -0.3, vjust= -0.5, parse = T, size = 4) +
  labs(y = "Bray-Curtis dissimilarity\n(All orders)", x = "Mean elevation (m)") +
  xlab("Mean elevation (m)") + 
  geom_smooth(method = "lm", colour = "grey20", size = 0.5) +
  theme(panel.background = element_blank())

# Define function to calculate correlations
calculateCorr <- function(x){
  temp <- cor.test(x$bc, x$avgElevation, exact = FALSE, method = "spearman")
  corr <- temp$estimate
  p <- temp$p.value
  data.frame(corr, p)
}


if(p < 0.05){
  label <- as.expression(bquote(rho==.(as.character(round(cor, 2)))~"*"))  
} else {
  label <- as.expression(bquote(rho==.(as.character(round(cor, 2)))))
}
# Plot bray-curtis (individual orders) with elevation (Fig 3c) ==========
elevation_dist_compare_melt <- 
  melt(elevation_dist_compare,
     id.var = c("X2", "X1", "avgElevation"),
     value.name = "bc",
     measure.vars = c("bc_araneae", "bc_hemiptera", "bc_lepidoptera", "bc_coleoptera", "bc_orthoptera", "bc_psocoptera"))
elevation_dist_compare_melt$variable <- factor(elevation_dist_compare_melt$variable,
  levels = c("bc_araneae", "bc_hemiptera", "bc_lepidoptera", "bc_coleoptera", "bc_orthoptera", "bc_psocoptera"),
  labels = c("Araneae", "Hemiptera", "Lepidoptera", "Coleoptera", "Orthoptera", "Psocoptera"))

elevation_dist_order_corr <- ddply(elevation_dist_compare_melt, .variables = .(variable), .fun = calculateCorr)
elevation_dist_taxon_plot <- 
  ggplot(data = elevation_dist_compare_melt,
         aes(y = bc, x = avgElevation, colour = variable)) + 
    geom_point(fill = "white", pch = 21, size = 3) +
    geom_smooth(method = "lm") +
    geom_text(aes(label = paste("rho==",
                                as.character(round(corr, 2)),
                                ifelse(p < 0.05, "~'*'", "")),
                  x = -Inf, y = -Inf),
              data = elevation_dist_order_corr,
              hjust = -0.1, vjust= -0.5, size = 4,
              colour = "black", parse = T) +
  labs(y="Bray-Curtis dissimilarity", x="Mean elevation (m)") + 
  facet_wrap(~variable, scales = "free") +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        strip.text = element_text(size = 8)) +
  scale_color_manual(values = berkcol3[c(1,5,7,9,11,14)])

# Plot bray-curtis (individual orders + presence-absence) with elevation (Fig S5c) ==========
elevation_distpa_compare_melt <- 
  melt(elevation_dist_compare,
       id.var = c("X2", "X1", "avgElevation"),
       value.name = "bcpa",
       measure.vars = c("bcpa_araneae", "bcpa_hemiptera", "bcpa_lepidoptera", "bcpa_coleoptera", "bcpa_orthoptera", "bcpa_psocoptera"))
elevation_distpa_compare_melt$variable <- factor(elevation_distpa_compare_melt$variable,
                                               levels = c("bcpa_araneae", "bcpa_hemiptera", "bcpa_lepidoptera", "bcpa_coleoptera", "bcpa_orthoptera", "bcpa_psocoptera"),
                                               labels = c("Araneae", "Hemiptera", "Lepidoptera", "Coleoptera", "Orthoptera", "Psocoptera"))

elevation_distpa_order_corr <- ddply(elevation_distpa_compare_melt, .variables = .(variable), .fun = calculateCorr)

elevation_distpa_taxon_plot <- 
  ggplot(data = elevation_distpa_compare_melt,
         aes(y = bcpa, x = avgElevation, colour = variable)) + 
  geom_point(fill = "white", pch = 21, size = 3) +
  geom_smooth(method = "lm") +
  geom_text(aes(label = paste("rho==",
                              as.character(round(corr, 2)),
                              ifelse(p < 0.05, "~'*'", "")),
                x = -Inf, y = -Inf),
            data = elevation_distpa_order_corr,
            hjust = -0.1, vjust= -0.5, size = 4,
            colour = "black", parse = T) +
  labs(y="Bray-Curtis dissimilarity", x="Mean elevation (m)") + 
  facet_wrap(~variable, scales = "free") +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        strip.text = element_text(size = 8)) +
  scale_color_manual(values = berkcol3[c(1,5,7,9,11,14)])

# Combine plots (Fig 3) =======
distance_plot <- 
  plot_grid(
    plot_grid(nmds_bray_plot + coord_fixed() + theme(legend.position = "right"),
              elevation_dist_plot,
              labels = c("a","b"), label_size = 18, rel_widths = c(1, 0.8)),
    elevation_dist_taxon_plot, nrow = 2, labels = c("","c"), label_size = 18, rel_heights = c(0.8, 1))
ggsave(distance_plot, filename = file.path(fig.dir, "distance_plot.pdf"), width = 10, height = 10)

# Bray-curtis for presence-absence (Fig S5) =======
distancepa_plot <- 
  plot_grid(
    plot_grid(nmds_bray_pa_plot + coord_fixed() + theme(legend.position = "right"),
              elevation_distpa_plot, labels = c("a","b"), label_size = 18, rel_widths = c(1, 0.8)),
    elevation_distpa_taxon_plot, nrow = 2, labels = c("","c"), label_size = 18, rel_heights = c(0.8, 1))
ggsave(distancepa_plot, filename = file.path(fig.dir, "distancepa_plot.pdf"), width = 10, height = 10)

## RESISTANCE ANALYSIS ====================
resist <- read_excel(file.path(res.dir, "jairo_resistance_09072020.xlsx"))
resist <- as.data.frame(resist)

resist$Model <- factor(resist$Model,
                       levels = c("Distance",
                                  "tair_ann",
                                  "rfgrid_mm_state_ann",
                                  "rfgrid_mm_state_ann.tair_ann",
                                  "Null"),
                       labels = c("Geographic distance",
                                  "Temp",
                                  "Precip",
                                  "Temp + Precip",
                                  "Null"))

groups <- c("Araneae",
            "Coleoptera",
            "Hemiptera",
            "Lepidoptera",
            "Orthoptera",
            "Psocoptera")
resist$Order <- factor(resist$Order, 
                       levels = groups)
resist$weight <- round(as.numeric(resist$weight), 3)
resist$AICc <- round(as.numeric(resist$AICc), 3) 
resist$k <- as.numeric(resist$k)
resist$delta.AICc <- round(as.numeric(resist$delta.AICc), 3) 
resist$R2m <- round(as.numeric(resist$R2m), 3)
resist$R2c <- round(as.numeric(resist$R2c), 3)
resist$LL <- round(as.numeric(resist$LL),3)
resist$weight2 <- ifelse(resist$weight == 0, "< 0.01", as.character(resist$weight))
resist$Transect <- factor(resist$Site, levels = c("Laupahoehoe", "Stainback"), labels= c("L", "S"))
write.csv(resist[c("Rarefied_ID", "Transect", "Order", "Model", "R2m", "R2c", "delta.AICc", "weight2")],
          file.path(res.dir, "jairo_model_results_clean.csv"),
          row.names = FALSE)

# Calculate average Akaike weight across rarefied datasets

## check that akaike weights sum to one
# resist_sum <- ddply(.data = resist, .variables = .(Site, Order, Rarefied_ID),
#                     .fun = summarise, avgAkaikeWeight = sum(weight))
resist_summary <- ddply(.data = resist,
                        .variables = .(Site, Order, Model),
                        .fun = summarise,
                        avgAkaikeWeight = mean(weight),
                        nDelta = sum(delta.AICc == 0))

resist_summary2 <- ddply(resist_summary,
                         .variables = .(Site, Order),
                         .fun = summarise,
                         Model = Model,
                         avgAkaikeWeight = avgAkaikeWeight,
                         nDelta = nDelta,
                         pos = cumsum(avgAkaikeWeight) - (0.5*avgAkaikeWeight) )

# Plot akaike weights of resistance models (Fig. 5) ==========
resist_bar <- 
  ggplot(data = resist_summary2) + 
  geom_bar(aes(fill = factor(Model), x= Order, weight = avgAkaikeWeight)) + 
  #geom_text(aes(label = nDelta, y = 1-pos, x = Order), size = 3) +
  facet_wrap(~Site) + 
  guides(fill = guide_legend(title = "Model")) +
  labs(x = NULL, y = "Akaike weight") +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_brewer(palette = "YlGnBu", direction = -1) + 
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12, colour = "grey20"))

ggsave(resist_bar, filename = file.path(fig.dir, "resistance_bar.pdf"), height = 4, width = 8)
ggsave(resist_bar, filename = file.path(fig.dir, "resistance_bar.png"), height = 4, width = 8)


## NICHE CONSERVATISM ANALYSIS  ================
# Plot histogram of coefficients
t_med_corr_t5 <- read.csv(file.path(res.dir, "t_med_corr_t5.csv"))
rf_med_corr_t5 <- read.csv(file.path(res.dir, "rf_med_corr_t5.csv"))

t_med_corr_t5$Variable <- "Niche position (temperature)"
rf_med_corr_t5$Variable <- "Niche position (precipitation)"
t_med_corr_t5$Variable2 <- "All OTUs"
rf_med_corr_t5$Variable2 <- "All OTUs"

t_med_corr_t5_subset <- read.csv(file.path(res.dir, "t_med_corr_t5_subset.csv"))
rf_med_corr_t5_subset <- read.csv(file.path(res.dir, "rf_med_corr_t5_subset.csv"))
t_med_corr_t5_subset$Variable <- "Niche position (temperature)"
rf_med_corr_t5_subset$Variable <- "Niche position (precipitation)"
t_med_corr_t5_subset$Variable2 <- "OTUs in all rarefied datasets"
rf_med_corr_t5_subset$Variable2 <- "OTUs in all rarefied datasets"

all_corr_t5 <- Reduce("rbind", list(t_med_corr_t5, rf_med_corr_t5, t_med_corr_t5_subset, rf_med_corr_t5_subset))
nicheposcorr_t5 <- ggplot(data = all_corr_t5) + geom_histogram(aes(x = cor, fill = ifelse(p<0.05, "1", "0")), bins = 20) +
  facet_grid(Variable2~Variable) + labs(y = "Frequency", x = "Correlation coefficient") +
  scale_fill_manual(values = c("grey70", "darkgreen")) +
  theme(panel.background = element_blank(),
        legend.position = "none")
ggsave(nicheposcorr_t5, filename = file.path(fig.dir, "nicheposcorr_t5.pdf"), height= 5, width = 9)

t_med_corr_t10 <- read.csv(file.path(res.dir, "t_med_corr_t10.csv"))
rf_med_corr_t10 <- read.csv(file.path(res.dir, "rf_med_corr_t10.csv"))
t_med_corr_t10$Variable <- "Niche position (temperature)"
rf_med_corr_t10$Variable <- "Niche position (precipitation)"
t_med_corr_t10$Variable2 <- "All OTUs"
rf_med_corr_t10$Variable2 <- "All OTUs"

t_med_corr_t10_subset <- read.csv(file.path(res.dir, "t_med_corr_t10_subset.csv"))
rf_med_corr_t10_subset <- read.csv(file.path(res.dir, "rf_med_corr_t10_subset.csv"))
t_med_corr_t10_subset$Variable <- "Niche position (temperature)"
rf_med_corr_t10_subset$Variable <- "Niche position (precipitation)"
t_med_corr_t10_subset$Variable2 <- "OTUs in all rarefied datasets"
rf_med_corr_t10_subset$Variable2 <- "OTUs in all rarefied datasets"

all_corr_t10 <- Reduce("rbind", list(t_med_corr_t10, rf_med_corr_t10, t_med_corr_t10_subset, rf_med_corr_t10_subset))
nicheposcorr_t10 <- ggplot(data = all_corr_t10) + geom_histogram(aes(x = cor, fill = ifelse(p<0.05, "1", "0")), bins = 20) +
  facet_grid(Variable2~Variable) + labs(y = "Frequency", x = "Correlation coefficient") +
  scale_fill_manual(values = c("grey70", "darkgreen")) +
  theme(panel.background = element_blank(),
        legend.position = "none")
ggsave(nicheposcorr_t10, filename = file.path(fig.dir, "nicheposcorr_t10.pdf"), height= 5, width = 9)




# 

t_med_t5 <- readRDS(file.path(res.dir, "t_med_t5.rds"))
rf_med_t5 <- readRDS(file.path(res.dir, "rf_med_t5.rds"))

test <- ldply(.data = t_med_t5, .fun = function(x)
  ddply(subset(x, V4 != "Coleoptera"), .variables = .(V4), nicheConsTest)
)
table(t_med_t5[[1]]$V4)
ggplot(data = test) +
  geom_histogram(aes(x= cor, fill = ifelse(p<0.05, "1", "0")), bins = 20) +
  facet_wrap(.~V4)

ggplot(data = t_med_t5[[1]]) +
  geom_point(aes(y = Laupahoehoe, x = Stainback)) +
  facet_wrap(.~V4)

