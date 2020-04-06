## Generate figures
# Authors: Jun Ying Lim

## Packages ====================
rm(list = ls())
library(readxl)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(viridis)
library(wesanderson)
library(ggnewscale)
library(plyr)

## Directories ====================
main.dir <- "~/Dropbox/Projects/2017/hawaiiCommunityAssembly/elevationAssemblyHawaii/"
res.dir <- file.path(main.dir, "results")
fig.dir <- file.path(main.dir, "figures")
data.dir <- file.path(main.dir, "data")
siteData <- read.csv(file.path(data.dir, "clim.final.csv"), stringsAsFactors = FALSE)
siteData <- subset(siteData, island == "BigIsland")

## Resistance analysis ====================
resist <- read_excel(file.path(res.dir, "jairo_resistance_10032020.xlsx"))
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

groups <- c(#"All",
            "Araneae",
            "Coleoptera",
            "Hemiptera",
            "Lepidoptera",
            "Orthoptera",
            "Psocoptera")
resist$Order <- factor(resist$Order, 
                       levels = rev(groups))



# Calculate average akaike weight across rarefied datasets
resist$weight <- as.numeric(resist$weight)

# Jairo's error?
resist_sum <- ddply(.data = subset(resist, !(Rarefied_ID == 92 & Site == "Laupahoehoe" & Order == "Coleoptera")), .variables = .(Site, Order, Rarefied_ID), .fun = summarize, avgAkaikeWeight = sum(weight) )

resist_summary <- ddply(.data = subset(resist, !(Rarefied_ID == 92 & Site == "Laupahoehoe" & Order == "Coleoptera")),
                        .variables = .(Site, Order, Model), .fun = summarize, avgAkaikeWeight = mean(weight) )

resist_sum2 <- ddply(resist_summary, .variables = .(Site, Order), .fun = summarize, sumAW = sum(avgAkaikeWeight))

# Plotting as a heatmap
resist_summary$weight2 <- cut(resist_summary$avgAkaikeWeight, breaks = seq(0, 1, 1/8))

resist_plot <- ggplot(data = resist) +
  geom_tile(aes(x = Model, y = Group, fill = weight2)) + 
  facet_wrap(.~Site) + 
  scale_fill_manual(values = brewer.pal(name = "YlGnBu", n = 9)[8:1]) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size = 12, colour = "grey20"),
        panel.background = element_blank(), 
        axis.title.x = element_blank())

ylcolors <- rev(brewer.pal(name = "YlGnBu", n = 8))

resist_legend <- ggplot(data = resist) +
  geom_tile(aes(x = Model, y = Group, fill = weight)) +
  facet_wrap(.~Site) + 
  guides(fill = 
           guide_colorbar(nbin = 8, title = "Akaike Weight", label = T, raster = F,
                          draw.ulim = T, draw.llim = T, ticks = F,
                          barwidth = unit(2,"in"),
                          barheight = unit(0.1, "in"))) +
  scale_fill_gradientn(colours = ylcolors,
                       breaks = c(0, 0.5, 1),
                       limits = c(0,1),
                       labels = c(0, 0.5, 1)) +
  theme(legend.position = "bottom")

resist_plot_comb <- plot_grid(resist_plot + theme(legend.position = "none"),
                              get_legend(resist_legend), 
                              nrow = 2,
                              rel_heights = c(0.9, 0.1))
ggsave(resist_plot_comb, filename = file.path(fig.dir, "resistance.pdf"), width = 8, height = 6)

# Plotting as stacked bars
resist_bar <- ggplot(data = resist_summary) + 
  geom_bar(aes(fill = factor(Model), x= Order, weight = avgAkaikeWeight)) + 
  facet_wrap(~Site) + 
  guides(fill = guide_legend(title = "Model")) +
  labs(x = NULL, y = "Akaike weight") +
  scale_y_continuous(expand = c(0,0)) +
  #scale_fill_viridis(discrete = T) +
  scale_fill_brewer(palette = "YlGnBu", direction = -1) + 
  #scale_fill_manual(values = wes_palette(name = "Cavalcanti1", n = 5)) +
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12, colour = "grey20"))

ggsave(resist_bar, filename = file.path(fig.dir, "resistance_bar.pdf"), height = 7, width = 10)

## Site patterns =========
rf_site_plot <- ggplot(data = subset(siteData, site %in%  c("Laupahoehoe", "Stainback"))) +
  geom_point(aes(y = rf_ann, x = elevation, fill = site ), pch = 21, size = 3, colour = "grey20" ) +
  scale_fill_manual(values = c("#003262", "#FDB515")) +
  labs(y = "Total annual rainfall (mm)", x = "Elevation (m)") +
  guides(fill = guide_legend(title = NULL, override.aes = list(size = 5))) +
  theme(panel.background = element_blank(),
        legend.position = c(1, 1),
        #legend.background = element_rect(colour = "grey20"),
        legend.justification = c(1,1),
        legend.text = element_text(size = 10))

temp_site_plot <- ggplot(data = subset(siteData, site %in%  c("Laupahoehoe", "Stainback"))) +
  geom_point(aes(y = t_ann, x = elevation, fill = site ), pch = 21, size = 3, colour = "grey20" ) +
  scale_fill_manual(values = c("#003262", "#FDB515")) +
  labs(y = expression("Mean annual temperature ("~degree*C~")" ), x = "Elevation (m)") +
  guides(fill = guide_legend(title = NULL, override.aes = list(size = 5))) +
  theme(panel.background = element_blank(),
        legend.position = c(1, 1),
        #legend.background = element_rect(colour = "grey20"),
        legend.justification = c(1,1),
        legend.text = element_text(size = 10))


## Nsp ================
library(plyr)
OTU_mat_nsp_site <- readRDS(file.path(res.dir, "OTU_nsp_site.rds"))
OTU_mat_nsp_site_avg <- ddply(OTU_mat_nsp_site, .variables = .(site_id), summarize, mean_nsp = mean(nsp))
OTU_mat_nsp_site_avg <- merge(x = OTU_mat_nsp_site_avg, y= siteData, by.x = "site_id", by.y = "site.id")

nsp_plot <- ggplot() +
  geom_point(aes(y = nsp, x = elevation, colour = site),
             alpha = 0.05, pch = 16, size = 3, data = OTU_mat_nsp_site) +
  geom_point(aes(y = mean_nsp, x = elevation, fill = site, colour = site), data = OTU_mat_nsp_site_avg,
             pch = 21, size = 3, colour = "grey20") +
  geom_smooth(aes(y = nsp, x = elevation, group = site, colour = site),
              se = TRUE, method = "loess", alpha = 0.1, size = 1,
              data = OTU_mat_nsp_site) +
  guides(colour = FALSE, fill = guide_legend(title = NULL, override.aes = list(size = 5))) +
    scale_fill_manual(values = c("#003262", "#FDB515")) +
  scale_colour_manual(values = c("#003262", "#FDB515")) +
  labs(x = "Elevation (m)",  y= "No. of unique\nOTUs") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.background = element_blank(),
        legend.position = c(0,1),
        legend.justification = c(0,1))

# GEOGRAPHIC MAP =====================
library(rgdal)
library(raster)
library(RStoolbox)
hillshade <- raster(file.path(data.dir, "MHI_digital_elevation_model_hillshade_GIS_data/Hawaii_Hillshade/Hawaii_Hillshade_45deg.img"))
siteSP <- SpatialPointsDataFrame(coords = siteData[c("longitude","latitude")], proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"), data = siteData[c("longitude","latitude", "site")])
siteSP_reproj <- spTransform(siteSP, CRSobj = CRS(proj4string(hillshade)))
siteSP_df <- cbind(siteSP_reproj@data["site"],siteSP_reproj@coords)

hillshade_fort <- fortify(hillshade, maxpixels = 100000) #197,920,112

map_plot <- ggplot() +
  geom_tile(aes(y = y, x = x, fill = Hawaii_Hillshade_45deg.1_BinValues), data = hillshade_fort) +
  coord_fixed() +
  scale_fill_gradient(low = grey(0.01), high = grey(1), na.value = "white") +
  ggnewscale::new_scale_fill() +
  geom_point(aes(y= latitude, x= longitude, fill = site, colour = site), data = siteSP_df,
             pch = 21, size = 3.5) +
  scale_colour_manual(values = c("white", "grey20")) +
  geom_label(aes(x = 285000, y= 2210000, label = "Laupahoehoe"), size = 3) +
  geom_label(aes(x = 270000, y= 2155000, label = "Stainback"), size = 3) +
  scale_fill_manual(values = c("#003262", "#FDB515")) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")
ggsave(file.path(fig.dir, "map.pdf"), map_plot, width = 7, height = 7)

## FIGURE 1  ================

fig1 <- plot_grid(map_plot, nsp_plot, rf_site_plot, temp_site_plot, labels = "auto")
ggsave(fig1, filename = file.path(fig.dir, "fig1.pdf"), height = 8, width = 8)

## NMDS  ================
library(viridis); library(directlabels)
site_nmds <- read.csv(file.path(res.dir, "site_nmds.csv"))
site_nmds_contour <- read.csv(file.path(res.dir, "site_nmds_contour.csv"))

nmds_plot <- ggplot() + 
  geom_point(aes(x = MDS1, y = MDS2, fill = site, shape = site), data = site_nmds, size = 3) +
  scale_colour_manual(values = c("#003262", "#FDB515")) +
  new_scale_color() +
  stat_contour(data = site_nmds_contour, aes(x, y, z = z, group = ..level.., colour = ..level..)) +
  #scale_colour_brewer(palette = "BrBG")
  scale_colour_viridis() +
  theme(legend.title = element_blank(),
        panel.background = element_blank(), 
        legend.position = c(1,1),
        legend.justification = c(1,1))
nmds_contour_plot <- direct.label(nmds_plot, "top.pieces")
ggsave(nmds_contour_plot, filename = file.path(fig.dir, "nmds_plot.pdf"), width = 6, height = 6)

elevation_dist <- read.csv(file.path(res.dir, "elevation_dist.csv"))

elevation_dist_plot <- ggplot(data = subset(elevation_dist, elevation_X1_class == elevation_X2_class & site_X1 != site_X2)) +
  geom_violin(aes(y = dist, x = elevation_X1_class, fill = elevation_X1_class)) +
  geom_point(aes(y = dist, x = elevation_X1_class)) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_blank()) +
  labs(y = "Pairwise Bray-Curtis distance", x = "Elevation bin (m)") +
  scale_fill_manual(values = wes_palette("IsleofDogs2", n = 4))

distance_plot <- plot_grid(nmds_contour_plot, elevation_dist_plot, labels = "auto", label_size = 18)
ggsave(distance_plot, filename = file.path(fig.dir, "distance_plot.pdf"), width = 10, height = 5)

## Niche position ================
t_med_site_t10 <- read.csv(file.path(res.dir, "t_med_res_prep_avg_t10.csv"))
t_med_site_t5 <- read.csv(file.path(res.dir, "t_med_res_prep_avg_t5.csv"))

rf_med_site_t10 <- read.csv(file.path(res.dir, "rf_med_res_prep_avg_t10.csv"))
rf_med_site_t5 <- read.csv(file.path(res.dir, "rf_med_res_prep_avg_t5.csv"))


t_med_plot_t10 <- ggplot(aes(y = Laupahoehoe, x = Stainback),data = t_med_site_t10) + 
  geom_point(aes(fill = V4),
             colour = "grey20", size = 3, pch = 21) +
  #geom_smooth(method = "lm") +
  #geom_abline(aes(slope = 1, intercept = 0)) +
  labs(y = "Niche position\n(Laupahoehoe)",
       x = "Niche position\n(Stainback)") +
  scale_fill_brewer(palette = "Spectral") +
  scale_colour_brewer(palette = "Spectral") +
  guides(fill = guide_legend(title = "Taxonomic order")) + 
  theme(plot.background = element_blank(),
        panel.background = element_blank())

rf_med_plot_t10 <- ggplot(data = rf_med_site_t10) + 
  geom_point(aes(y = Laupahoehoe, x = Stainback, fill = V4),
             colour = "grey20", size = 3, pch = 21) +
  #geom_abline(aes(slope = 1, intercept = 0)) +
  labs(y = "Niche position\n(Laupahoehoe)",
       x = "Niche position\n(Stainback)") +
  scale_fill_brewer(palette = "Spectral") +
  scale_colour_brewer(palette = "Spectral") +
  theme(plot.background = element_blank(),
        panel.background = element_blank())

t_med_plot_t5 <- ggplot(aes(y = Laupahoehoe, x = Stainback),data = t_med_site_t5) + 
  geom_point(aes(fill = V4),
             colour = "grey20", size = 3, pch = 21) +
  #geom_smooth(method = "lm") +
  #geom_abline(aes(slope = 1, intercept = 0)) +
  labs(y = "Niche position\n(Laupahoehoe)",
       x = "Niche position\n(Stainback)") +
  scale_fill_brewer(palette = "Spectral") +
  scale_colour_brewer(palette = "Spectral") +
  guides(fill = guide_legend(title = "Taxonomic order")) + 
  theme(plot.background = element_blank(),
        panel.background = element_blank())

rf_med_plot_t5 <- ggplot(data = rf_med_site_t5) + 
  geom_point(aes(y = Laupahoehoe, x = Stainback, fill = V4),
             colour = "grey20", size = 3, pch = 21) +
  #geom_abline(aes(slope = 1, intercept = 0)) +
  labs(y = "Niche position\n(Laupahoehoe)",
       x = "Niche position\n(Stainback)") +
  scale_fill_brewer(palette = "Spectral") +
  scale_colour_brewer(palette = "Spectral") +
  theme(plot.background = element_blank(),
        panel.background = element_blank())


## Climate breadth ================
t_breadth_site_t10 <- read.csv(file.path(res.dir, "t_bdt_res_prep_avg_t10.csv"))
rf_breadth_site_t10 <- read.csv(file.path(res.dir, "rf_bdt_res_prep_avg_t10.csv"))

t_breadth_site_t5 <- read.csv(file.path(res.dir, "t_bdt_res_prep_avg_t5.csv"))
rf_breadth_site_t5 <- read.csv(file.path(res.dir, "rf_bdt_res_prep_avg_t5.csv"))

t_breadth_plot_t10 <- ggplot(data = t_breadth_site_t10) + 
  geom_point(aes(y = Laupahoehoe, x = Stainback, fill = V4),
             colour = "grey20", size = 3, pch = 21) +
  labs(y = "Niche breadth\n(Laupahoehoe)",
       x = "Niche breadth\n(Stainback)") +
  scale_fill_brewer(palette = "Spectral") +
  scale_colour_brewer(palette = "Spectral") +
  theme(plot.background = element_blank(),
        panel.background = element_blank())

rf_breadth_plot_t10 <- ggplot(data = rf_breadth_site_t10) +
  geom_point(aes(y = Laupahoehoe, x = Stainback, fill = V4),
             colour = "grey20", size = 3, pch = 21) +
  labs(y = "Niche breadth\n(Laupahoehoe)",
       x = "Niche breadth\n(Stainback)") +
  scale_fill_brewer(palette = "Spectral") +
  scale_colour_brewer(palette = "Spectral") +
  theme(plot.background = element_blank(),
        panel.background = element_blank())

t_breadth_plot_t5 <- ggplot(data = t_breadth_site_t5) + 
  geom_point(aes(y = Laupahoehoe, x = Stainback, fill = V4),
             colour = "grey20", size = 3, pch = 21) +
  labs(y = "Niche breadth\n(Laupahoehoe)",
       x = "Niche breadth\n(Stainback)") +
  scale_fill_brewer(palette = "Spectral") +
  scale_colour_brewer(palette = "Spectral") +
  theme(plot.background = element_blank(),
        panel.background = element_blank())

rf_breadth_plot_t5 <- ggplot(data = rf_breadth_site_t5) +
  geom_point(aes(y = Laupahoehoe, x = Stainback, fill = V4),
             colour = "grey20", size = 3, pch = 21) +
  labs(y = "Niche breadth\n(Laupahoehoe)",
       x = "Niche breadth\n(Stainback)") +
  scale_fill_brewer(palette = "Spectral") +
  scale_colour_brewer(palette = "Spectral") +
  theme(plot.background = element_blank(),
        panel.background = element_blank())

## Combined
library(cowplot)

niche_cons_plot_t10 <- plot_grid(plot_grid(t_med_plot_t10 + geom_abline(aes(intercept = 0, slope = 1)) + theme(legend.position = "none"),
                                       rf_med_plot_t10 + geom_abline(aes(intercept = 0, slope = 1)) + theme(legend.position = "none"),
                                       t_breadth_plot_t10+ geom_abline(aes(intercept = 0, slope = 1)) + theme(legend.position = "none"),
                                       rf_breadth_plot_t10+ geom_abline(aes(intercept = 0, slope = 1)) + theme(legend.position = "none"), labels = "auto", scale = 0.9),
                             get_legend(t_med_plot_t10 + theme(legend.position = "bottom")), nrow = 2, rel_heights = c(0.9,0.1))

ggsave(niche_cons_plot_t10, filename=file.path(fig.dir, "niche_cons_t10.pdf"), height = 7, width = 7)

niche_cons_plot_t5 <- plot_grid(plot_grid(t_med_plot_t5 + geom_abline(aes(intercept = 0, slope = 1))+ theme(legend.position = "none"),
                                       rf_med_plot_t5 + geom_abline(aes(intercept = 0, slope = 1))+ theme(legend.position = "none"),
                                       t_breadth_plot_t5 + geom_abline(aes(intercept = 0, slope = 1))+ theme(legend.position = "none"),
                                       rf_breadth_plot_t5+ geom_abline(aes(intercept = 0, slope = 1)) + theme(legend.position = "none"), labels = "auto", scale = 0.9),
                             get_legend(t_med_plot_t5 + theme(legend.position = "bottom")), nrow = 2, rel_heights = c(0.9,0.1))
ggsave(niche_cons_plot_t5, filename=file.path(fig.dir, "niche_cons_t5.pdf"), height = 7, width = 7)


## Niche breadth taxonomic orders
rf_bdt_melt_t10 <- read.csv(file.path(res.dir, "rf_bdt_melt_t10.csv"))
rf_bdt_melt_t5 <- read.csv(file.path(res.dir, "rf_bdt_melt_t5.csv"))
t_bdt_melt_t10 <- read.csv(file.path(res.dir, "t_bdt_melt_t10.csv"))
t_bdt_melt_t5 <- read.csv(file.path(res.dir, "t_bdt_melt_t5.csv"))


t_breadth_plot_t5 <- ggplot(aes(y = t_breadth, x = V4), data = t_bdt_melt_t5) +
  geom_boxplot(aes(fill = V4)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~Site) +
  guides(fill = guide_legend(title = "Taxonomic order")) +
  labs(y = "Niche breadth\n(Mean annual temperature)") +
  scale_fill_brewer(palette = "Spectral") +
  theme(legend.position = "bottom",
        axis.text = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())

rf_breadth_plot_t5 <- ggplot(aes(y = rf_breadth, x = V4), data = rf_bdt_melt_t5) +
  geom_boxplot(aes(fill = V4)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~Site) +
  guides(fill = guide_legend(title = "Taxonomic order")) +
  labs(y = "Niche breadth\n(Total annual precipitation)") +
  scale_fill_brewer(palette = "Spectral") +
  theme(legend.position = "bottom",
        axis.text = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())

t_breadth_plot_t10 <- ggplot(aes(y = t_breadth, x = V4), data = t_bdt_melt_t10) +
  geom_boxplot(aes(fill = V4)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~Site) +
  guides(fill = guide_legend(title = "Taxonomic order")) +
  labs(y = "Niche breadth\n(Mean annual temperature)") +
  scale_fill_brewer(palette = "Spectral") +
  theme(legend.position = "bottom",
        axis.text = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())

rf_breadth_plot_t10 <- ggplot(aes(y = rf_breadth, x = V4), data = rf_bdt_melt_t10) +
  geom_boxplot(aes(fill = V4)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~Site) +
  guides(fill = guide_legend(title = "Taxonomic order")) +
  labs(y = "Niche breadth\n(Total annual precipitation)") +
  scale_fill_brewer(palette = "Spectral") +
  theme(legend.position = "bottom",
        axis.text = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())

niche_bdt_plot_t10 <- plot_grid(plot_grid(t_breadth_plot_t10 + theme(legend.position = "none"),
                                          rf_breadth_plot_t10+ theme(legend.position = "none"), labels = "auto"),
                                get_legend(t_breadth_plot_t10), rel_heights = c(0.9, 0.1), nrow = 2)

niche_bdt_plot_t5 <- plot_grid(plot_grid(t_breadth_plot_t5 + theme(legend.position = "none"),
                                         rf_breadth_plot_t5+ theme(legend.position = "none"), labels = "auto"),
                               get_legend(t_breadth_plot_t5), rel_heights = c(0.9, 0.1), nrow = 2)
ggsave(niche_bdt_plot_t10, filename = file.path(fig.dir, "niche_bdt_plot_t10.pdf"), height = 4, width = 6)
ggsave(niche_bdt_plot_t5, filename = file.path(fig.dir, "niche_bdt_plot_t5.pdf"), height = 4, width = 6)
