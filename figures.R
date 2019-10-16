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

## Directories ====================
main.dir <- "~/Dropbox/Projects/2017/hawaiiCommunityAssembly/elevationAssemblyHawaii/"
res.dir <- file.path(main.dir, "results")
fig.dir <- file.path(main.dir, "figures")
data.dir <- file.path(main.dir, "data")
siteData <- read.csv(file.path(data.dir, "clim.final.csv"), stringsAsFactors = FALSE)

## Resistance analysis ====================
resist <- read_excel(file.path(res.dir, "jairo_resistance.xlsx"))
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
groups <- c("All",
            "Araneae",
            "Coleoptera",
            "Hemiptera",
            "Lepidoptera",
            "Orthoptera",
            "Psocoptera")
resist$Group <- factor(resist$Group, 
                       levels = rev(groups))

resist$weight2 <- cut(resist$weight, breaks = seq(0, 1, 1/8))

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

# Plot as stacked bars
resist_bar <- ggplot(data = resist) + 
  geom_bar(aes(fill = factor(Model), x= Group, weight = weight)) + 
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

ggsave(resist_bar, filename = file.path(fig.dir, "resistance_bar.pdf"), height = )

## Site patterns =========
ggplot(data = subset(siteData, island == "BigIsland")) +
   geom_point(aes(y = rf_ann, x = elevation)) + facet_wrap(~ site)
 
ggplot(data = subset(siteData, island == "BigIsland")) +
   geom_point(aes(y = t_ann, x = elevation)) + facet_wrap(~ site)
 
ggplot(data = subset(siteData, island == "BigIsland")) +
   geom_point(aes(y = t_ann, x = rf_ann, color = site))

## Nsp ================
site_nsp <- read.csv(file.path(res.dir, "site_nsp.csv"))
site_nsp <- merge(site_nsp, siteData)
head(site_nsp)
ggplot(aes(y = nsp, x = elevation, colour = site), data = site_nsp) + geom_point() + geom_smooth(formula = y ~ poly(x,2), se = FALSE, method = "lm")
?geom_smooth

## Median climate ================
t_med_site <- read.csv(file.path(res.dir, "t_median_site.csv"))
rf_med_site <- read.csv(file.path(res.dir, "rf_median_site.csv"))
ggplot(data = t_med_site) + geom_point(aes(y = Laupahoehoe, x = Stainback, fill = V4),
                      colour = "grey20", size = 3, pch = 21) +
  geom_abline(aes(slope = 1, intercept = 0)) +
  #geom_smooth(aes(y = Laupahoehoe, x = Stainback, group = V4, colour = V4),
  #            method = "lm", se = FALSE) +
  labs(y = "Median annual temperature\n(Laupahoehoe)",
       x = "Median annual temperature\n(Stainback)") +
  scale_fill_brewer(palette = "Spectral") +
  scale_colour_brewer(palette = "Spectral") +
  theme(plot.background = element_blank(),
        panel.background = element_blank())

ggplot(data = rf_med_site) + geom_point(aes(y = Laupahoehoe, x = Stainback, fill = V4),
                                       colour = "grey20", size = 3, pch = 21) +
  geom_abline(aes(slope = 1, intercept = 0)) +
  #geom_smooth(aes(y = Laupahoehoe, x = Stainback, group = V4, colour = V4),
  #            method = "lm", se = FALSE) +
  labs(y = "Median annual precipitation\n(Laupahoehoe)",
       x = "Median annual precipitation\n(Stainback)") +
  scale_fill_brewer(palette = "Spectral") +
  scale_colour_brewer(palette = "Spectral") +
  theme(plot.background = element_blank(),
        panel.background = element_blank())

## Climate breadth ================
t_breadth_site <- read.csv(file.path(res.dir, "t_breadth_site.csv"))
rf_breadth_site <- read.csv(file.path(res.dir, "rf_breadth_site.csv"))

ggplot(data = t_breadth_site) + geom_point(aes(y = Laupahoehoe, x = Stainback, fill = V4),
                                       colour = "grey20", size = 3, pch = 21) +
  geom_abline(aes(slope = 1, intercept = 0)) +
  #geom_smooth(aes(y = Laupahoehoe, x = Stainback, group = V4, colour = V4),
  #            method = "lm", se = FALSE) +
  labs(y = "Range in annual temperature\n(Laupahoehoe)",
       x = "Range in annual temperature\n(Stainback)") +
  scale_fill_brewer(palette = "Spectral") +
  scale_colour_brewer(palette = "Spectral") +
  theme(plot.background = element_blank(),
        panel.background = element_blank())

ggplot(data = rf_breadth_site) + geom_point(aes(y = Laupahoehoe, x = Stainback, fill = V4),
                                           colour = "grey20", size = 3, pch = 21) +
  geom_abline(aes(slope = 1, intercept = 0)) +
  #geom_smooth(aes(y = Laupahoehoe, x = Stainback, group = V4, colour = V4),
  #            method = "lm", se = FALSE) +
  labs(y = "Range in annual precipitation\n(Laupahoehoe)",
       x = "Range in annual precipitation\n(Stainback)") +
  scale_fill_brewer(palette = "Spectral") +
  scale_colour_brewer(palette = "Spectral") +
  theme(plot.background = element_blank(),
        panel.background = element_blank())

## 
