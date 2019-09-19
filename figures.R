library(readxl)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

main.dir <- "~/Dropbox/Projects/2017/hawaiiCommunityAssembly/elevationAssemblyHawaii/"
res.dir <- file.path(main.dir, "results")
fig.dir <- file.path(main.dir, "figures")

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
