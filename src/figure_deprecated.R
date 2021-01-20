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
