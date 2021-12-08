
library(ggplot2)
library(RColorBrewer)

set.seed(12345)

data_list <- list()
for(i in 1:4){
  y <- paste0("OTU", 1:5)
  x <- factor(c("Site A", "Site B", "Site C"), levels = paste("Site", LETTERS[3:1]))
  test <- expand.grid(y,x)
  z <- sample(x = c(0,1), prob = c(0.2, 0.8), replace = T, size = 15)
  #z[z > 0] <- rpois(n = sum(z>0), lambda = 10)
  z[z > 0] <- round(rlnorm(n = sum(z>0), meanlog = 4), 0)
  test$z <- z
  data_list[[i]] <- test
}

plot_list <- list()
width = 3
height = 1.5
for(i in 1:4){
  plot_list[[i]] <- ggplot(data = subset(data_list[[i]], Var2 == "Site A")) +
    geom_tile(aes(y = Var2, x = Var1, alpha = z)) + 
    geom_text(aes(y = Var2, x = Var1, label = z), col = "white") + 
    scale_x_discrete(expand = c(0,0), position = "top") +
    scale_y_discrete(expand = c(0,0)) +
    scale_alpha(range(0.1, 1), limits = c(0, 287), trans = "sqrt") +
    coord_fixed() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank())
  ggsave(plot_list[[i]],
         filename = file.path("~/Dropbox/projects/2017/hawaiiCommunityAssembly/elevationAssemblyHawaii/figures", paste0("fig1_otu", i, ".pdf")),
         width = width, height = height)
}



sum(subset(data_list[[1]], Var2 == "Site A")$z)
sum(subset(data_list[[2]], Var2 == "Site A")$z)
sum(subset(data_list[[3]], Var2 == "Site A")$z)
sum(subset(data_list[[4]], Var2 == "Site A")$z)

## MOCK RAREFACTION
# nSpecimen_list <- c(15, 8, 12, 5)
# rarefied_list_list <- list()
# 
# for(j in 1:3){
#   rarefied_list <- list()
#   for(i in 1:4){
#     print(i)
#     df = subset(data_list[[i]], Var2 == "Site A")
#     rawReadAbund <- reshape2::acast(df, Var1~Var2, value.var = "z")
#     nSpecimen <- nSpecimen_list[i] 
#     rarefiedReadAbund <- vegan::rrarefy(rawReadAbund, sample = nSpecimen*7) 
#     rarefied_list[[i]] <- reshape2::melt(rarefiedReadAbund)
#   }
#   rarefied_list_list[[j]] <- rarefied_list
# }

nSpecimen_list <- c(3, 8, 12, 5)
rarefied_list <- list()
for(i in 1:4){
  print(i)
  #i = 1
  df = subset(data_list[[i]], Var2 == "Site A")
  rawReadAbund <- reshape2::acast(df, Var1~Var2, value.var = "z")
  nSpecimen <- nSpecimen_list[i]
  rarefiedReadAbund <- vegan::rrarefy(rawReadAbund, sample = nSpecimen*20)
  rarefied_list[[i]] <- reshape2::melt(rarefiedReadAbund)
}


# First rarefaction
for(i in 1:4){
  x <- ggplot(data = rarefied_list[[i]]  ) +
    geom_tile(aes(y = Var1, x = Var2, alpha = value)) + 
    geom_text(aes(y = Var1, x = Var2, label = value), col = "white") + 
    scale_x_discrete(expand = c(0,0), position = "top") +
    scale_y_discrete(expand = c(0,0)) +
    scale_alpha(range(0.1, 1), limits = c(0, 287), trans = "sqrt") +
    coord_fixed() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank())
  ggsave(x, filename = file.path("~/Dropbox/projects/2017/hawaiiCommunityAssembly/elevationAssemblyHawaii/figures", paste0("fig1_otu_rarefy",i, ".pdf")),
         width = width, height = height, bg = "transparent")
  
}


# combined <- do.call("rbind", rarefied_list_list[[1]])
# library(plyr)
# combined2 <- ddply(.data = combined, .variables = .(Var1,Var2), summarize, total = sum(value))
# x <- ggplot(data = combined2 ) +
#   geom_tile(aes(y = Var1, x = Var2, alpha = total)) + 
#   geom_text(aes(y = Var1, x = Var2, label = total), col = "white") + 
#   scale_x_discrete(expand = c(0,0), position = "top") +
#   scale_y_discrete(expand = c(0,0)) +
#   scale_alpha(range(0.1, 1), limits = c(0, 23)) +
#   coord_fixed() +
#   theme(legend.position = "none",
#         panel.grid = element_blank(),
#         panel.background = element_blank(),
#         axis.title = element_blank())
# ggsave(x, filename = "~/Desktop/x.pdf", width = 3, height = 1.5, bg = "transparent")

combined <- do.call("rbind", rarefied_list)
library(plyr)
combined2 <- ddply(.data = combined, .variables = .(Var1,Var2), summarize, total = sum(value))
x <- ggplot(data = combined2 ) +
  geom_tile(aes(y = Var1, x = Var2, alpha = total)) +
  geom_text(aes(y = Var1, x = Var2, label = total), col = "white") +
  scale_x_discrete(expand = c(0,0), position = "top") +
  scale_y_discrete(expand = c(0,0)) +
  scale_alpha(range(0.1, 1), limits = c(0, 287), trans = "sqrt") +
  coord_fixed() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank())
ggsave(x, filename = "~/Dropbox/projects/2017/hawaiiCommunityAssembly/elevationAssemblyHawaii/figures/fig1_otucomb.pdf", width = 3, height = 1.5, bg = "transparent")

