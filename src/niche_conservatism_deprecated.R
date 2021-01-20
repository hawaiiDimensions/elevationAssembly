## QUESTION 2: IS NICHE BREADTH CONSERVED ACROSS SITES? ========

t_bdt_t10 <- llply(.data = t_res_list_t10_subset,
                   .fun = nicheConsPrep, value = "temp_breadth", otu_tab = OTUtaxonData)
rf_bdt_t10 <- llply(.data = rf_res_list_t10_subset,
                    .fun = nicheConsPrep, value = "rf_breadth", otu_tab = OTUtaxonData)
saveRDS(t_bdt_t10, file.path(res.dir, "t_bdt_t10.rds"))
saveRDS(rf_bdt_t10, file.path(res.dir, "rf_bdt_t10.rds"))

t_bdt_corr_t10 <- ldply(.data = t_bdt_t10, .fun = nicheConsTest)
rf_bdt_corr_t10 <- ldply(.data = rf_bdt_t10, .fun = nicheConsTest)

write.csv(t_bdt_corr_t10, file.path(res.dir, "t_bdt_corr_t10.csv"), row.names = F)
write.csv(rf_bdt_corr_t10, file.path(res.dir, "rf_bdt_corr_t10.csv"), row.names = F)

round(mean(t_bdt_res_final_t10$cor), 2)
round(range(t_bdt_res_final_t10$cor), 2)
sum(t_bdt_res_final_t10$p < 0.05) / 100
round(range(t_bdt_res_final_t10$n), 2)

round(mean(rf_bdt_res_final_t10$cor), 2)
round(range(rf_bdt_res_final_t10$cor), 2)
sum(rf_bdt_res_final_t10$p < 0.05) / 100
round(range(rf_bdt_res_final_t10$n), 2)

# Re-run analysis with a lower site threshold
t_bdt_t5 <- llply(.data = t_res_list_t5_subset,
                  .fun = nicheConsPrep, value = "temp_breadth", otu_tab = OTUtaxonData)
rf_bdt_t5 <- llply(.data = rf_res_list_t5_subset,
                   .fun = nicheConsPrep, value = "rf_breadth", otu_tab = OTUtaxonData)
saveRDS(t_bdt_t5, file.path(res.dir, "t_bdt_t5.rds"))
saveRDS(rf_bdt_t5, file.path(res.dir, "rf_bdt_t5.rds"))

t_bdt_corr_t5 <- ldply(.data = t_bdt_t5, .fun = nicheConsTest)
rf_bdt_corr_t5 <- ldply(.data = rf_bdt_t5, .fun = nicheConsTest)

write.csv(t_bdt_corr_t5, file.path(res.dir, "t_bdt_corr_t5.csv"), row.names = F)
write.csv(rf_bdt_corr_t5, file.path(res.dir, "rf_bdt_corr_t5.csv"), row.names = F)

round(mean(t_bdt_res_final_t5$cor), 2)
round(range(t_bdt_res_final_t5$cor), 2)
sum(t_bdt_res_final_t5$p < 0.05) / 100
round(range(t_bdt_res_final_t5$n), 2)

round(mean(rf_bdt_res_final_t5$cor), 2)
round(range(rf_bdt_res_final_t5$cor), 2)
sum(rf_bdt_res_final_t5$p < 0.05) / 100
round(range(rf_bdt_res_final_t5$n), 2)



## MEHHHHH
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
