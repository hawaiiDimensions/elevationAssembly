require(pracma)
require(picante)

RA <- function(x, ra_max, ra_m, ra_v, x_min, x_max){
  # Gaussian function describing relationship between relative abundance and an environmental variable
  # 
  # Args:
  #     x = environmental variable
  #     ra_max = maximum relative abundance
  #     ra_v = standard deviation of gaussian function
  #     ra_m = mode (environmental value where maximum relative abundance is achieved)
  #     x_min = minimum (lower bound, relative abundance is 0 below or equal to this value)
  #     x_max = maximum (upper bound, relative abundance is 0 above or equal to this value)
  res <- vector()
  for(i in 1:length(x)){
    if(x[i] >= x_max | x[i] <= x_min){
        res[i] <- 0
    } else {
        res[i] <- ra_max * exp( -0.5 * ((x[i] - ra_m) / ra_v)^2 )  
    } 
  }
  return(res)
}

RA_simple <- function(x, ra_max, ra_m, ra_v){
  res <- vector()
  for(i in 1:length(x)){
      res[i] <- ra_max * exp( -0.5 * ((x[i] - ra_m) / ra_v)^2 )  
  }
  return(res)
}


# x <- seq(0, 4, by = 0.01)
# ra_max = 0.5; ra_m = 0.3; ra_v = 1; x_min = -2; x_max = 2
# plot(RA(x = x, ra_max = ra_max, ra_m = ra_m, ra_v = ra_v, x_min = x_min, x_max = x_max) ~ x)
# abline(v = 0.3)

# Define least squares
RA_RSS <- function(data, par){
  with(data, sum((y - RA( x, ra_max = par[1],  ra_m = par[2], ra_v = par[3], x_min = par[4], x_max = par[5]))^2))
}

RA_simple_RSS <- function(data, par){
  with(data, sum((y - RA_simple( x, ra_max = par[1],  ra_m = par[2], ra_v = par[3]))^2))
}

fitGaussianModel <- function(data){
  # For each OTU, fits relative abundance as a gaussian function of temperature and rainfall
  # Designed to be iteratively run on subsets of OTUs
  # Args:
  #   data = data.frame containing OTU data
  
  #data <- subset(data, nReads > 0) 
  data$relabund <- data$nReads / sum(data$nReads)
  rf_abund <- data.frame(y = data$relabund, x = data$rf_ann)
  t_abund <- data.frame(y = data$relabund,x = data$t_ann)
  rf_init_par <- c( ra_max = max(data$relabund),
                    ra_m = data$rf_ann[which.max(data$relabund)],
                    ra_v = sd(data$rf_ann[data$relabund>0]),
                    xmin = min(data$rf_ann[data$relabund>0]),
                    xmax = max(data$rf_ann[data$relabund>0]))
  t_init_par <- c( ra_max = max(data$relabund),
                    ra_m = data$t_ann[which.max(data$relabund)],
                    ra_v = sd(data$t_ann[data$relabund>0]),
                    xmin = min(data$t_ann[data$relabund>0]),
                    xmax = max(data$t_ann[data$relabund>0]))
  rf_gauss_RSS <- optim(par = rf_init_par, fn = RA_RSS, method = "BFGS",
                        data = rf_abund)
  t_gauss_RSS <- optim(par = t_init_par, fn = RA_RSS, method = "BFGS",
                       data = t_abund)
  names(rf_gauss_RSS$par) <- paste0("rf_", names(rf_gauss_RSS$par))
  names(t_gauss_RSS$par) <- paste0("t_", names(t_gauss_RSS$par))
  gauss_para <- data.frame(t(c(rf_gauss_RSS$par, t_gauss_RSS$par)))
  # Generate predicted relative abundances
  data$rf_predict <- RA(x = rf_abund$x,
                        rf_gauss_RSS$par[1],
                        rf_gauss_RSS$par[2],
                        rf_gauss_RSS$par[3],
                        rf_gauss_RSS$par[4],
                        rf_gauss_RSS$par[5])
  data$t_predict <- RA(x = t_abund$x,
                       t_gauss_RSS$par[1],
                       t_gauss_RSS$par[2],
                       t_gauss_RSS$par[3],
                       t_gauss_RSS$par[4],
                       t_gauss_RSS$par[5])
  return(list(gauss_para, data))
}

fitGaussianModel2 <- function(data){
  # For each OTU, fits relative abundance as a gaussian function of temperature and rainfall
  # Designed to be iteratively run on subsets of OTUs
  # Args:
  #   data = data.frame containing OTU data
  
  data <- subset(data, nReads > 0) 
  data$relabund <- data$nReads / sum(data$nReads)
  rf_abund <- data.frame(y = data$relabund, x = data$rf_ann)
  t_abund <- data.frame(y = data$relabund,x = data$t_ann)
  rf_init_par <- c( ra_max = max(data$relabund),
                    ra_m = data$rf_ann[which.max(data$relabund)],
                    ra_v = sd(data$rf_ann[data$relabund>0]))
  t_init_par <- c( ra_max = max(data$relabund),
                   ra_m = data$t_ann[which.max(data$relabund)],
                   ra_v = sd(data$t_ann[data$relabund>0]))
  rf_gauss_RSS <- optim(par = rf_init_par, fn = RA_simple_RSS, method = "BFGS",
                        data = rf_abund)
  t_gauss_RSS <- optim(par = t_init_par, fn = RA_simple_RSS, method = "BFGS",
                       data = t_abund)
  names(rf_gauss_RSS$par) <- paste0("rf_", names(rf_gauss_RSS$par))
  names(t_gauss_RSS$par) <- paste0("t_", names(t_gauss_RSS$par))
  gauss_para <- data.frame(t(c(rf_gauss_RSS$par, t_gauss_RSS$par)))
  # Generate predicted relative abundances
  data$rf_predict <- RA_simple(x = rf_abund$x,
                        rf_gauss_RSS$par[1],
                        rf_gauss_RSS$par[2],
                        rf_gauss_RSS$par[3])
  data$t_predict <- RA_simple(x = t_abund$x,
                       t_gauss_RSS$par[1],
                       t_gauss_RSS$par[2],
                       t_gauss_RSS$par[3])
  gauss_para$site <- data$site[1]
  gauss_para$OTU_ID <- data$OTU_ID[1]
  return(list("gauss_para" = gauss_para, "data" = data))
}

calcQuantiles <- function(x){
  # Designed to take output from fitGaussian model
  # Calculates quantiles based on normal distribution parameters estimated using the fitGaussian model
  rf_q3 <- x$rf_ra_m + x$rf_ra_v * sqrt(2) * erfinv(2 * 0.80 - 1)
  rf_q1 <- x$rf_ra_m + x$rf_ra_v * sqrt(2) * erfinv(2 * 0.20 - 1)
  t_q3 <- x$t_ra_m + x$t_ra_v * sqrt(2) * erfinv(2 * 0.80 - 1)
  t_q1 <- x$t_ra_m + x$t_ra_v * sqrt(2) * erfinv(2 * 0.20 - 1)
  rf_iqr <- rf_q3 - rf_q1
  t_iqr <- t_q3 - t_q1
  data.frame(rf_q1, rf_q3, t_q1, t_q3)
}

elevationNullModel <- function(mat_x, mat_y, site.data, nreps, env.vars, site.var){
  # Compares species-level niche position and niche breadth between two sites
  # 
  # Args:
  #     mat_x: matrix, site (rows) by species (column) abundances for site x
  #     mat_y: matrix, site (rows) by species (column) abundances for site y
  #     site.data: data.frame, environment data by site
  #     env.vars: chr, column name of target environmental variables in siteData
  #     site.var: chr, column name of site IDs (this should contain the rownames in mat_x and mat_y)
  #     nreps: int, number of randomizations for the null model
  # Returns:
  #     data.frame containing observed differences in median position
  
  # env.vars = c("rf_ann", "t_ann")
  # site.var= "site.id"
  # site.data = siteData
  # mat_x = OTU_gauss_subset_t_laup_mat
  # mat_y = OTU_gauss_subset_t_stbk_mat
  # nreps = 10

  # Calculate observed differences
  obs <- compareNiche(mat_x = mat_x,
                      mat_y = mat_y,
                      site.data = site.data,
                      site.var = site.var,
                      env.vars = env.vars)
  rand_res <- list()
  for(i in 1:nreps){
   rand_mat_x <- randomizeMatrix(mat_x, null.model = "frequency")
   rand_mat_y <- randomizeMatrix(mat_y, null.model = "frequency")
   rand <- compareNiche(mat_x = rand_mat_x,
                        mat_y = rand_mat_y,
                        site.data = site.data,
                        site.var = site.var,
                        env.vars = env.vars)
   rand_res[[i]] <- rand
  }
  rand_res_comb <- do.call("rbind", rand_res)
  target.col <- do.call("paste0", expand.grid(env.vars, c("_range_diff", "_median_diff")))
  
  res_list <- list()
  for(i in 1:nrow(obs)){
    temp <- subset(rand_res_comb, ID == obs$ID[i])
    temp_list <- list()
    for(j in target.col){
      #j = "rf_ann_range_diff"
      temp_list[[ paste0(j, "_z") ]] <- (obs[i,j] - mean(temp[,j])) / sd(temp[,j])
      temp_list[[ paste0(j, "_p") ]] <- 1- (rank(c(obs[i,j], temp[,j]))[1] / (nreps+1))
    }
    res_list[[i]] <- data.frame(temp_list, "ID" = obs$ID[i])
  }
  res <- merge(obs, do.call("rbind", res_list), by= "ID")
  return(res)
}



compareNiche <- function(mat_x, mat_y, site.data, env.vars, site.var){
  x_niche <- calcNiche(mat_x, site.data = site.data, env.vars = env.vars, site.var = site.var)
  y_niche <- calcNiche(mat_y, site.data = site.data, env.vars = env.vars, site.var = site.var)
  
  splist <- unique(as.vector(x_niche$ID))
  spres <- list()
  for(i in 1:length(splist)){
    x <- subset(x_niche, ID == splist[i])
    y <- subset(y_niche, ID == splist[i])
    compare <- list()
    for(j in 1:length(env.vars)){
      range_x <- x[paste0(env.vars[j], "_range")]
      range_y <- y[paste0(env.vars[j], "_range")]
      range_diff <- range_x - range_y
      
      median_x <- x[paste0(env.vars[j], "_median")]
      median_y <- y[paste0(env.vars[j], "_median")]
      median_diff <- median_x - median_y
      
      names(range_diff) <- paste0(env.vars[j], "_range_diff")
      names(range_x) <- paste0(env.vars[j], "_range_x")
      names(range_y) <- paste0(env.vars[j], "_range_y")
      names(median_diff) <- paste0(env.vars[j], "_median_diff")
      names(median_x) <- paste0(env.vars[j], "_median_x")
      names(median_y) <- paste0(env.vars[j], "_median_y")
      
      compare[[j]] <- data.frame(range_x,
                                 range_y,
                                 range_diff,
                                 median_x,
                                 median_y,
                                 median_diff)
    }
    spres[[i]] <- do.call("cbind",compare)
    spres[[i]]$ID <- splist[i]
  }
  return(do.call("rbind",spres))
}


calcNiche <- function(mat, site.data, env.vars, site.var){
  # 
  # env.vars = c("rf_ann", "t_ann")
  # site.var= "site.id"
  # site.data = siteData
  # mat = OTU_gauss_subset_t_laup_mat
  df <- melt(mat, varnames = c(site.var, "ID"), value.name = "Abund")
  df_site <- merge(df, siteData[c(site.var, env.vars)], by = site.var)
  
  # Calculate quantiles for each species
  ddply(.data = df_site, .variables = .(ID), .fun = function(x) {
    niche_para <- list()
    for(i in 1:length(env.vars)){
      temp <- quantile(rep(x[[env.vars[i]]], times = x$Abund), type = 8, probs = c(0.25, 0.5, 0.75))
      names(temp) <- paste(env.vars[i], c("q1", "median", "q3"), sep = "_")
      niche_para[[i]] <- data.frame(t(temp))
      niche_para[[i]][[paste(env.vars[i], "range", sep = "_")]] <- temp[3] - temp[1]
    }
    return(do.call("cbind",niche_para))
  })
}

