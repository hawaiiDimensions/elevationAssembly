

distmatrix2df <- function(matrix){
    # Converts a matrix into a dataframe with no duplicate pairwise values
    data.frame( t(combn(names(matrix),2)), dist=t(matrix)[lower.tri(matrix)] )  
}

mergeSiteData <- function(data, siteData, targetCol){
  # Merge pairwise
  
  #data <- betaAll; siteData <- siteData; targetCol <- c("elevation")
  
  site1match <- match(data$X1, siteData$site.id)
  site2match <- match(data$X2, siteData$site.id)
  data$ystart <- siteData[site1match,]$latitude
  data$xstart <- siteData[site1match,]$longitude
  data$yend <- siteData[site2match,]$latitude
  data$xend <- siteData[site2match,]$longitude
  data$site1 <- siteData[site1match,]$site
  data$site2 <- siteData[site2match,]$site
  data$geogDist <- sapply(X = 1:dim(data)[1], FUN = function(X){distVincentyEllipsoid(c(data$xstart[X], data$ystart[X]), c(data$xend[X], data$yend[X]))})
  
  for(i in targetCol){
    data[paste0(i, ".diff")] <- abs(siteData[site1match, i] - siteData[site2match, i])
  }
  return(data)
}


betaPairList <- function(comm){
  # Compute the nestedness and turnover components of beta-diversity
  betaPair <- beta.pair(comm)
  otuBetaTurnover <- as.data.frame(as.matrix(betaPair$beta.sim))
  otuBetaNest <- as.data.frame(as.matrix(betaPair$beta.sne))
  otuBetaTotal <- as.data.frame(as.matrix(betaPair$beta.sor))
  
  betaTurnoverDF <- distmatrix2df(otuBetaTurnover)
  betaNestDF <- distmatrix2df(otuBetaNest)
  betaTotalDF <- distmatrix2df(otuBetaTotal)
  
  names(betaTurnoverDF)[names(betaTurnoverDF) == "dist"] <- "turnover"
  names(betaNestDF)[names(betaNestDF) == "dist"] <- "nestedness"
  names(betaTotalDF)[names(betaTotalDF) == "dist"] <- "total"
  betaAll <- Reduce(list(betaTurnoverDF, betaNestDF, betaTotalDF), f = merge)
  return(betaAll)
}

betaPairNull <- function(comm, null, nsim){
  # Compute observed pairwise beta diversity
  #comm <- steinOtuPA; null <- "r1"; nsim <- 100
  obsBetaPair <- betaPairList(comm)
  
  # Compute
  nullSim <- simulate(nullmodel(x = comm, method = null), nsim = nsim)
  nullBetaPair <- alply(.data = nullSim, .margins = 3, .fun = betaPairList)
  
  turnoverList <- lapply(nullBetaPair, FUN = function(x){ x$turnover})
  nestednessList <- lapply(nullBetaPair, FUN = function(x){ x$nestedness})
  totalList <- lapply(nullBetaPair, FUN = function(x){ x$total})
  
  turnoverMatrix <- matrix(unlist(turnoverList), ncol = nsim, byrow = FALSE)
  nestednessMatrix <- matrix(unlist(nestednessList), ncol = nsim, byrow = FALSE)
  
  obsBetaPair$turnoverNullMean <- rowMeans(turnoverMatrix)
  obsBetaPair$turnoverNullSD <- apply(MARGIN = 1, X = turnoverMatrix, FUN = sd)
  obsBetaPair$turnoverObsRank <- sapply(1:dim(obsBetaPair)[1], function(X){rank(c(obsBetaPair$turnover[X], turnoverMatrix[X,]), ties.method = "last")[1]})
  obsBetaPair$turnoverPvalue <- obsBetaPair$turnoverObsRank / (nsim + 1)
  
  obsBetaPair$nestednessNullMean <- rowMeans(nestednessMatrix)
  obsBetaPair$nestednessNullSD <- apply(MARGIN = 1, X = nestednessMatrix, FUN = sd)
  obsBetaPair$nestednessObsRank <- sapply(1:dim(obsBetaPair)[1], function(X){rank(c(obsBetaPair$nestedness[X], nestednessMatrix[X,]), ties.method = "last")[1]})
  obsBetaPair$nestednessPvalue <- obsBetaPair$nestednessObsRank / (nsim + 1)
  
  return(obsBetaPair)
}



