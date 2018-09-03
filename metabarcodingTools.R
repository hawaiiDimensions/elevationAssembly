

collapseOTU <- function(x, ref, collapse.by, to.collapse, subset.by){
  # Arguments:
  #     x = matrix, e.g., zotu table to be collapsed; sites as columns, otus as rows
  #     ref = data.frame, reference table by which otu table
  #     collapse.by = string, specifies column in ref data.frame IDs to collapse by in input matrix, x (e.g., species)
  #     to.collapse = string, specifies column in ref data.frame IDs to be collapsed in input matrix, x (e.g., zOTU)
  #     subset.by = string, specifies column in ref data.frame IDs to be subset. If specified, output will be a list of matrices
  # Returns:
  #    matrix, or list of matrices, that have been collapsed by specified parameters
  
  if(missing(collapse.by)){
    stop("Missing argument: 'collapse.by' ")
  } else if(!collapse.by %in% names(ref)){
    stop("collapse.by not found in reference dataframe, ref")
  }
  if(missing(to.collapse)){
    stop("Missing argument: 'to.collapse' ")
  } else if(!to.collapse %in% names(ref)) {
    stop("to.collapse not found in reference dataframe, ref")
  }
  if(!missing(subset.by)){
    if(!subset.by %in% names(ref)){
      stop("subset.by not found in reference dataframe, ref")
    }
  }
  
  collapseByList <- unique(ref[[collapse.by]])
  resList <- list()
  for(i in 1:length(collapseList)){
    tocollapseList <- ref[ref[collapse.by]==collapseByList[i],][[to.collapse]]
    if(length(tocollapseList) == 1){
      resList[[i]] <- x[rownames(x) %in% tocollapseList,]
    } else {
      resList[[i]] <- colSums(x[rownames(x) %in% tocollapseList,])
    }
  }
  res <- do.call("rbind", resList)
  rownames(res) <- collapseByList
  
  subset.by = "V4"
  if(missing(subset.by)){
    return(res)
  } else{
    subsetres <- list()
    subsetList <- unique(as.vector(ref[[subset.by]]))
    for(j in 1:length(subsetList)){
      subsetCollapsed <- unique(ref[ref[subset.by]==subsetList[j],][[collapse.by]])
      subsetres[[subsetList[j]]] <- res[subsetCollapsed,]
    }
    return(subsetres)
  }
}



rarefyOTUbySpecimenCount <- function(df, counts, readsPerIndividual, nReps){
  # Rarefies category by specimen counts. This function is designed to be iterated using `ddply`
  # Args:
  #   counts: dataframe containing specimen counts
  #   readsPerIndividual: number of reads per individual you would like to rarefy with
  #   nReps: number of randomizations
  # Returns:
  #   
  
  # Test parameters
  # readsPerIndividual = 10; nReps = 10
  # df = subset(dataARFmelt, SizeCategory == "0to2" & Site_ID == "680")
  # counts = specimenCounts
  
  # Extract the number of specimens for site and size category
  nSpecimen <- subset(counts, Site_ID == unique(df$Site_ID) &
                        SizeCategory == unique(df$SizeCategory))$Count
  
  rawReadAbund <- acast(df, X.OTU.ID~Site_SizeCategory, value.var = "nReads")
  #rarefiedReadAbund <- list()
  #for(i in 1:nReps){
  #rarefiedReadAbund[[i]] <- rrarefy(rawReadAbund, sample = nSpecimen*readsPerIndividual) 
  #}
  rarefiedReadAbund <- rrarefy(rawReadAbund, sample = nSpecimen*readsPerIndividual) 
  temp <- melt(rarefiedReadAbund,
               value.name = "rarefiedReadAbund",
               varnames = c("Site_SizeCategory", "X.OTU.ID"))
  
  res <- merge(df, temp, by = c("X.OTU.ID", "Site_SizeCategory"))
  return(res)
}


matchDist <- function(xdis, ydis, by = "x"){
  # Matches and orders rows and columns
  #xdis = testData_betadist
  #ydis = climDist
  xdist_mat <- as.matrix(xdis)
  ydist_mat <- as.matrix(ydis)
  if(by == "x"){
    ydist_mat <- ydist_mat[rownames(xdist_mat), colnames(xdist_mat)]
    return(as.dist(ydist_mat))
  } else {
    xdist_mat <- xdist_mat[rownames(ydist_mat), colnames(ydist_mat)]
    return(as.dist(xdist_mat))
  }
}


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



