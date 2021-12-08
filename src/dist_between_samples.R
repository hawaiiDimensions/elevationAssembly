## ECOLOGICAL ASSEMBLY OF THE HAWAIIAN ARTHROPOD COMMUNITIES
# Find distance between sampling points for each transect
# Author: Jun Ying Lim

library(geosphere)

## IMPORT DATA ============
main.dir <- "~/Dropbox/projects/2017/hawaiiCommunityAssembly/elevationAssemblyHawaii/"
src.dir <- file.path(main.dir, "src")
data.dir <- file.path(main.dir, "data")
fig.dir <- file.path(main.dir, "figures")
res.dir <- file.path(main.dir, "results")

siteData <- read.csv(file.path(data.dir, "clim.final.csv"), stringsAsFactors = FALSE)

distGeo(p1 = c(1,0), p2 = c(0,0)) # units in meters

laup <- subset(siteData, site == "Laupahoehoe" & elevation > 800)

laup_dist <- matrix(NA, nrow = nrow(laup), ncol = nrow(laup))
for(i in 1:nrow(laup)){
  for(j in 1:nrow(laup)){
    laup_dist[i,j] <- distGeo(laup[c("longitude", "latitude")][i,], laup[c("longitude", "latitude")][j,] )
  }
}
min(laup_dist[laup_dist>0]) # minimum distance is 104 meters



stnbk <- subset(siteData, site == "Stainback" & elevation > 800)

stnbk_dist <- matrix(NA, nrow = nrow(stnbk), ncol = nrow(stnbk))
for(i in 1:nrow(stnbk)){
  for(j in 1:nrow(stnbk)){
    stnbk_dist[i,j] <- distGeo(stnbk[c("longitude", "latitude")][i,], stnbk[c("longitude", "latitude")][j,] )
  }
}
min(stnbk_dist[stnbk_dist>0]) # minimum distance is 152 meters
