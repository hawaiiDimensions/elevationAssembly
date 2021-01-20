# MAP OF TRANSECTS =====================
library(rgeos)
library(maptools)

## Directories ====================
main.dir <- "~/Dropbox/Projects/2017/hawaiiCommunityAssembly/elevationAssemblyHawaii/"
res.dir <- file.path(main.dir, "results")
fig.dir <- file.path(main.dir, "figures")
data.dir <- file.path(main.dir, "data")
siteData <- read.csv(file.path(data.dir, "clim.final.csv"), stringsAsFactors = FALSE)
siteData <- subset(siteData, island == "BigIsland" & elevation > 800)

# 
siteSP <- SpatialPointsDataFrame(coords = siteData[c("longitude","latitude")], proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"), data = siteData[c("longitude","latitude", "site")])


writeOGR(siteSP, dsn = file.path(main.dir, "gis"), layer = "transect", driver = "ESRI Shapefile")

# 
hillshade.dir <- c("~/Dropbox/maps/hawaii/MHI_digital_elevation_model_hillshade_GIS_data/")
hillshade <- raster(file.path(hillshade.dir, "Hawaii_Hillshade/Hawaii_Hillshade_45deg.img"))
# 
# siteSP <- SpatialPointsDataFrame(coords = siteData[c("longitude","latitude")], proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"), data = siteData[c("longitude","latitude", "site")])
# siteSP_reproj <- spTransform(siteSP, CRSobj = CRS(proj4string(hillshade)))
# siteSP_df <- cbind(siteSP_reproj@data["site"],siteSP_reproj@coords)
# 
# hillshade_fort <- fortify(hillshade, maxpixels = 100000) #197,920,112 actual pixels so this is a coarse graining
# 
# map_plot <- ggplot() +
#   geom_tile(aes(y = y, x = x, fill = Hawaii_Hillshade_45deg.1_BinValues), data = hillshade_fort) +
#   coord_fixed() +
#   scale_fill_gradient(low = grey(0.01), high = grey(1), na.value = "white") +
#   ggnewscale::new_scale_fill() +
#   geom_point(aes(y= latitude, x= longitude, fill = site, colour = site), data = siteSP_df,
#              pch = 21, size = 3.5) +
#   scale_colour_manual(values = c("white", "grey20")) +
#   geom_label(aes(x = 285000, y= 2210000, label = "Laupahoehoe"), size = 3) +
#   geom_label(aes(x = 270000, y= 2155000, label = "Stainback"), size = 3) +
#   scale_fill_manual(values = c("#003262", "#FDB515")) +
#   theme(plot.background = element_blank(),
#         panel.background = element_blank(),
#         axis.title = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         legend.position = "none")
# ggsave(file.path(fig.dir, "map.pdf"), map_plot, width = 5, height = 5)
# ggsave(file.path(fig.dir, "map.jpg"), map_plot, width = 5, height = 5)
# 
# # plot(hillshade, axes = FALSE, legend = FALSE, box = FALSE, col = gray.colors(20, start = 0, end = 1))
# # points(subset(siteSP_df, site == "Laupahoehoe")$latitude, subset(siteSP_df, site == "Laupahoehoe")$longitude, col = "#003262")
# # points(subset(siteSP_df, site == "Stainback")$latitude, subset(siteSP_df, site == "Stainback")$longitude, col = "#FDB515")
# 
# 

geologic <- readOGR("~/Dropbox/projects/2014/Hawaii/Maps/Haw_geo/Haw_St_geo_20070426_region.shp")
#geologic <- spTransform(geologic, CRS = CRS(proj4string(hillshade)))

mloa_comb <- unionSpatialPolygons(geologic[geologic$VOLCANO == "mloa",],
                                  rep(1, length(mloa)))
hual_comb <- unionSpatialPolygons(geologic[geologic$VOLCANO == "hual",],
                                  rep(2, length(hual)))
koha_comb <- unionSpatialPolygons(geologic[geologic$VOLCANO == "koha",],
                                  rep(3, length(koha)))
mkea_comb <- unionSpatialPolygons(geologic[geologic$VOLCANO == "mkea",],
                                  rep(4, length(mkea)))
kila_comb <- unionSpatialPolygons(geologic[geologic$VOLCANO == "kila",],
                                  rep(5, length(kila)))
spl = list(kila_comb, mkea_comb, koha_comb, hual_comb, mloa_comb)
joined = SpatialPolygons(Srl = lapply(spl, function(x){x@polygons[[1]]}))

joined2 = SpatialPolygonsDataFrame(Sr = joined, data = data.frame(ID = 1:5,
                                                                  volcano = c("mloa", "hual", "koha", "mkea", "kila")))

writeOGR(joined2, dsn = file.path(main.dir, "gis"), layer = "volcano", driver = "ESRI Shapefile")

