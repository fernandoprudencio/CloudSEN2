library(tidyverse)
library(mapedit)
library(mapview)
library(raster)
library(tmap)
library(rgee)
library(sf)
source("utils.R")

data("World")
ee_Initialize("csaybar@gmail.com")

# 1. Create LandUse -------------------------------------------------------
ee_geom <- ee$Geometry$Rectangle(
  c(-179, -89, 179, 89),
  geodesic = FALSE,
  proj = "EPSG:4326"
)

ee_task <- ee_image_to_asset(
  image = ee_create_landuse(10000),
  description = "land_use_1km",
  region = ee_geom,
  assetId = "users/csaybar/world_land_use",
  overwrite = TRUE,
  maxPixels = 10**13
)

ee_task$start()
ee_monitoring()

land_use_metadata <- tibble(
  class = c("Barren", "Tropical Forest", "Temperated Forest",
            "Grass/Crop", "shrubland", "Snow/Ice/Lichen",
            "Urban", "Water", "Wetlands", "Ocean"),
  color = c("B4B4B4", "8DB400", "A0DC00", "F096FF",
            "FFBB22", "00E1FF", "FA0000", "0032C8",
            "0096A0", "0000FF"),
  value = 1:10
)

# 2. Select Potential points ---------------------------------------------------
# data("World")
# World <- World[!World$name == "Antarctica",]$geometry
# World <- st_transform(World, 4326)
# World <- World %>%
#   st_transform(4326) %>%
#   st_buffer(.00001) %>%
#   st_union() %>%
#   st_combine()
# sf_as_ee(x = World,
#          via = "getInfo_to_asset",
#          assetId = "users/csaybar/world")
# ee_world <- ee$FeatureCollection("users/csaybar/world")
# potential_points <- ee$FeatureCollection$randomPoints(ee_world$geometry(), 48000)
# task <- ee_table_to_asset(
#   collection = potential_points,
#   description = "potential_points",
#   assetId = "users/csaybar/potential_cloudsen2",
#   overwrite = TRUE
# )
# task$start()
# ee_monitoring()


# 3. Select Potential points by class ------------------------------------------

## baren
landuse_world <- ee$Image("users/csaybar/world_land_use")
barren <- landuse_world$eq(1)
map_base <- Map$addLayer(barren$updateMask(barren), list(min = 1, max = 10), name = "baren") +
  Map$addLayer(landuse_world, list(min = 0, max = 10, palette = land_use_metadata$color))
barren_points <- mapedit::editMap(ee_as_mapview(map_base))

cloudsen2 <- barren_points$drawn['geometry']
cloudsen2$type = "baren"



## Tropical Forest
tforest <- landuse_world$eq(2)
map_base <- Map$addLayer(tforest$updateMask(tforest), list(min = 1, max = 10), name = "tf") +
  Map$addLayer(landuse_world, list(min = 0, max = 10, palette = land_use_metadata$color))  +
  mapview(cloudsen2)
tforest_points <- mapedit::editMap(ee_as_mapview(map_base))
tf_db <- tforest_points$drawn['geometry']
tf_db$type = "Tropical Forest"
cloudsen2 <- rbind(cloudsen2, tf_db)



## Temparated Forest
tdforest <- landuse_world$eq(3)
map_base <- Map$addLayer(tdforest$updateMask(tdforest), list(min = 1, max = 10), name = "tdf") +
  Map$addLayer(landuse_world, list(min = 0, max = 10, palette = land_use_metadata$color))  +
  mapview(cloudsen2)
tdforest_points <- mapedit::editMap(ee_as_mapview(map_base))
tdf_db <- tdforest_points$drawn['geometry']
tdf_db$type = "Temparated Forest"
cloudsen2 <- rbind(cloudsen2, tdf_db)


## Grass/Crop
crop <- landuse_world$eq(4)
map_base <- Map$addLayer(crop$updateMask(crop), list(min = 1, max = 10), name = "crop") +
  Map$addLayer(landuse_world, list(min = 0, max = 10, palette = land_use_metadata$color))  +
  mapview(cloudsen2)
crop_points <- mapedit::editMap(ee_as_mapview(map_base))
crop_db <- crop_points$drawn['geometry']
crop_db$type = "Grass/Crop"
cloudsen2 <- rbind(cloudsen2, crop_db)


## shrubland
shrubland <- landuse_world$eq(5)
map_base <- Map$addLayer(shrubland$updateMask(shrubland), list(min = 1, max = 10), name = "shrubland") +
  Map$addLayer(landuse_world, list(min = 0, max = 10, palette = land_use_metadata$color))  +
  mapview(cloudsen2)
shrubland_points <- mapedit::editMap(ee_as_mapview(map_base))
shrubland_db <- shrubland_points$drawn['geometry']
shrubland_db$type = "Shrubland"
cloudsen2 <- rbind(cloudsen2, shrubland_db)


## Snow
snow <- landuse_world$eq(6)
map_base <- Map$addLayer(snow$updateMask(snow), list(min = 1, max = 10), name = "snow") +
  Map$addLayer(landuse_world, list(min = 0, max = 10, palette = land_use_metadata$color))  +
  mapview(cloudsen2)
snow_points <- mapedit::editMap(ee_as_mapview(map_base))
snow_db <- snow_points$drawn['geometry']
snow_db$type = "Snow"
cloudsen2 <- rbind(cloudsen2, snow_db)

## Urban
urban <- landuse_world$eq(7)
map_base <- Map$addLayer(urban$updateMask(urban), list(min = 1, max = 10), name = "urban") +
  Map$addLayer(landuse_world, list(min = 0, max = 10, palette = land_use_metadata$color))  +
  mapview(cloudsen2)

urban_points <- mapedit::editMap(ee_as_mapview(map_base))
urban_db <- urban_points$drawn['geometry']
urban_db$type = "Urban"
cloudsen2 <- rbind(cloudsen2, urban_db)


## Water
water <- landuse_world$eq(8)
map_base <- Map$addLayer(water$updateMask(water), list(min = 1, max = 10), name = "water") +
  Map$addLayer(landuse_world, list(min = 0, max = 10, palette = land_use_metadata$color))  +
  mapview(cloudsen2)

water_points <- mapedit::editMap(ee_as_mapview(map_base))
water_db <- water_points$drawn['geometry']
water_db$type = "Water"
cloudsen2 <- rbind(cloudsen2, water_db)


## Wetlands
wetlands <- landuse_world$eq(9)
map_base <- Map$addLayer(wetlands$updateMask(wetlands), list(min = 1, max = 10), name = "wetlands") +
  Map$addLayer(landuse_world, list(min = 0, max = 10, palette = land_use_metadata$color))  +
  mapview(cloudsen2)

wetlands_points <- mapedit::editMap(ee_as_mapview(map_base))
wetlands_db <- wetlands_points$drawn['geometry']
wetlands_db$type = "Wetlands"
cloudsen2 <- rbind(cloudsen2, wetlands_db)

write_sf(cloudsen2, "data/cloudsen2.geojson")
