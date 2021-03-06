#' Pass a COPERNICUS 100 m Land use ee$Image to a multiband ee$Image where
#' each band is a land-use class (<0-1>, spatial one-hot encoding).
#'
#' It is an auxiliary data that helps to determinate manually the points in
#' CLOUDSEN1000
ee_create_landuse <- function(resolution = 1000) {
  # Read land use product
  ee_product <- "COPERNICUS/Landcover/100m/Proba-V/Global"
  lc_100 <- ee$ImageCollection(ee_product) %>%
    ee$ImageCollection$select("discrete_classification") %>%
    ee$ImageCollection$first()

  #barren
  lc_class_barren <- 60
  lc_barren <- lc_100$eq(lc_class_barren)$multiply(1)

  # tropical forest
  lc_class_tropicalforest <- 112
  lc_tropicalforest <- lc_100$eq(lc_class_tropicalforest)$multiply(2)

  # temperated forest
  lc_class_temparateforest <- c(111,113, 115)
  lc_temparateforest <- lc_100 %>%
    ee$Image$eq(lc_class_temparateforest) %>%
    ee$Image$reduce(ee$Reducer$anyNonZero()) %>%
    ee$Image$multiply(3)

  # grass
  lc_class_grass <- c(30, 40)
  lc_grass <- lc_100 %>%
    ee$Image$eq(lc_class_grass) %>%
    ee$Image$reduce(ee$Reducer$anyNonZero()) %>%
    ee$Image$multiply(4)

  # shrubland
  lc_class_shrubland <- 20
  lc_shrubland <- lc_100$eq(lc_class_shrubland)$multiply(5)

  # snow
  lc_class_snow <- c(100, 70)
  lc_snow <- lc_100 %>%
    ee$Image$eq(lc_class_snow) %>%
    ee$Image$reduce(ee$Reducer$anyNonZero()) %>%
    ee$Image$multiply(6)


  # urban
  lc_class_urban <- 50
  lc_urban <- lc_100$eq(lc_class_urban)$multiply(7)

  # water
  lc_class_water <- 80
  lc_water <- lc_100$eq(lc_class_water)$multiply(8)


  # wetlands
  lc_class_wetlands <- 90
  lc_wetlands <- lc_100$eq(lc_class_wetlands)$multiply(9)

  # ocean
  lc_class_ocean <- 200
  lc_ocean <- lc_100$eq(lc_class_ocean)$multiply(10)

  # bringing all together!
  new_landuse <- ee$Image$cat(
      lc_urban,
      lc_water,
      lc_snow,
      lc_grass,
      lc_barren,
      lc_shrubland,
      lc_tropicalforest,
      lc_temparateforest,
      lc_wetlands,
      lc_wetlands
  )$reduce(ee$Reducer$max())

  # Reduce Land Use to (1km)
  new_landuse %>%
    ee$Image$reduceResolution(
      reducer = ee$Reducer$mode(),
      bestEffort = TRUE
    ) %>%
  ee$Image$reproject(
    crs = "EPSG:4326",
    scale = resolution
  )
}


#' Create 1000x1000 SENTINEL-2 (level-1C) chips to feed a DL model
dataset_creator_thumbnail <- function(cloudsen2_row,
                                      bands = "B.*|probability|SCL",
                                      scale = 100,
                                      img_dim = c(1000, 1000),
                                      data_range = c("2019-01-01", "2020-07-31"),
                                      output = "results/") {
  # 1. Create a point which represent the center of the chip
  point <- ee$Geometry$Point(cloudsen2_row$geometry[[1]])

  # 2. Create a S2 ImageCollection filtering by space and time.
  s2Sr <- ee$ImageCollection("COPERNICUS/S2_SR") %>%
    ee$ImageCollection$filterBounds(point) %>%
    ee$ImageCollection$filterDate(data_range[1], data_range[2])

  # 2. Create a S2 ImageCollection filtering by space and time.
  s2_level1c <- ee$ImageCollection("COPERNICUS/S2") %>%
    ee$ImageCollection$filterBounds(point) %>%
    ee$ImageCollection$filterDate(data_range[1], data_range[2])

  # 3. Create a CRS
  crs_kernel <- s2Sr$first()$select(0)$projection()$getInfo()$crs
  point_utm <- st_transform(cloudsen2_row$geometry[1], crs_kernel)
  ee_point <- ee$Geometry$Point(point_utm[[1]], proj = crs_kernel)

  # 4. Create a S2_CLOUD_PROBABILITY ImageCollection filtering by space and time.
  s2Clouds <- ee$ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY') %>%
    ee$ImageCollection$filterBounds(ee_point)%>%
    ee$ImageCollection$filterDate(data_range[1], data_range[2])

  # 5. Merge S2 and S2_CLOUD_PROBABILITY
  # 'add_S2cloud' is a function which convert the 'cloud' property to a band.
  s2SrWithCloudMask_1 <- ee$Join$saveFirst('cloud_mask')$apply(
    primary = s2Sr,
    secondary = s2Clouds,
    condition = ee$Filter$equals(leftField = 'system:index', rightField = 'system:index')
  ) %>% ee$ImageCollection() %>% ee$ImageCollection$map(add_S2cloud)

  s2SrWithCloudMask <- ee$Join$saveFirst('cloud_mask')$apply(
    primary = s2_level1c,
    secondary = s2SrWithCloudMask_1,
    condition = ee$Filter$equals(leftField = 'system:index', rightField = 'system:index')
  ) %>% ee$ImageCollection() %>% ee$ImageCollection$map(add_S2cloud2)

  # Create a function to estimate the mean cloud propability
  # cloud_prob_fn <- cloud_fun_creator_s2cloudness(scale = 1000)
  cloud_prob_fn <- cloud_fun_creator_sen2cor(scale = 10,
                                             geom = ee_point$buffer(10*512)$bounds())
  # 6. Apply 'cloud_fun_creator' to s2SrWithCloudMask
  # cloud probabibily with a convolves and dilation
  s2SrWithCloudMask %>%
    ee$ImageCollection$map(cloud_prob_fn) %>%
    ee$ImageCollection$aggregate_array("cprob") %>%
    ee$List$getInfo() %>%
    unlist -> cloud_percentage

  # 7. Identify the image with the cloud coverage closer to the desire cloud percentage
  if (!all(cloud_percentage == 0)) {
    cloud_percentage[cloud_percentage == 0] <- -999
  }
  ideal_position <- which.min(abs(cloud_percentage*100 - cloudsen2_row$potential_probability))[1]

  # 8. Get Image ID
  s2_chip_id <- ee_get(s2SrWithCloudMask, ideal_position - 1)$first()$
    get("system:id")$getInfo() %>% basename()

  # 9. Create a folder to save resutls
  dir_id <- sprintf("%s/%s",output, s2_chip_id)
  dir_results <- sprintf("%s/results", dir_id, showWarnings = FALSE)
  dir.create(dir_id, showWarnings = FALSE)
  dir.create(dir_results, showWarnings = FALSE)

  # 10. Select an image considering a position
  s2_img_data <- ee_get(s2SrWithCloudMask, ideal_position - 1) %>%
    ee$ImageCollection$first()

  ## 11. SAVE image values
  s2_cdi <- ee$Algorithms$Sentinel2$CDI(s2_img_data)
  s2_ref <- s2_img_data %>%
    ee$Image$select(bands) %>%
    ee$Image$addBands(s2_cdi) %>%
    ee$Image$reproject(crs = crs_kernel, scale = 100)

  # Classify bands according our method
  # Cloud mask according to Zupanc et al. 2019
  boxcar1 <- ee$Kernel$square(radius = 4, units = 'pixels')
  boxcar2 <- ee$Kernel$square(radius = 2, units = 'pixels')
  s2cloudness_prob <- s2_ref$select("probability") %>%
    ee$Image$convolve(boxcar1) %>%
    ee$Image$focal_max(kernel = boxcar2, iterations = 1) %>%
    ee$Image$gte(70) %>%
    ee$Image$rename("s2cloudness_mask")

  # Cloud mask according to our sen2cor
  sen2cor_scl <- sen2cor_reclass(s2_ref$select("SCL"))$rename("sen2cor_mask")
  s2_ref_f <- s2_ref$addBands(sen2cor_scl)$addBands(s2cloudness_prob) %>%
    ee$Image$float()

  ref <- ee_as_raster(
    image = s2_ref_f,
    scale = 100,
    via = "drive",
    crs = crs_kernel
  )
  band_names <- names(ref)
  ref[is.na(ref)] <- 0
  names(ref) <- band_names

  # Save bands
  ## group1 - predictors
  bands_group1 <- names(ref)[grepl("B.*|cdi", names(ref))]
  ref_group1 <- ref[[bands_group1]]/10000
  ref_group1[["cdi"]] <- ref_group1[["cdi"]]*10000
  writeRaster(ref_group1,  sprintf("%s/predictors.tif", dir_results))

  ## group2 - cloud values
  bands_group2 <- names(ref)[grepl("SCL|probability", names(ref))]
  sen2cor_scl <- ref[[bands_group2[1]]]
  sen2cloudness_cprop <- ref[[bands_group2[2]]]
  writeRaster(sen2cor_scl,
              sprintf("%s/sen2cor_SCL.tif", dir_results))
  writeRaster(sen2cloudness_cprop,
              sprintf("%s/sen2cloudness_prop.tif", dir_results))

  ## group3 - cloud mask
  writeRaster(ref[["sen2cor_mask"]],
              sprintf("%s/sen2cor_mask.tif", dir_results))
  writeRaster(ref[["s2cloudness_mask"]],
              sprintf("%s/s2cloudness_mask.tif", dir_results))

  ## 12. Obtain metadata (% cloud, azimuth)
  s2_ref_cloud <- cloud_percentage[ideal_position] * 100
  s2_ref_azimuth_angle <- s2_img_data$get("MEAN_SOLAR_AZIMUTH_ANGLE")$getInfo()
  s2_ref_zenith_angle <- s2_img_data$get("MEAN_SOLAR_ZENITH_ANGLE")$getInfo()

  ## Calculating angle shadow
  azimuth_radians <- (s2_ref_azimuth_angle + 180)*pi/180
  zenith_radians <- (s2_ref_azimuth_angle * pi)/180
  shadowCastedDistance <- tan(zenith_radians)*100
  x <- sin(azimuth_radians) * shadowCastedDistance * -1
  y <- cos(azimuth_radians) * shadowCastedDistance * -1
  angle_shadow <- atan2(y,x)*(180/pi)

  title <- sprintf("Cloud:%s, dir_shadow:%s",
                   round(s2_ref_cloud,2),
                   round(angle_shadow,2))
  message(title)

  # 13. Save the REF image (interactive)
  rmap01 <- viewRGB(
    x = ref[[c("B2","B3","B4")]],
    quantiles = c(0.02, 0.98),
    na.color = "#BEBEBE00",
    layer.name = sprintf("%s",title)
  )

  # B10
  rmap02 <- mapview(
    x = ref[[c("B10")]],
    map =  rmap01@map,
    quantiles = c(0, 1),
    layer.name = "B10"
  )

  # CDI
  rmap03 <- mapview(
    x = ref[[c("cdi")]],
    map =  rmap02@map,
    quantiles = c(0, 1),
    layer.name = "CDI"
  )

  # cloud probability
  rmap04 <- mapview(
    x = ref[[c("sen2cor_mask")]],
    map =  rmap03@map,
    col.regions = c("#007B00","#FFB000","#606060"),
    at = c(-1, 0, 1, 2),
    layer.name = "sen2cor_mask",
    na.color = "#BEBEBE00",
    legend = TRUE
  )
  rmap05 <- mapview(
    x = ref[[c("s2cloudness_mask")]],
    map =  rmap04@map,
    col.regions = c("#007B00","#FFB000"),
    at = c(-1, 0, 1),
    layer.name = "sen2cloudness_mask",
    na.color = "#BEBEBE00",
    legend = TRUE
  )

  mapshot(rmap05, url = sprintf("%s/map.html",normalizePath(dir_id)))
  unlink(sprintf("%s/map_files",normalizePath(dir_id)), recursive = TRUE)

  # 14. Save static values
  # RGB -REAL COLOR
  cloudsen2_plot(ref_group1[[c("B4","B3","B2")]], sprintf("%s/RGB.png",dir_id))

  # NRG - FALSE COLOR
  cloudsen2_plot(ref_group1[[c("B8","B4","B3")]], sprintf("%s/NRG.png",dir_id))

  # NDVI
  final_stack_ndvi <- (ref_group1[["B8"]] - ref_group1[["B4"]])/(ref_group1[["B8"]] + ref_group1[["B4"]])
  cloudsen2_plot(
    cloud_brick = final_stack_ndvi,
    linStretchVec = FALSE,
    output = sprintf("%s/NDVI.png",dir_id),
    pal = viridis_pal(100),
    limits = c(-0.1, 0.4)
  )

  # NDSI
  final_stack_ndsi <- (ref_group1[["B3"]] - ref_group1[["B11"]])/(ref_group1[["B3"]] + ref_group1[["B11"]])
  cloudsen2_plot(
    cloud_brick = final_stack_ndsi,
    output = sprintf("%s/NDSI.png",dir_id),
    linStretchVec = FALSE,
    pal = viridis_pal(100),
    limits = c(0.4, 1)
  )

  # CDI
  cdi <- ref_group1[["cdi"]]
  cloudsen2_plot(
    cloud_brick = cdi,
    output = sprintf("%s/CDI.png", dir_id),
    linStretchVec = FALSE,
    pal = viridis_pal(100),
    limits = c(-1, -0.4)
  )

  # CIRRUS
  cirrus <- ref_group1[["B10"]]
  cloudsen2_plot(
    cloud_brick = cirrus,
    output = sprintf("%s/CIRRUS.png",dir_id),
    linStretchVec = FALSE,
    pal = viridis_pal(100),
    limits = c(0.005, 0.01)
  )

  # POTENTIAL CLOUD
  pcloud <- ref[["probability"]]/100
  cloudsen2_plot(
    cloud_brick = pcloud,
    output = sprintf("%s/s2cloudness_prob.png",dir_id),
    linStretchVec = FALSE,
    pal = viridis_pal(100),
    limits = c(0.4, 1)
  )

  # Target
  cloudsen2_plot(
    cloud_brick = ref[["sen2cor_mask"]],
    output = sprintf("%s/sen2cor_mask.png",dir_id),
    linStretchVec = FALSE,
    pal = c("#007B00","#00FFFF", "#7D7D7D"),
    limits = c(0, 2)
  )
  cloudsen2_plot(
    cloud_brick = ref[["s2cloudness_mask"]],
    output = sprintf("%s/s2cloudness_mask.png",dir_id),
    linStretchVec = FALSE,
    pal = c("#007B00","#00FFFF"),
    limits = c(0, 1)
  )
  tile_center <- cloudsen2_row$geometry[[1]]
  list(
    landcover = cloudsen2_row$type,
    bands = names(ref),
    source = s2_chip_id,
    pointref_x = tile_center[1],
    pointref_y = tile_center[2],
    avg_cloud_percentage = s2_ref_cloud,
    mean_solar_azimuth_angle = s2_ref_azimuth_angle
  ) -> chip_metadata

  write_json(chip_metadata, sprintf("%s/metadata.json", dir_id))
  new_dir_id <- sprintf("%s/%s_CC_%s_%s", output, index, round(s2_ref_cloud,2), s2_chip_id)
  file.rename(dir_id, new_dir_id)
  invisible(TRUE)
}

#' Create 511x511 SENTINEL-2 (level-1C) chips to feed a DL model
dataset_creator_chips <- function(cloudsen2_row,
                                  bands = "B.*|probability|SCL",
                                  kernel_size = c(255, 255),
                                  data_range = c("2019-01-01", "2020-07-31"),
                                  output = "results/") {
  # 1. Create a point which represent the center of the chip
  point <- ee$Geometry$Point(cloudsen2_row$geometry[[1]])

  # 2. Create a S2 ImageCollection with all the necessary images.
  list_data <- tsen2_dataset_creator(point, data_range)
  s2SrWithCloudMask <- list_data$ic
  crs_kernel <- list_data$metadata$crs
  point_utm <- list_data$metadata$local_point
  ee_point <- list_data$metadata$ee_point

  # 3. Create a function to estimate the mean cloud propability
  ee_new_kernel <- ee_point$buffer(10*255)$bounds()
  cloud_prob_fn <- cloud_fun_creator_sen2cor(scale = 10,
                                             geom = ee_new_kernel)

  # 4. Apply 'cloud_prob_fn' to s2SrWithCloudMask
  s2SrWithCloudMask %>%
    ee$ImageCollection$map(cloud_prob_fn) %>%
    ee$ImageCollection$aggregate_array("cprob") %>%
    ee$List$getInfo() %>%
    unlist -> cloud_percentage

  # 5. Identify the image with the cloud coverage closer to the desire
  # cloud percentage.
  if (!all(cloud_percentage == 0)) {
    cloud_percentage[cloud_percentage == 0] <- -999
  }
  ideal_position <- which.min(abs(cloud_percentage*100 - cloudsen2_row$potential_probability))[1]
  # cloud_percentage[ideal_position]*100
  # cloudsen2_row$potential_probability
  # 6. Pick the Image
  s2_img_data <- ee_get(s2SrWithCloudMask, ideal_position - 1) %>%
    ee$ImageCollection$first()
  s2_chip_id <- s2_img_data$get("system:id")$getInfo() %>% basename()

  # 7. Create folders to save results
  dir_id <- sprintf("%s/%s",output, s2_chip_id)
  dir_results <- sprintf("%s/results", dir_id, showWarnings = FALSE)
  dir.create(dir_id, showWarnings = FALSE)
  dir.create(dir_results, showWarnings = FALSE)

  # 8. Save image values
  s2_cdi <- ee$Algorithms$Sentinel2$CDI(s2_img_data)
  s2_ref <- s2_img_data %>%
    ee$Image$select(bands) %>%
    ee$Image$addBands(s2_cdi)
  # ee$Image$reproject(crs = crs_kernel, scale = 10)

  # 9. Classify bands according our method
  ## Cloud mask according to Zupanc et al. 2019
  boxcar1 <- ee$Kernel$square(radius = 4, units = 'pixels')
  boxcar2 <- ee$Kernel$square(radius = 2, units = 'pixels')
  s2cloudness_prob <- s2_ref$select("probability") %>%
    ee$Image$convolve(boxcar1) %>%
    ee$Image$focal_max(kernel = boxcar2, iterations = 1) %>%
    ee$Image$gte(70) %>%
    ee$Image$rename("s2cloudness_mask")

  ## Cloud mask according to our sen2cor
  sen2cor_scl <- sen2cor_reclass(s2_ref$select("SCL"))$rename("sen2cor_mask")
  s2_ref_f <- s2_ref$addBands(sen2cor_scl)$addBands(s2cloudness_prob) %>%
    ee$Image$float()
  # Map$centerObject(ee_new_kernel)
  # Map$addLayer(s2_ref_f$select("sen2cor_mask")) +  Map$addLayer(ee_new_kernel, name = "map_3")
  # demo1 <- ee_as_raster(s2_ref_f$select("sen2cor_mask"), ee_new_kernel,via = "drive")

  # 10. Create a reference image
  ## It helps us to understand better the context of the chip
  ## From ee to local
  ref <- ee_as_raster(
    image = s2_ref_f$select(c("B4","B3","B2")),
    region = ee_point$buffer(kernel_size[1]*40)$bounds(),
    scale = 50,
    via = "drive"
  )

  ## Obtain metadata (% cloud, azimuth)
  s2_ref_cloud <- cloud_percentage[ideal_position] * 100
  s2_ref_azimuth_angle <- s2_img_data$get("MEAN_SOLAR_AZIMUTH_ANGLE")$getInfo()
  s2_ref_zenith_angle <- s2_img_data$get("MEAN_SOLAR_ZENITH_ANGLE")$getInfo()

  ### Calculating angle shadow
  azimuth_radians <- (s2_ref_azimuth_angle + 180)*pi/180
  zenith_radians <- (s2_ref_azimuth_angle * pi)/180
  shadowCastedDistance <- tan(zenith_radians)*100
  x <- sin(azimuth_radians) * shadowCastedDistance * -1
  y <- cos(azimuth_radians) * shadowCastedDistance * -1
  angle_shadow <- atan2(y,x)*(180/pi)

  # Save the REF image (static)
  kroi <- ee_as_sf(x = ee_point$buffer(10*256)$bounds())$geometry %>%
    st_cast("MULTILINESTRING")
  kroi_utm <- st_transform(kroi, crs_kernel)

  png(sprintf("%s/ref_RGB.png", dir_id), 800, 600)
  plotRGB(ref/1000*255, stretch='hist')
  plot(kroi_utm, lwd=1, col = "red", add = TRUE, lty=2)
  dev.off()

  # Save the REF image (interactive)
  Rmap <- viewRGB(
    x = ref[[c(3,2,1)]],
    quantiles = c(0.02, 0.98),
    na.color = "#BEBEBE00",
    layer.name = "Reference"
  )

  ## Create a TARGET image
  # Retrieve data only on the kernel (511x511)
  s2_img_array <- s2_ref_f %>%
    ee$Image$addBands(ee$Image$pixelLonLat()) %>%
    ee$Image$reproject(crs = "EPSG:4326", scale = 10) %>%
    ee$Image$neighborhoodToArray(
      kernel = ee$Kernel$rectangle(kernel_size[1], kernel_size[2], "pixels")
    ) %>%
    ee$Image$sampleRegions(ee$FeatureCollection(point)) %>%
    ee$FeatureCollection$getInfo()

  # Convert data from list to data_frame
  message("Processing image: ", s2_chip_id)
  band_names <- names(s2_img_array$features[[1]]$properties)
  extract_fn <- function(x) as.numeric(unlist(s2_img_array$features[[1]]$properties[x]))
  image_as_df <- do.call(cbind,lapply(band_names, extract_fn))
  colnames(image_as_df) <- band_names
  image_as_tibble <- as_tibble(image_as_df)

  # Convert data from data_frame to stack
  coordinates(image_as_tibble) <- ~longitude+latitude
  sf_to_stack <- function(x) rasterFromXYZ(image_as_tibble[x])
  final_stack <- stack(lapply(names(image_as_tibble), sf_to_stack))
  crs(final_stack) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  names(final_stack) <- band_names[!band_names %in% c("latitude", "longitude")]

  B_bands <- grepl("B", names(final_stack))

  bstack <- final_stack[[names(final_stack)[B_bands]]]*0.0001
  ostack <- final_stack[[names(final_stack)[!B_bands]]]
  final_stack <- stack(bstack, ostack)

  # Save bands
  ## group1 - predictors
  bands_group1 <- names(final_stack)[grepl("B.*|cdi", names(final_stack))]
  ref_group1 <- final_stack[[bands_group1]]
  writeRaster(ref_group1,  sprintf("%s/predictors.tif", dir_results))

  ## group2 - cloud values
  bands_group2 <- names(final_stack)[grepl("SCL|probability", names(final_stack))]
  sen2cor_scl <- final_stack[[bands_group2[1]]]
  sen2cloudness_cprop <- final_stack[[bands_group2[2]]]
  writeRaster(sen2cor_scl,
              sprintf("%s/sen2cor_SCL.tif", dir_results))
  writeRaster(sen2cloudness_cprop,
              sprintf("%s/sen2cloudness_prop.tif", dir_results))

  ## group3 - cloud mask
  writeRaster(final_stack[["sen2cor_mask"]],
              sprintf("%s/sen2cor_mask.tif", dir_results))
  writeRaster(final_stack[["s2cloudness_mask"]],
              sprintf("%s/s2cloudness_mask.tif", dir_results))

  rmap_target <- viewRGB(
    x = final_stack[[c(5,6,7)]],
    map =  Rmap@map,
    quantiles = c(0.02, 0.98),
    layer.name = "target_RGB"
  )

  # B10
  rmap2_target <- mapview(
    x = final_stack[[c("B10")]],
    map =  rmap_target@map,
    quantiles = c(0, 1),
    layer.name = "B10"
  )

  # CDI
  rmap3_target <- mapview(
    x = final_stack[[c("cdi")]],
    map =  rmap2_target@map,
    quantiles = c(0, 1),
    layer.name = "CDI"
  )

  # cloud probability
  rmap4_target <- mapview(
    x = final_stack[[c("probability")]],
    map =  rmap3_target@map,
    quantiles = c(0, 1),
    layer.name = "sen2cloudness_probability"
  )

  mapshot(rmap4_target, url = sprintf("%s/map.html",normalizePath(dir_id)))
  unlink(sprintf("%s/map_files",normalizePath(dir_id)), recursive = TRUE)

  # RGB -REAL COLOR
  cloudsen2_plot(
    cloud_brick = final_stack[[c("B4","B3","B2")]],
    output = sprintf("%s/RGB.png",dir_id),
    expect_dim = c(511, 511)
  )

  # NRG - FALSE COLOR
  cloudsen2_plot(
    cloud_brick = final_stack[[c("B8","B4","B3")]],
    output = sprintf("%s/NRG.png",dir_id),
    expect_dim = c(511, 511)
  )

  # NDVI
  final_stack_ndvi <- (final_stack[["B8"]] - final_stack[["B4"]])/(final_stack[["B8"]] + final_stack[["B4"]])
  cloudsen2_plot(
    cloud_brick = final_stack_ndvi,
    expect_dim = c(511, 511),
    linStretchVec = FALSE,
    output = sprintf("%s/NDVI.png",dir_id),
    pal = viridis_pal(100),
    limits = c(-0.1, 0.4)
  )

  # NDSI
  final_stack_ndsi <- (final_stack[["B3"]] - final_stack[["B11"]])/(final_stack[["B3"]] + final_stack[["B11"]])
  cloudsen2_plot(
    cloud_brick = final_stack_ndsi,
    expect_dim = c(511, 511),
    output = sprintf("%s/NDSI.png",dir_id),
    linStretchVec = FALSE,
    pal = viridis_pal(100),
    limits = c(0.4, 1)
  )

  # CDI
  cdi <- final_stack[["cdi"]]
  cloudsen2_plot(
    cloud_brick = cdi,
    output = sprintf("%s/CDI.png", dir_id),
    expect_dim = c(511, 511),
    linStretchVec = FALSE,
    pal = viridis_pal(100),
    limits = c(-1, -0.4)
  )

  # CIRRUS
  cirrus <- final_stack[["B10"]]
  cloudsen2_plot(
    cloud_brick = cirrus,
    output = sprintf("%s/CIRRUS.png",dir_id),
    expect_dim = c(511, 511),
    linStretchVec = FALSE,
    pal = viridis_pal(100),
    limits = c(0.005, 0.01)
  )

  # POTENTIAL CLOUD
  pcloud <- final_stack[["probability"]]/100
  cloudsen2_plot(
    cloud_brick = pcloud,
    output = sprintf("%s/s2cloudness_prob.png",dir_id),
    expect_dim = c(511, 511),
    linStretchVec = FALSE,
    pal = viridis_pal(100),
    limits = c(0.4, 1)
  )

  # Target
  cloudsen2_plot(
    cloud_brick = final_stack[["sen2cor_mask"]],
    output = sprintf("%s/sen2cor_mask.png",dir_id),
    expect_dim = c(511, 511),
    linStretchVec = FALSE,
    pal = c("#007B00","#00FFFF", "#7D7D7D"),
    limits = c(0, 2)
  )
  cloudsen2_plot(
    cloud_brick = final_stack[["s2cloudness_mask"]],
    output = sprintf("%s/s2cloudness_mask.png",dir_id),
    expect_dim = c(511, 511),
    linStretchVec = FALSE,
    pal = c("#007B00","#00FFFF"),
    limits = c(0, 1)
  )

  tile_center <- cloudsen2_row$geometry[[1]]
  tibble(
    landcover = cloudsen2_row$type,
    source = s2_chip_id,
    tile_center_x = tile_center[1],
    tile_center_y = tile_center[2],
    avg_cloud_percentage = s2_ref_cloud,
    mean_solar_azimuth_angle = s2_ref_azimuth_angle,
    mean_solar_zenith_angle = s2_ref_zenith_angle
  ) -> chip_metadata

  write_json(chip_metadata, sprintf("%s/metadata.json", dir_id))
  new_dir_id <- sprintf("%s/%s_CC_%s_%s", output, index, round(s2_ref_cloud,2), s2_chip_id)
  file.rename(dir_id, new_dir_id)
  invisible(TRUE)
}

add_S2cloud <- function(img) {
  img %>%
    ee$Image$get("cloud_mask") %>%
    ee$Image() %>%
    ee$Image$select("probability") %>%
    img$addBands()
}

add_S2cloud2 <- function(img) {
  image_merge <- ee$Image(img %>% ee$Image$get("cloud_mask"))
  s2cloudness <- image_merge$select("probability")
  sen2cor <- image_merge$select("SCL")
  img$addBands(s2cloudness)$addBands(sen2cor)
}


# Determine the potential cloud probability desired
get_prob_by_class <- function(local_cloudsen2_points) {
  colnames <- sprintf("pcloud_%02d",1:5)
  pbaren_cloud <- local_cloudsen2_points %>% filter(value == 1) %>% nrow()
  baren_td <- gen_rcloudpoints(pbaren_cloud)  %>%
    `colnames<-`(colnames)

  pgrass_cloud <- local_cloudsen2_points %>% filter(value == 2) %>% nrow()
  grass_td <- gen_rcloudpoints(pgrass_cloud)  %>%
    `colnames<-`(colnames)

  pshrubland_cloud <- local_cloudsen2_points %>% filter(value == 3) %>% nrow()
  shrubland_td <- gen_rcloudpoints(pshrubland_cloud)  %>%
    `colnames<-`(colnames)

  psnow_cloud <- local_cloudsen2_points %>% filter(value == 4) %>% nrow()
  snow_td <- gen_rcloudpoints(psnow_cloud)  %>%
    `colnames<-`(colnames)

  ptmforest_cloud <- local_cloudsen2_points %>% filter(value == 5) %>% nrow()
  tmforest_td <- gen_rcloudpoints(ptmforest_cloud)  %>%
    `colnames<-`(colnames)

  ptrforest_cloud <- local_cloudsen2_points %>% filter(value == 6) %>% nrow()
  trforest_td <- gen_rcloudpoints(ptrforest_cloud)  %>%
    `colnames<-`(colnames)

  purban_cloud <- local_cloudsen2_points %>% filter(value == 7) %>% nrow()
  urban_td <- gen_rcloudpoints(purban_cloud)  %>%
    `colnames<-`(colnames)

  pwater_cloud <- local_cloudsen2_points %>% filter(value == 8) %>% nrow()
  water_td <- gen_rcloudpoints(pwater_cloud)  %>%
    `colnames<-`(colnames)

  pwetlands_cloud <- local_cloudsen2_points %>% filter(value == 9) %>% nrow()
  wetlands_td <- gen_rcloudpoints(pwetlands_cloud)  %>%
    `colnames<-`(colnames)

  complete_pprob <- rbind(baren_td,grass_td, shrubland_td, snow_td, tmforest_td,
                          trforest_td, urban_td, water_td, wetlands_td)

  cloudsen2_type <- local_cloudsen2_points$type
  cloudsen2_value <- local_cloudsen2_points$value
  fn_dataset <- cbind(id = seq_along(cloudsen2_value), type = cloudsen2_type, value = cloudsen2_value, complete_pprob)
  st_sf(fn_dataset, geometry = local_cloudsen2_points$geometry)
  # list(
  #   baren_prob = baren_td,
  #   grass_prob = grass_td,
  #   shrubland_prob = shrubland_td,
  #   snow_prob = snow_td,
  #   tmforest_prob = tmforest_td,
  #   trforest_prob = trforest_td,
  #   urban_prob = urban_td,
  #   water_prob = water_td,
  #   wetlands_prob = wetlands_td
  # )
}

# Cloud probability (range of interest) to pick up images
gen_rcloudpoints <- function(n) {
  # groups_n <- floor(c(0.05,0.3,0.3,0.3,0.05)*n)
  cloud_ppoints <- list()
  for (index in seq_len(n)) {
    groups_n <- rep(1, 5)
    # if (sum(groups_n) != n) {
    #   difff <- n - sum(groups_n)
    #   random_add <- sample(5,1)
    #   groups_n[random_add] <- groups_n[random_add] + difff
    # }
    cloud_ppoints[[index]] <- c(
      runif(n = groups_n[1], min = 0, max = 10),  # clear
      runif(n = groups_n[2], min = 10, max = 25), # almost clear
      runif(n = groups_n[3], min = 25, max = 45), # low-cloudy
      runif(n = groups_n[4], min = 45, max = 65), # mid-cloudy
      runif(n = groups_n[5], min = 65, max = 100)) # cloudy
  }
  cloud_ppoints %>%
    as.data.frame() %>%
    `colnames<-`(NULL) %>%
    t %>%
    as.data.frame()
}

cloud_fun_creator_s2cloudness <- function(scale,
                                          cloud_threshold = 40,
                                          convolution_kernel_radius = 4,
                                          dilatation_kernel_radius = 8) {
  function(img) {
    boxcar1 <- ee$Kernel$square(radius = convolution_kernel_radius,
                                units = 'pixels',
                                normalize = TRUE)
    boxcar2 <- ee$Kernel$square(radius = dilatation_kernel_radius,
                                units = 'pixels',
                                normalize = TRUE)
    cloud_prob <- img$select("probability")$
      convolve(boxcar1)$
      focal_max(kernel = boxcar2, iterations = 1)$
      gt(cloud_threshold)
    new_img <- img$select("probability")$updateMask(cloud_prob)$unmask(0)
    cprob <- new_img$reduceRegion(
      reducer = ee$Reducer$mean(),
      scale = scale
    )
    img$set(list(cprob = cprob$get("probability")))
  }
}

cloud_fun_creator_sen2cor <- function(scale, geom) {
  function(img) {
    cloud_mask <- sen2cor_reclass(img)$eq(1)
    cloud_mask_f <- cloud_mask$updateMask(cloud_mask)$unmask(0)
    # Map$addLayer(cloud_mask_f) + Map$addLayer(geom,name = "map_01")
    cprob <- cloud_mask_f$reduceRegion(
      reducer = ee$Reducer$mean(),
      geometry = geom,
      scale = scale
    )
    img$set(list(cprob = cprob$get("remapped")))
  }
}

#' Create a Sentinel2 Dataset with s2 bands, sen2cloudness
#' cloud probability, sen2cloudness cloud mask,   Scene
#' Classification map and sen2cor cloud mask
#' @param point ee$Geometry used to filter collections by bounds.
#' @param data_range range of dates used to filter collections by dates.
tsen2_dataset_creator <- function(point, data_range) {
  # 1. Create a S2 ImageCollection and filter by space and time.
  s2Sr <- ee$ImageCollection("COPERNICUS/S2_SR") %>%
    ee$ImageCollection$filterBounds(point) %>%
    ee$ImageCollection$filterDate(data_range[1], data_range[2])

  # 2. Create a S2 ImageCollection and filter by space and time.
  s2_level1c <- ee$ImageCollection("COPERNICUS/S2") %>%
    ee$ImageCollection$filterBounds(point) %>%
    ee$ImageCollection$filterDate(data_range[1], data_range[2])

  # 3. Create a CRS
  crs_kernel <- s2Sr$first()$select(0)$projection()$getInfo()$crs
  point_utm <- st_transform(cloudsen2_row$geometry[1], crs_kernel)
  ee_point <- ee$Geometry$Point(point_utm[[1]], proj = crs_kernel)

  # 4. Create a S2_CLOUD_PROBABILITY ImageCollection filtering by space and time.
  s2Clouds <- ee$ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY') %>%
    ee$ImageCollection$filterBounds(ee_point)%>%
    ee$ImageCollection$filterDate(data_range[1], data_range[2])

  # 5. Merge S2 and S2_CLOUD_PROBABILITY
  # 'add_S2cloud' is a function which convert the 'cloud' property to a band.
  s2SrWithCloudMask_1 <- ee$Join$saveFirst('cloud_mask')$apply(
    primary = s2Sr,
    secondary = s2Clouds,
    condition = ee$Filter$equals(leftField = 'system:index', rightField = 'system:index')
  ) %>% ee$ImageCollection() %>% ee$ImageCollection$map(add_S2cloud)

  s2SrWithCloudMask <- ee$Join$saveFirst('cloud_mask')$apply(
    primary = s2_level1c,
    secondary = s2SrWithCloudMask_1,
    condition = ee$Filter$equals(leftField = 'system:index', rightField = 'system:index')
  ) %>% ee$ImageCollection() %>% ee$ImageCollection$map(add_S2cloud2)

  list(ic = s2SrWithCloudMask,
       metadata = list(
         crs = crs_kernel,
         local_point = point_utm,
         ee_point = ee_point
       )
  )
}

sen2cor_reclass <- function(img) {
  s2_scl <- img$select("SCL")
  # 4,5,6,11 -> clear
  # 8,9,10 -> cloud
  # 2, 3 -> cloud shadows
  # 1, 7 -> no data
  # OBS: In Earth Engine "no data" values are masked out.
  s2_scl$remap(
    c(4, 5, 6, 11, 8, 9, 10, 2, 3, 1, 7),
    c(0, 0, 0, 0, 1, 1, 1, 2, 2, 1, 1)
  )
}

sen2cor_visparams <- function(){
  list(min = 0, max = 3, palette = c("#10d22c", "#00ffff", "#868686", "#ff0004"))
}



# Create cloudsen2 viz
cloudsen2_plot <- function(cloud_brick, output,
                           expect_dim = c(1000, 1000),
                           linStretchVec = TRUE,
                           pal = inferno_pal(100),
                           limits = NULL) {
  RGB <- as.matrix(getValues(cloud_brick))
  if (linStretchVec) {
    RGB_lin <- apply(RGB, 2, linStretchVec)
  } else {
    RGB_lin <- RGB*255
  }
  create_mtx <- function(x){
    matrix(RGB_lin[,x], nrow=nrow(cloud_brick), ncol=ncol(cloud_brick), byrow=T)
  }
  new_RGB <- simplify2array(lapply(seq_len(ncol(RGB)), create_mtx))
  png(output,expect_dim[1],expect_dim[2])
  if (ncol(RGB) == 1) {
    new_RGB[,,1] <- map2color(new_RGB[,,1]/255, pal,limits = limits)
    grid.raster(new_RGB[,,1])
  } else {
    grid.raster(new_RGB/255)
  }
  dev.off()
}

# Linear Stretch obtained from raster::plotRGB
linStretchVec <- function (x, quantile = c(0.02, 0.98)) {
  v <- stats::quantile(x, quantile, na.rm = TRUE)
  temp <- (255 * (x - v[1]))/(v[2] - v[1])
  temp[temp < 0] <- 0
  temp[temp > 255] <- 255
  return(temp)
}

#' Map over a serie of values
#' Function provided by Dave X
#' https://stackoverflow.com/questions/15006211/how-do-i-generate-a-mapping-from-numbers-to-colors-in-r
map2color <- function(x,pal,limits=NULL) {
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

# my palettes :)
inferno_pal <- function (n) {grDevices::hcl.colors(n, palette = "Inferno")}
viridis_pal <- function (n) {grDevices::hcl.colors(n, palette = "viridis")}


select_dataset_thumbnail_creator <- function(cloudsen2_row,
                                             n_images = 50,
                                             kernel_size = c(255, 255),
                                             data_range = c("2019-01-01", "2020-07-31"),
                                             output = "results/") {
  # 1. Create a point which represent the center of the chip
  point <- ee$Geometry$Point(cloudsen2_row$geometry[[1]])

  # 2. Create a S2 ImageCollection with all the necessary images.
  s2Sr <- ee$ImageCollection("COPERNICUS/S2_SR") %>%
    ee$ImageCollection$filterBounds(point) %>%
    ee$ImageCollection$filterDate(data_range[1], data_range[2])

  # 3. Create a S2 ImageCollection with all the necessary images.
  img_crs <- s2Sr$first()$select(0)$projection()$getInfo()[["crs"]]

  # 4. Donwload the images
  images_available <- s2Sr$size()$getInfo()
  if (images_available < n_images) {
    n_images <- images_available
  }
  images_position <- sample(images_available, n_images)

  # 4. Create folders to save results
  dir_id <- sprintf("%s/point_%04d",output, cloudsen2_row$id)
  dir.create(dir_id, showWarnings = FALSE)
  n_images <- seq_along(images_position)

  for (r_index in seq_along(images_position)) {
    # 4.1 Select the image to download
    img_to_download <- ee_get(s2Sr, images_position[r_index] - 1)$first()
    img_id <- img_to_download$id()$getInfo()
    s2_img_array <- img_to_download %>%
      ee$Image$select(c("B4", "B3", "B2")) %>%
      ee$Image$addBands(ee$Image$pixelCoordinates(projection = img_crs)) %>%
      # ee$Image$reproject(crs = img_crs, scale = 10) %>%
      ee$Image$neighborhoodToArray(
        kernel = ee$Kernel$rectangle(kernel_size[1], kernel_size[2], "pixels")
      ) %>%
      ee$Image$sampleRegions(ee$FeatureCollection(point),
                             projection = img_crs,
                             scale = 10) %>%
      ee$FeatureCollection$getInfo()

    # Convert data from list to data_frame
    message(
      sprintf("Processing point [%s] image [%s] ... please wait",
              cloudsen2_row$id, r_index)
    )
    band_names <- names(s2_img_array$features[[1]]$properties)
    extract_fn <- function(x) as.numeric(unlist(s2_img_array$features[[1]]$properties[x]))
    image_as_df <- do.call(cbind,lapply(band_names, extract_fn))
    colnames(image_as_df) <- band_names
    image_as_tibble <- as_tibble(image_as_df)
    coordinates(image_as_tibble) <- ~x+y
    sf_to_stack <- function(x) rasterFromXYZ(image_as_tibble[x])
    final_stack <- stack(lapply(names(image_as_tibble), sf_to_stack))
    crs(final_stack) <- st_crs(img_crs)$proj4string
    png(sprintf("%s/%s.png", dir_id, img_id), 1000, 1000)
    max_value <- max(maxValue(final_stack))
    plotRGB(final_stack/max_value, r = 3, g = 2, b = 1, scale = 1)
    dev.off()
  }

  # Apologize my messy code! :3
  xy <- cloudsen2_row[["geometry"]][[1]] %>% as.numeric()
  x <- xy[1]
  y <- xy[2]
  st_geometry(cloudsen2_row) <- NULL
  list_w_data <- cloudsen2_row %>%
    dplyr::select(starts_with("pcloud")) %>%
    as.list()
  names(list_w_data) <- sprintf("PUT_HERE_ID_%02d", 1:5)
  list_w_data$x <- x
  list_w_data$y <- y
  list_w_data$comments <- "PUT_HERE_YOUR_COMMENT"
  jsonlite::write_json(
      x = list_w_data,
      path = sprintf("%s/cprob_%04d.json", dir_id, cloudsen2_row$id),
      pretty = TRUE,
      auto_unbox = TRUE
  )
}
