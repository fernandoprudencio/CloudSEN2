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


#' Create 511x511 SENTINEL-2 (level-1C) chips to feed a DL model
dataset_creator <- function(cloudsen2_row,
                            bands = "B.*|probability",
                            kernel_size = c(255, 255),
                            data_range = c("2015-01-01", "2016-12-31"),
                            output = "results/") {

  # 1. Create a point which represent the center of the chip
  point <- ee$Geometry$Point(cloudsen2_row$geometry[[1]])

  # Create a S2 ImageCollection filtering by space and time.
  s2Sr <- ee$ImageCollection("COPERNICUS/S2") %>%
    ee$ImageCollection$filterBounds(point) %>%
    ee$ImageCollection$filterDate(data_range[1], data_range[2])

  # Create a S2_CLOUD_PROBABILITY ImageCollection filtering by space and time.
  s2Clouds <- ee$ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY') %>%
    ee$ImageCollection$filterBounds(point)%>%
    ee$ImageCollection$filterDate(data_range[1], data_range[2])

  # Merge S2 and S2_CLOUD_PROBABILITY
  # 'add_S2cloud' is a function which convert the 'cloud' property
  # to a band.
  s2SrWithCloudMask <- ee$Join$saveFirst('cloud_mask')$apply(
    primary = s2Sr,
    secondary = s2Clouds,
    condition = ee$Filter$equals(leftField = 'system:index', rightField = 'system:index')
  ) %>% ee$ImageCollection() %>% ee$ImageCollection$map(add_S2cloud)

  # Create a function to estimate the mean cloud propability
  cloud_prob_fn <- cloud_fun_creator(point = point, kernel = 2555)

  # Apply 'cloud_fun_creator' to s2SrWithCloudMask
  s2SrWithCloudMask %>%
    ee$ImageCollection$map(cloud_prob_fn) %>%
    ee$ImageCollection$aggregate_array("cprob") %>%
    ee$List$getInfo() %>%
    unlist -> cloud_percentage

  lowest_position <- which.min(abs(cloud_percentage)) - 1
  ideal_position <- which.min(abs(cloud_percentage - cloudsen2_row$potential_probability)) - 1
  to_dowload <- c(lowest_position, ideal_position)

  # Create a folder
  s2_chip_id <- ee_get(s2SrWithCloudMask, ideal_position)$first()$
    get("system:id")$getInfo() %>% basename()
  dir_id <- sprintf("%s/%s",output, s2_chip_id)
  dir.create(dir_id, showWarnings = FALSE)
  folder_names <- c("reference", "target")

  for (index in seq_along(to_dowload)) {
    # Select an Image considering a position
    s2_img_data <- ee_get(s2SrWithCloudMask, to_dowload[index]) %>%
      ee$ImageCollection$first()

    # Retrieve data
    s2_img_array <- s2_img_data %>%
      ee$Image$select(bands) %>%
      ee$Image$addBands(ee$Algorithms$Sentinel2$CDI(s2_img_data)) %>%
      ee$Image$addBands(ee$Image$pixelLonLat()) %>%
      ee$Image$reproject(crs = "EPSG:4326", scale = 10) %>%
      ee$Image$neighborhoodToArray(
        kernel = ee$Kernel$rectangle(kernel_size[1], kernel_size[2], "pixels")
      ) %>%
      ee$Image$sampleRegions(ee$FeatureCollection(point)) %>%
      ee$FeatureCollection$getInfo()

    # Convert data from list to data_frame
    s2_chip_id <- basename(s2_img_data$get("system:id")$getInfo())
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
    final_stack <- final_stack/10000
    final_stack[["cdi"]] <- final_stack[["cdi"]]*10000
    final_stack[["probability"]] <- final_stack[["probability"]]*10000

    # Create a folder
    dir_id_f <- sprintf("%s/%s",dir_id, folder_names[index])
    dir.create(dir_id_f, showWarnings = FALSE)

    # Create and save a mapview object
    if (index == 2) {
      r_aux <- mapview(
        x = final_stack[[c("B10", "cdi", "probability")]],
        map = cleanRmap@map,
        layer.name = c("cirrus", "cdi", "cloud_prob")
      )
    } else {
      r_aux <- mapview(
        x = final_stack[[c("B10", "cdi", "probability")]],
        layer.name = c("cirrus", "cdi", "cloud_prob")
      )
    }

    Rmap <- viewRGB(
      x = final_stack[[c(5,6,7)]],
      map =  r_aux@map,
      quantiles = c(0, 1),
      layer.name = sprintf("RGB_%s", folder_names[index])
    )

    if (index == 1) {
      cleanRmap <- viewRGB(
        x = final_stack[[c(5,6,7)]],
        quantiles = c(0, 1),
        layer.name = sprintf("RGB_%s", folder_names[index])
      )
    }
    mapshot(Rmap, url = sprintf("%s/map.html",normalizePath(dir_id_f)))

    # RGB -REAL COLOR
    cloudsen2_plot(final_stack[[c("B4","B3","B2")]], sprintf("%s/RGB.png",dir_id_f))

    # NRG - FALSE COLOR
    cloudsen2_plot(final_stack[[c("B8","B4","B3")]], sprintf("%s/NRG.png",dir_id_f))

    # NDVI
    final_stack_ndvi <- (final_stack[["B8"]] - final_stack[["B4"]])/(final_stack[["B8"]] + final_stack[["B4"]])
    cloudsen2_plot(
      cloud_brick = final_stack_ndvi,
      linStretchVec = FALSE,
      output = sprintf("%s/NDVI.png",dir_id_f),
      pal = viridis_pal(100),
      limits = c(-0.1, 0.4)
    )

    # NDSI
    final_stack_ndsi <- (final_stack[["B3"]] - final_stack[["B11"]])/(final_stack[["B3"]] + final_stack[["B11"]])
    cloudsen2_plot(
      cloud_brick = final_stack_ndsi,
      output = sprintf("%s/NDSI.png",dir_id_f),
      linStretchVec = FALSE,
      pal = viridis_pal(100),
      limits = c(0.4, 1)
    )

    # CDI
    cdi <- final_stack[["cdi"]]
    cloudsen2_plot(
      cloud_brick = cdi,
      output = sprintf("%s/CDI.png",dir_id_f),
      linStretchVec = FALSE,
      pal = viridis_pal(100),
      limits = c(-1, -0.4)
    )

    # CIRRUS
    cirrus <- final_stack[["B10"]]
    cloudsen2_plot(
      cloud_brick = cirrus,
      output = sprintf("%s/CIRRUS.png",dir_id_f),
      linStretchVec = FALSE,
      pal = viridis_pal(100),
      limits = c(0.005, 0.01)
    )

    # POTENTIAL CLOUD
    pcloud <- final_stack[["probability"]]/100
    cloudsen2_plot(
      cloud_brick = pcloud,
      output = sprintf("%s/pcloud.png",dir_id_f),
      linStretchVec = FALSE,
      pal = viridis_pal(100),
      limits = c(0.4, 1)
    )

    # Target
    target <- pcloud > 0.7
    cloudsen2_plot(
      cloud_brick = target,
      output = sprintf("%s/01_s2cloudless_TARGET.png",dir_id_f),
      linStretchVec = FALSE,
      pal = c("#00FF00","#FFB000"),
      limits = c(0, 1)
    )
    cloudsen2_plot(
      cloud_brick = target,
      output = sprintf("%s/01_manual_TARGET.png",dir_id_f),
      linStretchVec = FALSE,
      pal = c("#00FF00","#FFB000"),
      limits = c(0, 1)
    )

    tile_center <- cloudsen2_row$geometry[[1]]
    tibble(
      landcover = cloudsen2_row$type,
      source = s2_chip_id,
      tile_center_x = tile_center[1],
      tile_center_y = tile_center[2],
      avg_cloud_percentage = mean(getValues(pcloud))
    ) -> chip_metadata
    write_csv(chip_metadata, sprintf("%s/metadata.csv",dir_id_f))
  }
}

add_S2cloud <- function(img) {
  img %>%
    ee$Image$get("cloud_mask") %>%
    ee$Image() %>%
    ee$Image$select("probability") %>%
    img$addBands()
}

# Determine the potential cloud probability desired
get_prob_by_class <- function(local_cloudsen2_points) {
  pbaren_cloud <- local_cloudsen2_points %>% filter(value == 1) %>% nrow()
  baren_td <- gen_rcloudpoints(pbaren_cloud)

  pgrass_cloud <- local_cloudsen2_points %>% filter(value == 2) %>% nrow()
  grass_td <- gen_rcloudpoints(pgrass_cloud)

  pshrubland_cloud <- local_cloudsen2_points %>% filter(value == 3) %>% nrow()
  shrubland_td <- gen_rcloudpoints(pshrubland_cloud)

  psnow_cloud <- local_cloudsen2_points %>% filter(value == 4) %>% nrow()
  snow_td <- gen_rcloudpoints(psnow_cloud)

  ptmforest_cloud <- local_cloudsen2_points %>% filter(value == 5) %>% nrow()
  tmforest_td <- gen_rcloudpoints(ptmforest_cloud)

  ptrforest_cloud <- local_cloudsen2_points %>% filter(value == 6) %>% nrow()
  trforest_td <- gen_rcloudpoints(ptrforest_cloud)

  purban_cloud <- local_cloudsen2_points %>% filter(value == 7) %>% nrow()
  urban_td <- gen_rcloudpoints(purban_cloud)

  pwater_cloud <- local_cloudsen2_points %>% filter(value == 8) %>% nrow()
  water_td <- gen_rcloudpoints(pwater_cloud)

  pwetlands_cloud <- local_cloudsen2_points %>% filter(value == 9) %>% nrow()
  wetlands_td <- gen_rcloudpoints(pwetlands_cloud)

  c(
    baren_prob = baren_td,
    grass_prob = grass_td,
    shrubland_prob = shrubland_td,
    snow_prob = snow_td,
    tmforest_prob = tmforest_td,
    trforest_prob = trforest_td,
    urban_prob = urban_td,
    water_prob = water_td,
    wetlands_prob = wetlands_td
  )

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
  groups_n <- floor(c(0.05,0.3,0.3,0.3,0.05)*n)
  if (sum(groups_n) != n) {
    difff <- n - sum(groups_n)
    random_add <- sample(5,1)
    groups_n[random_add] <- groups_n[random_add] + difff
  }
  c(runif(n = groups_n[1], min = 0, max = 10),
    runif(n = groups_n[2], min = 10, max = 20), # almost clear
    runif(n = groups_n[3], min = 20, max = 40), # low-cloudy
    runif(n = groups_n[4], min = 40, max = 65), # mid-cloudy
    runif(n = groups_n[5], min = 65, max = 100)) # cloudy
}

cloud_fun_creator <- function(img, point, kernel = 2555) {
  function(img) {
    ee_rect_kernel <- point$buffer(kernel)$bounds()
    cloud_prob <- img$select("probability")
    cprob <- cloud_prob$reduceRegion(
      reducer = ee$Reducer$mean(),
      geometry = ee_rect_kernel
    )
    img$set(list(cprob = cprob$get("probability")))
  }
}


# Create cloudsen2 viz
cloudsen2_plot <- function(cloud_brick, output,
                           kernel_size = c(255, 255),
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
  expect_dim <- kernel_size*2+1
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
linStretchVec <- function (x) {
  v <- stats::quantile(x, c(0, 1), na.rm = TRUE)
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

