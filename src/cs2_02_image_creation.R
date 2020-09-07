#' Create a Cloud Detection dataset
#'
#' This script helps to select the sentinel-2 tiles which will be used
#' to create the train/test/val dataset of the CLOUDSEN2.
#'
#' @author Cesar Aybar <csaybar.github.io>
#'

library(tidyverse)
library(jsonlite)
library(mapview)
library(mapedit)
library(raster)
library(scales)
library(stars)
library(grid)
library(rgee)
library(png)
library(sf)
library(sp)

set.seed(101)
source("src/utils.R")

ee_Initialize("csaybar")

# Load potential points
local_cloudsen2_points <- read_sf("data/cloudsen2.geojson") %>% arrange(type)
create_potential_prob <- get_prob_by_class(local_cloudsen2_points) # potential points
local_cloudsen2_points$potential_probability <-  create_potential_prob

# index <- 30
for (index in 1:1000) {
  cloudsen2_row <- local_cloudsen2_points[index,]
  dataset_creator_chips(
    cloudsen2_row = cloudsen2_row,
    bands = "B.*|probability|SCL",
    kernel_size = c(255, 255),
    data_range = c("2019-01-01", "2020-07-31"),
    output = "results/"
  )
}
