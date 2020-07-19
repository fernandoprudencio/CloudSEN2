#' Create a Cloud Detection dataset
#'
#' This script helps to select the sentinel-2 tiles which will be used
#' to create the train/test/val dataset of the CLOUDSEN2.
#'
#' @author Cesar Aybar <csaybar.github.io>
#'

library(tidyverse)
library(mapview)
library(mapedit)
library(raster)
library(scales)
library(stars)
library(rgee)
library(png)
library(sf)
library(sp)

set.seed(101)
source("src/utils.R")

ee_Initialize("csaybar", gcs = TRUE)

# Load potential points
local_cloudsen2_points <- read_sf("data/cloudsen2.geojson") %>% arrange(type)
create_potential_prob <- get_prob_by_class(local_cloudsen2_points) # potential points
local_cloudsen2_points$potential_probability <-  create_potential_prob

cloudsen2_row <- local_cloudsen2_points[20,]
dataset_creator(cloudsen2_row)

