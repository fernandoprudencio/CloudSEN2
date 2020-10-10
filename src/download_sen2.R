#' @description Download Sentinel 2 dataset
#' @author Fernando Prudencio

#' for linux distribution
#'   for more information, you have to see:
#'     http://sen2r.ranghetti.info/articles/installation.html#on-linux-systems
#'   put on terminal:
#'     sudo apt install -y r-base gdal-bin aria2 libpython2-dev libudunits2-dev
#'       libgdal-dev libjq-dev libprotobuf-dev protobuf-compiler libv8-dev
#'       libssl-dev libcairo2-dev

#' DEPENDENCIES
# GDAL, PROJ and GEOS (required by rgdal and rgeos)
# jq (required by jqr)
# protobuf (required by protolite)
# V8 (required by geojsonio)
# openssl (required by openssl)
# cairo (required by gdtools)
# libcurl (required by curl)
# NetCDF (required by ncdf4)
# libxml2 (required by XML)

rm(list = ls())

#' INSTALL PACKAGES
pkg <- c("sf", "raster", "sen2r")

sapply(
  pkg,
  function(x) {
    is.there <- x %in% rownames(installed.packages())
    if (is.there == FALSE) {
      install.packages(x, dependencies = T)
    }
  }
)

#' LOAD PACKAGE
library(sf)
library(raster)
library(sen2r)
library(mapedit)
library(tidyverse)

#' LOAD VECTOR DATA
# sf.area <- st_read(
#   dsn = "data/vector/gpkg/study_area.gpkg",
#   layer = "study_area_GEOWGS84", quiet = T, as_tibble = T
# )
sf.area <- mapedit::drawFeatures()

#' TIME PERIOD
ts <- as.Date(c("2020-01-01", "2020-09-15"))

#' LOGIN
#'   you must use your own ESA credentials
write_scihub_login(username = "geographyteam", password = "Geografia2019")

#' LIST OF SENTINEL DATASET
data <- s2_list(sf.area, time_interval = ts, max_cloud = 25, level = "L1C")

#' OBTAIN DATETIME OF SENTINEL DATASET
#'   however, you can also run View(data)
safe_getMetadata(data, "sensing_datetime")

#' DOWNLOAD SENTINEL DATASET
date()
s2_download(data[19], outdir = "data/raster/")
date()

#' TIME OF DOWNLOADING
#'   start time: Thu Sep 24 13:51:53 202
#'   end time: Thu Sep 24 14:11:33 2020
#'   time lapse: 20min to 17.1 Mbps
