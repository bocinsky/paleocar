library(FedData)

##### MVNP Spatial Polygon

## Load spatial polygon for the boundary of Mesa Verde National Park in southwestern Colorado:
dir.create("./data-raw/MVNP", showWarnings = F, recursive = T)
# download shapefile directory
FedData::download_data("http://nrdata.nps.gov/programs/Lands/meve_tracts.zip", destdir="./data-raw/MVNP")
# uncompress shapefile directory
utils::unzip("./data-raw/MVNP/meve_tracts.zip", exdir="./data-raw/MVNP/meve_tracts")
# read shapefile
mvnp <- rgdal::readOGR("./data-raw/MVNP/meve_tracts", layer='MEVE_boundary')

unlink("./data-raw/MVNP", recursive=T)

rgdal::writeOGR(mvnp, dsn = "./data-raw/", layer = "mvnp", driver = "ESRI Shapefile", overwrite_layer = TRUE)

devtools::use_data(mvnp, overwrite = TRUE)

##### ITRDB DATA

## Set the calibration period
# Here, we use a 60 year period ending at 1983 
# to maximize the number of dendro series.
calibration.years <- 1924:1983

## Set the retrodiction years (AD)
prediction.years <- 1:2000

## Set a spatial buffer around the area you wish to reconstruct
## from which to grab tree-ring chronologies. This will probably
## be in degrees.
tree.buffer <- 10

## Get Tree-ring data from the ITRDB for 10-degree buffer around MVNP
# create a 10-degree buffer around the four corner states
treePoly <- suppressWarnings(rgeos::gBuffer(mvnp, width=tree.buffer, quadsegs=1000))
# extract the re-whitened residuals of the tree-ring chronologies
itrdb <- FedData::get_itrdb(template=treePoly, label="MVNP_PLUS_10DEG", raw.dir = "./data-raw/ITRDB/RAW/ITRDB/", extraction.dir = "./data-raw/ITRDB/EXTRACTIONS/ITRDB/", recon.years=prediction.years, calib.years=calibration.years, measurement.type="Ring Width", chronology.type="ARSTND", force.redo = TRUE)

unlink("./data-raw/ITRDB", recursive = T)

Encoding(levels(itrdb$metadata$NAME)) <- "latin1"
levels(itrdb$metadata$NAME) <- iconv(
  levels(itrdb$metadata$NAME), 
  "latin1", 
  "UTF-8"
)

Encoding(levels(itrdb$metadata$CONTRIBUTOR)) <- "latin1"
levels(itrdb$metadata$CONTRIBUTOR) <- iconv(
  levels(itrdb$metadata$CONTRIBUTOR), 
  "latin1", 
  "UTF-8"
)

devtools::use_data(itrdb, overwrite = TRUE)


##### PRISM WATER-YEAR PRECIPITATION
FedData::pkg_test("parallel")
# Load all the auxillary functions
source("./data-raw/annualize_prism_monthly.R")
## Read the monthly PRISM net precipitation values
mvnp_prism_ppt_monthly <- raster::brick("./data-raw/mvnp_prism_ppt_monthly.nc4")
raster::projection(mvnp_prism_ppt_monthly) <- raster::crs("EPSG:4326")
mvnp_prism <- annualize_prism_monthly(prism.brick = mvnp_prism_ppt_monthly, months = c(-2:9), fun = 'sum')
mvnp_prism <- mvnp_prism[[paste0("X",calibration.years)]]
devtools::use_data(mvnp_prism, overwrite = TRUE)
