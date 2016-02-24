# PaleoCAR Tester

## Install and load PaleoCAR package
# install.packages("devtools")
# devtools::install_github("bocinsky/PaleoCAR")
library(PaleoCAR)
library(raster)
library(sp)

# Suppress use of scientific notation
options(scipen=999)

# Force Raster to load large rasters into memory
rasterOptions(chunksize=2e+07,maxmemory=2e+08)

# Set a directory for testing
testDir <- "~/PaleoCAR Test"

dir.create(testDir, showWarnings=F, recursive=T)
setwd(testDir)

######## BEGIN PARAMETERS

# @BEGIN main
# @PARAM calibration.years
# @PARAM prediction.years
# @IN mvnp
# @IN itrdb
# @IN mvnp_prism
# @OUT result_NEE_pdf  @URI file:result_NEE.pdf

## Set the calibration period
# Here, we use a 60 year period ending at 1983 
# to maximize the number of dendro series.
calibration.years <- 1924:1983

## Set the retrodiction years (AD)
prediction.years <- 500:1400

######## END PARAMETERS

## Load spatial polygon for the boundary of Mesa Verde National Park in southwestern Colorado:
data("mvnp")

## Get Tree-ring data from the ITRDB for 10-degree buffer around MVNP
data("itrdb")

## Get 1/3 arc-second PRISM data for the VEPII north study area
data("mvnp_prism")

## Run PaleoCAR over the MVNP brick
# @BEGIN paleoCAR.batch
# @PARAM calibration.years
# @PARAM prediction.years
# @PARAM out.dir
# @PARAM meanVar
# @PARAM floor
# @PARAM ceiling
# @PARAM asInt
# @PARAM force.redo
# @PARAM verbose
# @IN predictands  @AS mvnp_prism
# @IN label  @AS "mvnp_prism"
# @IN chronologies  @AS itrdb
# @OUT list  @AS mvnp_recon
mvnp_recon <- paleoCAR.batch(predictands = mvnp_prism, label = "mvnp_prism", chronologies = itrdb, calibration.years = calibration.years, prediction.years=prediction.years, out.dir=".", meanVar="chained", floor=0, ceiling=NULL, asInt=T, force.redo=T, verbose=T)
# @END paleoCAR.batch
# @END main

## Some plots...

# Mean landscape through time
plot(mean(mvnp_recon$recon))

# Mean across space
plot(y=cellStats(mvnp_recon$recon,mean),x=as.numeric(gsub("X","",names(mvnp_recon$recon))),ylab = "Precipitation (mm)", xlab = "Year AD", type='l')







