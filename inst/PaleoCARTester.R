## Trying this for the Zuni/Cibola area. Will try other areas later.
# Set the working directory to the location of R-scripts
setwd("~/Dropbox/SKOPE/PaleoCARTest")
# setwd("/Users/Bocinsky/IMPORTANT/LTVTP/R")

# Load the functions for all analyses below
# install.packages("devtools")
# devtools::install_github("cran/wild1")
devtools::install_github("bocinsky/FedData")
devtools::install_github("bocinsky/PaleoCAR")
library(PaleoCAR)

# Suppress use of scientific notation
options(scipen=999)

# Force Raster to load large rasters into memory
rasterOptions(chunksize=2e+07,maxmemory=2e+08)

######## END PREAMBLE


## [INTENTIONALLY LEFT BLANK]


######## BEGIN PARAMETERS
## Set the master data directory
## This should be the location of PRISM data
# MASTER.DATA <- "/Users/bocinsky/Desktop/DATA/"
# MASTER.DATA <- "/Users/bocinsky/Desktop/LTVTP/DATA/"
# MASTER.DATA <- "/Volumes/BOCINSKY_DATA/DATA/"


## Set the calibration period
# Here, I use a 60 year period ending at 1983 
# to maximize the number of dendro series.
calibration.years <- 1924:1983
# calibration.years <- 1910:1919

## Set the retrodiction years
# Here, we are setting it for 1--2000,
# for consistency with Bocinsky & Kohler 2014
prediction.years <- 1:2000




######## END PARAMETERS


## [INTENTIONALLY LEFT BLANK]


#### TEMPORARY: Select smaller sub-area for testing
ZuniCibola.annual.prcp.full <- brick("./DATA/ZuniCibola_PRISM_annual_prcp.tif")
ZuniCibola.annual.prcp.10s <- crop(ZuniCibola.annual.prcp.full,extent(c(-109,-108.91,34,34.08)))
ZuniCibola.annual.prcp.32s <- crop(ZuniCibola.annual.prcp.full,extent(c(-109,-108.73,34,34.26)))
ZuniCibola.annual.prcp.100s <- crop(ZuniCibola.annual.prcp.full,extent(c(-109,-108.16,34,34.83)))
names(ZuniCibola.annual.prcp.full) <- calibration.years
names(ZuniCibola.annual.prcp.10s) <- calibration.years
names(ZuniCibola.annual.prcp.32s) <- calibration.years
names(ZuniCibola.annual.prcp.100s) <- calibration.years

ITRDB.data <- read.csv("./DATA/ITRDB_DATA.csv")

## [INTENTIONALLY LEFT BLANK]
ZuniCibola.full <- paleoCAR.batch(predictands=ZuniCibola.annual.prcp.full, label="ZuniCibola.annual.prcp.full", chronologies=ITRDB.data, calibration.years=calibration.years, prediction.years=prediction.years, out.dir="./OUTPUT/", meanVarMatch=T, floor=NULL, ceiling=NULL, asInt=F, force.redo=F, verbose=T)
ZuniCibola.10s <- paleoCAR.batch(predictands=ZuniCibola.annual.prcp.10s, label="ZuniCibola.annual.prcp.10s", chronologies=ITRDB.data, calibration.years=calibration.years, prediction.years=prediction.years, out.dir="./OUTPUT/", meanVarMatch=T, floor=NULL, ceiling=NULL, asInt=F, force.redo=F, verbose=T)
ZuniCibola.32s <- paleoCAR.batch(predictands=ZuniCibola.annual.prcp.32s, label="ZuniCibola.annual.prcp.32s", chronologies=ITRDB.data, calibration.years=calibration.years, prediction.years=prediction.years, out.dir="./OUTPUT/", meanVarMatch=T, floor=NULL, ceiling=NULL, asInt=F, force.redo=F, verbose=T)
ZuniCibola.100s <- paleoCAR.batch(predictands=ZuniCibola.annual.prcp.100s, label="ZuniCibola.annual.prcp.100s", chronologies=ITRDB.data, calibration.years=calibration.years, prediction.years=prediction.years, out.dir="./OUTPUT/", meanVarMatch=T, floor=NULL, ceiling=NULL, asInt=F, force.redo=F, verbose=T)



plot(ZuniCibola.10s[['recon']][[1923]])
plot(ZuniCibola.10s[['errors']][[1923]])
ZuniCibola.10s.errorMeans <- (ZuniCibola.10s[['errors']])/(ZuniCibola.10s[['recon']])
ZuniCibola.10s.errorMeans[is.na(ZuniCibola.10s.errorMeans)] <- 0
ZuniCibola.10s.errorMeans <- calc(ZuniCibola.10s.errorMeans,function(x){x[is.na(x)] <- 0; return(x)})

plot(calc(ZuniCibola.10s[['models']][['predictands']],mean))
plot((ZuniCibola.10s[['errors']][[1923]])/calc(ZuniCibola.10s[['models']][['predictands']],mean))
plot(abs((ZuniCibola.10s[['errors']][[1923]])/(ZuniCibola.10s[['recon']][[1923]])), zlim=c(0,0.1))

ZuniCibola.10s.sd <- raster::calc(ZuniCibola.10s[['recon']],sd)
ZuniCibola.10s.mean <- raster::calc(ZuniCibola.10s[['recon']],mean)
ZuniCibola.10s.z <- raster::calc(ZuniCibola.10s[['recon']],scale)

ZuniCibola.10s.z <- calc(ZuniCibola.10s.z,function(x){x[x>0] <- 1; return(x)})
ZuniCibola.10s.z <- calc(ZuniCibola.10s.z,function(x){x[x<0] <- -1; return(x)})


