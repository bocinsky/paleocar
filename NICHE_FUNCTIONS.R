setwd("/Users/bocinsky/IMPORTANT/WSU/RESEARCH/MyPublications/BOCINSKY_KOHLER_2014_NCC/R")
source('./UTILITY_FUNCTIONS.R')

pkgTest("raster")
pkgTest("rgdal")
pkgTest("maps")
pkgTest("RColorBrewer")

# Read in the data
VEPIIN.h2o_year.prcp.retro.unique <- brick("../OUTPUT/VEPIIN_PRISM_h20_year_prcp_unique_recons.tif")
VEPIIN.h2o_year.prcp.retro.union <- brick("../OUTPUT/VEPIIN_PRISM_h20_year_prcp_union_recons.tif")

VEPIIN.grow.GDD.retro.unique <- brick("../OUTPUT/VEPIIN_PRISM_grow_GDD_unique_recons.tif")
VEPIIN.grow.GDD.retro.union <- brick("../OUTPUT/VEPIIN_PRISM_grow_GDD_union_recons.tif")

VEPIIS.h2o_year.prcp.retro.unique <- brick("../OUTPUT/VEPIIS_PRISM_h20_year_prcp_unique_recons.tif")
VEPIIS.h2o_year.prcp.retro.union <- brick("../OUTPUT/VEPIIS_PRISM_h20_year_prcp_union_recons.tif")

VEPIIS.grow.GDD.retro.unique <- brick("../OUTPUT/VEPIIS_PRISM_grow_GDD_unique_recons.tif")
VEPIIS.grow.GDD.retro.union <- brick("../OUTPUT/VEPIIS_PRISM_grow_GDD_union_recons.tif")


# Calculate the PRCP and GDD NICHES

# PRCP niche is everything above 30 cm
calc.h2o_year.prcp.niche <- function(x) {
  x[x<300] <- 0
  x[x>=300] <- 1
  
  return(x)
}

VEPIIN.h2o_year.prcp.retro.unique.niche <- calc(VEPIIN.h2o_year.prcp.retro.unique, fun=calc.h2o_year.prcp.niche)
VEPIIN.h2o_year.prcp.retro.union.niche <- calc(VEPIIN.h2o_year.prcp.retro.union, fun=calc.h2o_year.prcp.niche)
VEPIIS.h2o_year.prcp.retro.unique.niche <- calc(VEPIIS.h2o_year.prcp.retro.unique, fun=calc.h2o_year.prcp.niche)
VEPIIS.h2o_year.prcp.retro.union.niche <- calc(VEPIIS.h2o_year.prcp.retro.union, fun=calc.h2o_year.prcp.niche)


# GDD niche is everything > 1800 GDD
calc.grow.GDD.niche <- function(x) {
  x[x<1800] <- 0
  x[x>=1800] <- 1
  
  return(x)
}

VEPIIN.grow.GDD.retro.unique.retro.niche <- calc(VEPIIN.grow.GDD.retro.unique, fun=calc.grow.GDD.niche)
VEPIIN.grow.GDD.retro.union.retro.niche <- calc(VEPIIN.grow.GDD.retro.union, fun=calc.grow.GDD.niche)
VEPIIS.grow.GDD.retro.unique.retro.niche <- calc(VEPIIS.grow.GDD.retro.unique, fun=calc.grow.GDD.niche)
VEPIIS.grow.GDD.retro.union.retro.niche <- calc(VEPIIS.grow.GDD.retro.union, fun=calc.grow.GDD.niche)


## Calculate the total niche, or the product of the PRCP and GDD niches
VEPIIN.unique.niche <- VEPIIN.h2o_year.prcp.retro.unique.niche * VEPIIN.grow.GDD.retro.unique.retro.niche
VEPIIN.union.niche <- VEPIIN.h2o_year.prcp.retro.union.niche * VEPIIN.grow.GDD.retro.union.retro.niche
VEPIIS.unique.niche <- VEPIIS.h2o_year.prcp.retro.unique.niche * VEPIIS.grow.GDD.retro.unique.retro.niche
VEPIIS.union.niche <- VEPIIS.h2o_year.prcp.retro.union.niche * VEPIIS.grow.GDD.retro.union.retro.niche


## Write all the niches to disk
writeRaster(VEPIIN.h2o_year.prcp.retro.unique.niche,"../OUTPUT/VEPIIN_h20_year_prcp_unique_niche.tif", datatype="LOG1S", options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND", "NBITS=1"), overwrite=T, setStatistics=FALSE)
writeRaster(VEPIIS.h2o_year.prcp.retro.unique.niche,"../OUTPUT/VEPIIS_h20_year_prcp_unique_niche.tif", datatype="LOG1S", options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND", "NBITS=1"), overwrite=T, setStatistics=FALSE)
writeRaster(VEPIIN.grow.GDD.retro.unique.retro.niche,"../OUTPUT/VEPIIN_grow_GDD_unique_niche.tif", datatype="LOG1S", options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND", "NBITS=1"), overwrite=T, setStatistics=FALSE)
writeRaster(VEPIIS.grow.GDD.retro.unique.retro.niche,"../OUTPUT/VEPIIS_grow_GDD_unique_niche.tif", datatype="LOG1S", options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND", "NBITS=1"), overwrite=T, setStatistics=FALSE)

writeRaster(VEPIIN.h2o_year.prcp.retro.union.niche,"../OUTPUT/VEPIIN_h20_year_prcp_union_niche.tif", datatype="LOG1S", options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND", "NBITS=1"), overwrite=T, setStatistics=FALSE)
writeRaster(VEPIIS.h2o_year.prcp.retro.union.niche,"../OUTPUT/VEPIIS_h20_year_prcp_union_niche.tif", datatype="LOG1S", options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND", "NBITS=1"), overwrite=T, setStatistics=FALSE)
writeRaster(VEPIIN.grow.GDD.retro.union.retro.niche,"../OUTPUT/VEPIIN_grow_GDD_union_niche.tif", datatype="LOG1S", options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND", "NBITS=1"), overwrite=T, setStatistics=FALSE)
writeRaster(VEPIIS.grow.GDD.retro.union.retro.niche,"../OUTPUT/VEPIIS_grow_GDD_union_niche.tif", datatype="LOG1S", options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND", "NBITS=1"), overwrite=T, setStatistics=FALSE)

writeRaster(VEPIIN.unique.niche,"../OUTPUT/VEPIIN_unique_niche.tif", datatype="LOG1S", options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND", "NBITS=1"), overwrite=T, setStatistics=FALSE)
writeRaster(VEPIIS.unique.niche,"../OUTPUT/VEPIIS_unique_niche.tif", datatype="LOG1S", options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND", "NBITS=1"), overwrite=T, setStatistics=FALSE)
writeRaster(VEPIIN.union.niche,"../OUTPUT/VEPIIN_union_niche.tif", datatype="LOG1S", options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND", "NBITS=1"), overwrite=T, setStatistics=FALSE)
writeRaster(VEPIIS.union.niche,"../OUTPUT/VEPIIS_union_niche.tif", datatype="LOG1S", options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND", "NBITS=1"), overwrite=T, setStatistics=FALSE)




