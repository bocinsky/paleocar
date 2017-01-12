library(raster)
testthat::context("paleocar batch tests")

testthat::test_that("The Daymet tiles are available at the correct URL", {
  # Suppress use of scientific notation
  options(scipen=999)
  
  # Force Raster to load large rasters into memory
  raster::rasterOptions(chunksize=2e+07,maxmemory=2e+08)
  
  ## Set the calibration period
  # Here, we use a 60 year period ending at 1983 
  # to maximize the number of dendro series.
  calibration.years <- 1924:1983
  
  ## Set the retrodiction years (AD)
  prediction.years <- 500:1400
  
  ######## END PARAMETERS
  
  ## Load spatial polygon for the boundary of Mesa Verde National Park in southwestern Colorado:
  data(mvnp)
  
  ## Get Tree-ring data from the ITRDB for 10-degree buffer around MVNP
  data(itrdb)
  
  ## Get 1/3 arc-second PRISM data for the VEPII north study area
  data(mvnp_prism)
  
  ## Run PaleoCAR over the MVNP brick
  testthat::expect_warning(
    mvnp_recon <- paleocar_batch(predictands = mvnp_prism,
                               label = "mvnp_prism",
                               chronologies = itrdb,
                               calibration.years = calibration.years,
                               prediction.years=prediction.years,
                               out.dir="./",
                               meanVar="chained",
                               floor=0,
                               ceiling=NULL,
                               asInt=T,
                               force.redo=T,
                               verbose=T)
    )
  
  testthat::expect_is(mvnp_recon,"list")
  
  testthat::expect_named(mvnp_recon)
  
  testthat::expect_is(mvnp_recon$recon, "RasterBrick")
  
  unlink("mvnp_prism.models.rds")
  unlink("mvnp_prism.recon.tif")
  unlink("mvnp_prism.carscores.Rds")
})
