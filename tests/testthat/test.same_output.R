library(raster)
testthat::context("paleocar tests")

testthat::test_that("Test that the paleocar_models and paleocar_models_simple functions give same results", {
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
  
  # Get a vector for a single cell
  mvnp_prism_single <- mvnp_prism[1][1,]
  
  ## Run paleocar_models over the MVNP vector
  paleocar_models.out <- paleocar_models(predictands = mvnp_prism_single,
                  chronologies = itrdb,
                  calibration.years = calibration.years,
                  prediction.years = prediction.years)
  
  ## Run paleocar_models_simple over the MVNP vector
  paleocar_models_simple.out <- paleocar_models_simple(predictands = mvnp_prism_single,
                  chronologies = itrdb,
                  calibration.years = calibration.years,
                  prediction.years = prediction.years)
  
  testthat::expect_equal(paleocar_models.out,
                         paleocar_models_simple.out)
  
})
