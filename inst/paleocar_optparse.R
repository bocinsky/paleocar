#!/usr/bin/env Rscript
# USAGE: Rscript paleocar_optparse.R [options]
# MANUAL: Rscript paleocar_optparse.R --help

# Install the devtools package if not installed already
# Devtools' install functions do SHA1 matching
if(!suppressWarnings(require("devtools", quietly = TRUE)))
  suppressWarnings({
    utils::install.packages("devtools",
                            dependencies = TRUE,
                            repos = "http://cran.rstudio.com",
                            quiet = TRUE)
  })

# Python-style argument parsing
devtools::install_cran("optparse",
                       repos = "http://cran.rstudio.com",
                       quiet = TRUE)

# Use optparse for Python-style argument parsing
# Define a set of options; only output file and force_redo for now
option_list = list(
  optparse::make_option(c("--label"),
                        type = "character",
                        default = "test",
                        help = "A character label for the reconstruction, for saving [default = %default]"),
  
  optparse::make_option(c("-o", "--output_dir"),
                        type = "character",
                        default = "./out/",
                        help = "Output directory [default = %default]",
                        metavar = "character"),
  
  optparse::make_option(c("--predictands"),
                        type = "character",
                        default = "./mvnp_prism.rds",
                        help = "Path to an R data set (*.rds) containing the predictand raster (i.e., modern climate data) [default = %default]"),
  
  optparse::make_option(c("--chronologies"),
                        type = "character",
                        default = "./itrdb.rds",
                        help = "Path to an R data set (*.rds) containing the tree ring chronologies [default = %default]"),
  
  optparse::make_option(c("-c", "--clean"),
                        action = "store_true",
                        default = FALSE,
                        help = "Delete output directory and re-run all analyses? [default= %default]"),
  
  optparse::make_option(c("-v", "--verbose"),
                        action = "store_true",
                        default = TRUE,
                        help = "Print extra output [default= %default]"),
  
  optparse::make_option(c("-q", "--quietly"),
                        action = "store_false",
                        dest = "verbose",
                        help = "Print little output [default= %default]")
)

opt_parser = optparse::OptionParser(option_list = option_list);
opt = optparse::parse_args(opt_parser);

## Install and load paleocar and FedData packages
devtools::install_github(c("bocinsky/FedData","bocinsky/paleocar"),
                         quiet = opt$verbose)
devtools::install_cran(c("readr","stringr"),
                       repos = "http://cran.rstudio.com",
                       quiet = opt$verbose)
library(FedData)
library(paleocar)
library(raster)

# Force Raster to load large rasters into memory
raster::rasterOptions(chunksize=2e+07,maxmemory=2e+08)

# Delete the output directory if requested, then create it
opt$output_dir <- ifelse(stringr::str_sub(opt$output_dir,-1) != "/",
                         stringr::str_c(opt$output_dir,"/"),
                         opt$output_dir) 
if(opt$clean) unlink(opt$output_dir,
                     recursive = TRUE,
                     force = TRUE)
dir.create(opt$output_dir,
           showWarnings = FALSE,
           recursive = TRUE)
# a function that builds output paths
out <- function(...){
  stringr::str_c(opt$output_dir,...)
}

######## BEGIN PARAMETERS
##### DECLARE OTHER PARAMETERS HERE! #####

## Set the calibration period
# Here, we use a 60 year period ending at 1983 
# to maximize the number of dendro series.
calibration.years <- 1924:1983

######## END PARAMETERS

# Load the chronologies
tryCatch(chronologies <- readr::read_rds(opt$chronologies), 
         error = function(e){
           stop("Chronologies file does not exist!")
         })

# Load the predictands
tryCatch(predictands <- readr::read_rds(opt$predictands), 
         error = function(e){
           stop("Predictands file does not exist!")
         })


## Run paleocar over the predictands brick
recon <- paleocar(chronologies = chronologies,
                  predictands = predictands,
                  calibration.years = calibration.years,
                  label = opt$label,
                  out.dir = opt$output_dir,
                  force.redo = opt$clean,
                  verbose = opt$verbose)

message("paleocar reconstruction complete. Output is in '", opt$output_dir,"'.")

# # Example of reading output prediction raster brick
# test.predictions <- readr::read_rds(paste0(opt$output_dir,"/test.prediction.Rds"))
# test.predictions
# plot(test.predictions)
# 
# # Example of reading output uncertainty raster brick
# test.uncertainty <- readr::read_rds(paste0(opt$output_dir,"/test.uncertainty.Rds"))
# test.uncertainty
# plot(test.uncertainty)
