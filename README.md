paleocar
========

[![Build Status](https://api.travis-ci.org/bocinsky/paleocar.png)](https://travis-ci.org/bocinsky/paleocar)

`paleocar` is an *R* package implementing functions to perform spatio-temporal paleoclimate reconstruction from tree-rings using the CAR (Correlation Adjusted corRelation) approach of Zuber and Strimmer (2011). It is optimized for speed and memory use.

This is based on the approach used in Bocinsky and Kohler (2014): 

Bocinsky, R. K. and Kohler, T. A. (2014). A 2,000-year reconstruction of the rain-fed maize agricultural niche in the US Southwest. *Nature Communications*, 5:5618. doi: [10.1038/ncomms6618](http://www.nature.com/ncomms/2014/141204/ncomms6618/full/ncomms6618.html).

The primary difference is that here model selection is performed by minimizing the corrected Akaike's Information Criterion.

A more recent reference would be Bocinsky et al. (2016):

Bocinsky, R. K., Rush, J., Kintigh, K. W., and Kohler, T. A. (2016). Exploration and exploitation in the macrohistory of the pre-Hispanic Pueblo Southwest. *Science Advances*, 2:[e1501532](http://advances.sciencemag.org/content/2/4/e1501532).

This package has been built and tested on a source (Homebrew) install of *R* on macOS 10.12 (Sierra), and has been successfully run on Ubuntu 14.04.5 LTS (Trusty), Ubuntu 16.04.1 LTS (Xenial) and binary installs of *R* on Mac OS 10.12 and Windows 10.

### Development
+ [Kyle Bocinsky](http://bocinsky.io) - Crow Canyon Archaeological Center, Cortez, CO

### Install `paleocar`
+ Development version from GitHub:
```r
install.packages("devtools")
devtools::install_github("bocinsky/paleocar")
library(paleocar)
```
+ Linux (Ubuntu 14.04.5 or 16.04.1):

    First, in terminal:
    ```r
    sudo add-apt-repository ppa:ubuntugis/ppa -y
    sudo apt-get update -q
    sudo apt-get install libssl-dev libcurl4-openssl-dev netcdf-bin libnetcdf-dev gdal-bin libgdal-dev
    ```
    Then, in R:
    ```r
    update.packages("survival")
    install.packages("devtools")
    devtools::install_github("bocinsky/paleocar")
    library(paleocar)
    ```

### Demonstration
This demo script is available in the `/inst` folder at the location of the installed package.

#### Load `paleocar` and set a working directory
```r
library(paleocar)
library(magrittr) # The magrittr package enables piping in R.

# Set a directory for testing
testDir <- "~/paleocar Test"
# and create it if necessary
dir.create(testDir, showWarnings=F, recursive=T)
setwd(testDir)
  
```

#### Load test datasets
`paleocar` ships with test files defining a study area (Mesa Verde National Park), and pre-extracted data from the International Tree Ring Databank using the [`FedData` package](https://github.com/bocinsky/FedData). See the `data-raw/data.R` script (or the documentation for `FedData`) to learn how to download these data.
```r
# Load spatial polygon for the boundary of Mesa Verde National Park (MVNP) in southwestern Colorado:
data(mvnp)
  
# Get Tree-ring data from the ITRDB for 10-degree buffer around MVNP
data(itrdb)
  
# Get 1/3 arc-second PRISM gridded data for the MVNP north study area
data(mvnp_prism)
  
```

#### Run `paleocar`
```r
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
```

#### Plot results
```r
mvnp_recon$recon %>%
  mean() %>%
  raster::plot()
```
  