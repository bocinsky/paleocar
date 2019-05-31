
# paleocar

[![Build
Status](https://api.travis-ci.org/bocinsky/paleocar.png)](https://travis-ci.org/bocinsky/paleocar)
<!-- [![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/paleocar)](https://github.com/metacran/cranlogs.app) -->
[![cran
version](http://www.r-pkg.org/badges/version/paleocar)](https://cran.r-project.org/package=paleocar)

`paleocar` is an *R* package implementing functions to perform
spatio-temporal paleoclimate reconstruction from tree-rings using the
CAR (Correlation Adjusted corRelation) approach of Zuber and Strimmer as
implemented in the [`care`
package](https://CRAN.R-project.org/package=care) for *R*. It is
optimized for speed and memory use.

This is based on the approach used in Bocinsky and Kohler (2014):

Bocinsky, R. K. and Kohler, T. A. (2014). A 2,000-year reconstruction of
the rain-fed maize agricultural niche in the US Southwest. *Nature
Communications*, 5:5618. doi:
[10.1038/ncomms6618](http://www.nature.com/ncomms/2014/141204/ncomms6618/full/ncomms6618.html).

The primary difference between the latest version of `paleocar` and that
presented in Bocinsky and Kohler (2014) is, here, model selection is
performed by minimizing the corrected Akaike’s Information Criterion.

A more recent reference would be Bocinsky et al. (2016):

Bocinsky, R. K., Rush, J., Kintigh, K. W., and Kohler, T. A. (2016).
Exploration and exploitation in the macrohistory of the pre-Hispanic
Pueblo Southwest. *Science Advances*,
2:[e1501532](http://advances.sciencemag.org/content/2/4/e1501532).

This package has been built and tested on a source (Homebrew) install of
*R* on macOS 10.12 (Sierra), and has been successfully run on Ubuntu
14.04.5 LTS (Trusty), Ubuntu 16.04.1 LTS (Xenial) and binary installs of
*R* on Mac OS 10.12 and Windows 10.

### Development

  - [Kyle Bocinsky](http://bocinsky.io) - Crow Canyon Archaeological
    Center, Cortez, CO

### Install `paleocar`

  - Development version from GitHub:

<!-- end list -->

``` r
install.packages("devtools")
devtools::install_github("bocinsky/paleocar")
library(paleocar)
```

  - Linux (Ubuntu 14.04.5 or 16.04.1):

First, in terminal:

``` bash
sudo add-apt-repository ppa:ubuntugis/ppa -y
sudo apt-get update -q
sudo apt-get install libssl-dev libcurl4-openssl-dev netcdf-bin libnetcdf-dev gdal-bin libgdal-dev
```

Then, in R:

``` r
update.packages("survival")
install.packages("devtools")
devtools::install_github("bocinsky/paleocar")
library(paleocar)
```

### Demonstration

This demo script is available in the `/inst` folder at the location of
the installed package.

#### Load `paleocar` and set a working directory

``` r
library(paleocar)
```

    ## Registered S3 methods overwritten by 'ggplot2':
    ##   method         from 
    ##   [.quosures     rlang
    ##   c.quosures     rlang
    ##   print.quosures rlang

    ## Registered S3 method overwritten by 'xts':
    ##   method     from
    ##   as.zoo.xts zoo

    ## Registered S3 method overwritten by 'quantmod':
    ##   method            from
    ##   as.zoo.data.frame zoo

    ## Registered S3 methods overwritten by 'forecast':
    ##   method             from    
    ##   fitted.fracdiff    fracdiff
    ##   residuals.fracdiff fracdiff

``` r
library(magrittr) # The magrittr package enables piping in R.
library(ggplot2)

# Set a directory for testing
testDir <- "./paleocar_test/"
# and create it if necessary
dir.create(testDir, showWarnings=F, recursive=T)
```

#### Load test datasets

`paleocar` ships with test files defining a study area (Mesa Verde
National Park), and pre-extracted data from the International Tree Ring
Databank using the [`FedData`
package](https://github.com/bocinsky/FedData). See the `data-raw/data.R`
script (or the documentation for `FedData`) to learn how to download
these
data.

``` r
# Load spatial polygon for the boundary of Mesa Verde National Park (MVNP) in southwestern Colorado:
data(mvnp)

# Get Tree-ring data from the ITRDB for 10-degree buffer around MVNP
data(itrdb)

# Get 1/3 arc-second PRISM gridded data for the MVNP north study area (water-year [October--September] precipitation, in millimeters)
data(mvnp_prism)
```

#### Run `paleocar`

`paleocar` can be run for either single location given by a vector of
annualized climate data, a matrix of locations, or over gridded climate
data such as PRISM in raster format. There are three primary functions:

  - `paleocar_models()` calculates the CAR-ranked linear models for all
    reconstructions
  - `predict_paleocar_models()` generates climate predictions over a
    specified prediction period, and
  - `uncertainty_paleocar_models()` generates an estimate of model
    uncertainty over a specified prediction period.

Finally, the `paleocar()` method is a convenience wrapper that runs all
three of these functions and returns a list with their output. See the
documentation for each function for details.

##### `paleocar` reconstruction for a single location

`paleocar` may be run for a single location by providing a vector of
annualized values to be reconstructed. Simply provide a numeric vector
the same length as your calibration years as the `predictands`
parameter.

``` r
# Extract a vector of annualized climate data (the first cell in the raster)
mvnp_prism.vector <- mvnp_prism[1][1,]

test.vector <- paleocar_models(predictands = mvnp_prism.vector,
                               chronologies = itrdb,
                               calibration.years = 1924:1983,
                               prediction.years = 1:2000,
                               verbose = T)
```

    ## Calculating PaleoCAR models
    ## 
    ## Prepare data and calculate CAR scores: 0 minutes
    ## 
    ## Calculating models of with 1 input vectors.
    ## Define models: 0.02 minutes
    ## Calculate 5 linear models: 0 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.03 minutes
    ## 123 cell-years remaining
    ## 
    ## Calculating models of with 2 input vectors.
    ## Define models: 0.02 minutes
    ## Calculate 7 linear models: 0 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.03 minutes
    ## 115 cell-years remaining
    ## 
    ## Calculating models of with 3 input vectors.
    ## Define models: 0.02 minutes
    ## Calculate 7 linear models: 0 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.03 minutes
    ## 41 cell-years remaining
    ## 
    ## Calculating models of with 4 input vectors.
    ## Define models: 0.02 minutes
    ## Calculate 6 linear models: 0 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.02 minutes
    ## 13 cell-years remaining
    ## 
    ## Calculating models of with 5 input vectors.
    ## Define models: 0.01 minutes
    ## Calculate 2 linear models: 0 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.02 minutes
    ## 3 cell-years remaining
    ## 
    ## Calculating models of with 6 input vectors.
    ## Define models: 0.01 minutes
    ## Calculate 1 linear models: 0 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.02 minutes
    ## 
    ## Total Modeling Time: 0.1283067 minutes
    ## 
    ## Optimizing models: 0 minutes

    ## Warning: distinct() does not fully support columns of type `list`.
    ## List elements are compared by reference, see ?distinct for details.
    ## This affects the following columns:
    ## - `model`

``` r
# Generate predictions and uncertainty (and plot timeseries of each)                             
test.prediction <- predict_paleocar_models(models = test.vector,
                                           prediction.years = 600:1299)

test.prediction %>%
  ggplot(aes(x = year,
             y = Prediction)) +
  geom_ribbon(aes(ymin = Prediction - `PI Deviation`,
                  ymax = Prediction + `PI Deviation`),
              color = NA,
              fill = "dodgerblue") +
  geom_line(size = 0.2)
```

![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

##### `paleocar` reconstruction for multiple locations using the same set of predictors (in this case, tree-ring chronologies)

Running `paleocar` on a matrix of locations (`predictands`) will
generate reconstructions that select from the same set of predictors
(`chronologies`). The matrix must be formatted such that each location
is in a column, and each row is a year of data. Note that the number of
rows of the matrix must be the same as the number of years provided to
`calibration.years`.

``` r
# Extract a matrix of annualized climate data (all cells in the raster)
mvnp_prism.matrix <- mvnp_prism %>%
  raster::as.matrix() %>% 
  t()

test.matrix <- paleocar_models(predictands = mvnp_prism.matrix,
                               chronologies = itrdb,
                               calibration.years = 1924:1983,
                               prediction.years = 1:1985,
                               verbose = T)
```

    ## Calculating PaleoCAR models
    ## 
    ## Prepare data and calculate CAR scores: 0.1 minutes
    ## 
    ## Calculating models of with 1 input vectors.
    ## Define models: 0.03 minutes
    ## Calculate 9 linear models: 0.03 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.06 minutes
    ## 69264 cell-years remaining
    ## 
    ## Calculating models of with 2 input vectors.
    ## Define models: 0.03 minutes
    ## Calculate 24 linear models: 0.05 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.08 minutes
    ## 64246 cell-years remaining
    ## 
    ## Calculating models of with 3 input vectors.
    ## Define models: 0.03 minutes
    ## Calculate 34 linear models: 0.05 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.09 minutes
    ## 47452 cell-years remaining
    ## 
    ## Calculating models of with 4 input vectors.
    ## Define models: 0.03 minutes
    ## Calculate 36 linear models: 0.04 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.08 minutes
    ## 24085 cell-years remaining
    ## 
    ## Calculating models of with 5 input vectors.
    ## Define models: 0.02 minutes
    ## Calculate 27 linear models: 0.02 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.05 minutes
    ## 10839 cell-years remaining
    ## 
    ## Calculating models of with 6 input vectors.
    ## Define models: 0.02 minutes
    ## Calculate 12 linear models: 0.01 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.04 minutes
    ## 
    ## Total Modeling Time: 0.4005835 minutes
    ## 
    ## Optimizing models: 0.04 minutes

    ## Warning: distinct() does not fully support columns of type `list`.
    ## List elements are compared by reference, see ?distinct for details.
    ## This affects the following columns:
    ## - `model`

``` r
# Generate predictions and uncertainty (and plot location means in uncertainty)
test.prediction <- predict_paleocar_models(models = test.matrix,
                                           prediction.years = 600:1299)

test.prediction %>%
  dplyr::mutate(cell = as.factor(cell)) %>%
  dplyr::filter(cell %in% c(1,200,400,600)) %>%
  ggplot(aes(x = year,
             y = Prediction)) +
  geom_ribbon(aes(ymin = Prediction - `PI Deviation`,
                  ymax = Prediction + `PI Deviation`,
                  fill = cell),
              color = NA) +
  geom_line(size = 0.2) +
  facet_wrap(~cell, nrow = 2) +
  xlab("Year CE")
```

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

##### `paleocar` reconstruction over a grid

Paleocar can also be performed over a gridded climate dataset such as
PRISM, so long as it is a `RasterStack` or `RasterBrick` as defined in
the [`raster` package for
*R*](https://CRAN.R-project.org/package=raster). Results will be
returned in `RasterBrick` format.

``` r
# Print to show format
mvnp_prism
```

    ## class      : RasterStack 
    ## dimensions : 24, 26, 624, 60  (nrow, ncol, ncell, nlayers)
    ## resolution : 0.008333333, 0.008333333  (x, y)
    ## extent     : -108.5542, -108.3375, 37.15417, 37.35417  (xmin, xmax, ymin, ymax)
    ## crs        : +proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs 
    ## names      : X1924, X1925, X1926, X1927, X1928, X1929, X1930, X1931, X1932, X1933, X1934, X1935, X1936, X1937, X1938, ... 
    ## min values :   286,   360,   387,   499,   248,   434,   259,   289,   417,   239,   231,   324,   304,   377,   368, ... 
    ## max values :   498,   602,   615,   745,   417,   739,   437,   420,   690,   434,   364,   628,   588,   612,   720, ...

``` r
test.raster <- paleocar_models(predictands = mvnp_prism,
                               chronologies = itrdb,
                               calibration.years = 1924:1983,
                               prediction.years = 600:1299,
                               verbose = T)
```

    ## Calculating PaleoCAR models
    ## 
    ## Prepare data and calculate CAR scores: 0.1 minutes
    ## 
    ## Calculating models of with 1 input vectors.
    ## Define models: 0.01 minutes
    ## Calculate 5 linear models: 0.02 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.04 minutes
    ## 11856 cell-years remaining
    ## 
    ## Calculating models of with 2 input vectors.
    ## Define models: 0.01 minutes
    ## Calculate 12 linear models: 0.04 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.05 minutes
    ## 6838 cell-years remaining
    ## 
    ## Calculating models of with 3 input vectors.
    ## Define models: 0.01 minutes
    ## Calculate 16 linear models: 0.03 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.05 minutes
    ## 2361 cell-years remaining
    ## 
    ## Calculating models of with 4 input vectors.
    ## Define models: 0.01 minutes
    ## Calculate 14 linear models: 0.02 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.03 minutes
    ## 615 cell-years remaining
    ## 
    ## Calculating models of with 5 input vectors.
    ## Define models: 0.01 minutes
    ## Calculate 7 linear models: 0.01 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.02 minutes
    ## 105 cell-years remaining
    ## 
    ## Calculating models of with 6 input vectors.
    ## Define models: 0.01 minutes
    ## Calculate 3 linear models: 0 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.02 minutes
    ## 
    ## Total Modeling Time: 0.2213191 minutes
    ## 
    ## Optimizing models: 0.01 minutes

    ## Warning: distinct() does not fully support columns of type `list`.
    ## List elements are compared by reference, see ?distinct for details.
    ## This affects the following columns:
    ## - `model`

``` r
# Generate predictions and errors
test.raster.predictions <- predict_paleocar_models(models = test.raster,
                                                   prediction.years = 600:1299)

test.raster.predictions$Prediction %>%
  raster::mean() %>%
  raster::plot()
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

##### `paleocar()` convenience wrapper

The `paleocar()` convenience wrapper returns a list containing the
`models`, `reconstructions`, and `uncertainty`. The `paleocar()` method
also automatically saves the output of `predict_paleocar_models()` and
`errors_paleocar_models()`. Pass variables through this function to
other ones (e.g., `meanVar = "chained"`).

``` r
# Generate models and perform the reconstruction and error predictions.

mvnp_models <- paleocar_models(predictands = mvnp_prism,
                       label = "mvnp_prism",
                       chronologies = itrdb,
                       calibration.years = 1924:1983,
                       prediction.years = 600:1299,
                       out.dir = testDir,
                       force.redo = T,
                       verbose = T)
```

    ## Calculating PaleoCAR models
    ## 
    ## Prepare data and calculate CAR scores: 0.1 minutes
    ## 
    ## Calculating models of with 1 input vectors.
    ## Define models: 0.01 minutes
    ## Calculate 5 linear models: 0.02 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.04 minutes
    ## 11856 cell-years remaining
    ## 
    ## Calculating models of with 2 input vectors.
    ## Define models: 0.01 minutes
    ## Calculate 12 linear models: 0.03 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.05 minutes
    ## 6838 cell-years remaining
    ## 
    ## Calculating models of with 3 input vectors.
    ## Define models: 0.01 minutes
    ## Calculate 16 linear models: 0.03 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.05 minutes
    ## 2361 cell-years remaining
    ## 
    ## Calculating models of with 4 input vectors.
    ## Define models: 0.01 minutes
    ## Calculate 14 linear models: 0.02 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.03 minutes
    ## 615 cell-years remaining
    ## 
    ## Calculating models of with 5 input vectors.
    ## Define models: 0.01 minutes
    ## Calculate 7 linear models: 0.01 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.02 minutes
    ## 105 cell-years remaining
    ## 
    ## Calculating models of with 6 input vectors.
    ## Define models: 0.01 minutes
    ## Calculate 3 linear models: 0 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.02 minutes
    ## 
    ## Total Modeling Time: 0.2187503 minutes
    ## 
    ## Optimizing models: 0.01 minutes

    ## Warning: distinct() does not fully support columns of type `list`.
    ## List elements are compared by reference, see ?distinct for details.
    ## This affects the following columns:
    ## - `model`

``` r
mvnp_recon <- paleocar(models = mvnp_models,
                       predictands = mvnp_prism,
                       label = "mvnp_prism",
                       chronologies = itrdb,
                       calibration.years = 1924:1983,
                       prediction.years = 600:1299,
                       out.dir = testDir,
                       force.redo = T,
                       verbose = T)
```

    ## 
    ## Calculating all models
    ## Calculating PaleoCAR models
    ## 
    ## Prepare data and calculate CAR scores: 0.1 minutes
    ## 
    ## Calculating models of with 1 input vectors.
    ## Define models: 0.01 minutes
    ## Calculate 5 linear models: 0.02 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.04 minutes
    ## 11856 cell-years remaining
    ## 
    ## Calculating models of with 2 input vectors.
    ## Define models: 0.01 minutes
    ## Calculate 12 linear models: 0.03 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.05 minutes
    ## 6838 cell-years remaining
    ## 
    ## Calculating models of with 3 input vectors.
    ## Define models: 0.01 minutes
    ## Calculate 16 linear models: 0.03 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.05 minutes
    ## 2361 cell-years remaining
    ## 
    ## Calculating models of with 4 input vectors.
    ## Define models: 0.01 minutes
    ## Calculate 14 linear models: 0.02 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.03 minutes
    ## 615 cell-years remaining
    ## 
    ## Calculating models of with 5 input vectors.
    ## Define models: 0.01 minutes
    ## Calculate 7 linear models: 0.01 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.02 minutes
    ## 105 cell-years remaining
    ## 
    ## Calculating models of with 6 input vectors.
    ## Define models: 0.01 minutes
    ## Calculate 3 linear models: 0 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.02 minutes
    ## 
    ## Total Modeling Time: 0.2210963 minutes
    ## 
    ## Optimizing models: 0.01 minutes

    ## Warning: distinct() does not fully support columns of type `list`.
    ## List elements are compared by reference, see ?distinct for details.
    ## This affects the following columns:
    ## - `model`

    ## 
    ## Generating prediction

    ## 
    ## The entire reconstruction took 0.64 minutes

``` r
mvnp_recon <- paleocar(predictands = mvnp_prism,
                       label = "mvnp_prism",
                       chronologies = itrdb,
                       calibration.years = 1924:1983,
                       prediction.years = 600:1299,
                       out.dir = testDir,
                       force.redo = T,
                       verbose = T)
```

    ## 
    ## Calculating all models
    ## Calculating PaleoCAR models
    ## 
    ## Prepare data and calculate CAR scores: 0.1 minutes
    ## 
    ## Calculating models of with 1 input vectors.
    ## Define models: 0.01 minutes
    ## Calculate 5 linear models: 0.02 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.04 minutes
    ## 11856 cell-years remaining
    ## 
    ## Calculating models of with 2 input vectors.
    ## Define models: 0.01 minutes
    ## Calculate 12 linear models: 0.04 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.06 minutes
    ## 6838 cell-years remaining
    ## 
    ## Calculating models of with 3 input vectors.
    ## Define models: 0.02 minutes
    ## Calculate 16 linear models: 0.03 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.05 minutes
    ## 2361 cell-years remaining
    ## 
    ## Calculating models of with 4 input vectors.
    ## Define models: 0.01 minutes
    ## Calculate 14 linear models: 0.02 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.04 minutes
    ## 615 cell-years remaining
    ## 
    ## Calculating models of with 5 input vectors.
    ## Define models: 0.01 minutes
    ## Calculate 7 linear models: 0.01 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.03 minutes
    ## 105 cell-years remaining
    ## 
    ## Calculating models of with 6 input vectors.
    ## Define models: 0.01 minutes
    ## Calculate 3 linear models: 0 minutes
    ## Clean linear models: 0 minutes
    ## Total modeling time: 0.02 minutes
    ## 
    ## Total Modeling Time: 0.2312807 minutes
    ## 
    ## Optimizing models: 0.01 minutes

    ## Warning: distinct() does not fully support columns of type `list`.
    ## List elements are compared by reference, see ?distinct for details.
    ## This affects the following columns:
    ## - `model`

    ## 
    ## Generating prediction

    ## 
    ## The entire reconstruction took 0.61 minutes

``` r
# Examine the structure of the output
str(mvnp_recon, 
    max.level = 2)
```

    ## List of 2
    ##  $ models     :List of 5
    ##   ..$ models               :Classes 'tbl_df', 'tbl' and 'data.frame':    3115 obs. of  7 variables:
    ##   ..$ predictands          :Formal class 'RasterStack' [package "raster"] with 11 slots
    ##   ..$ predictor.matrix     : num [1:60, 1:120] 1.315 0.883 1.354 1.011 1.354 ...
    ##   .. ..- attr(*, "dimnames")=List of 2
    ##   ..$ reconstruction.matrix: num [1:700, 1:120] NA NA NA NA NA NA NA NA NA NA ...
    ##   .. ..- attr(*, "dimnames")=List of 2
    ##   ..$ carscores            :Classes 'tbl_df', 'tbl' and 'data.frame':    74880 obs. of  3 variables:
    ##  $ predictions:List of 4
    ##   ..$ cell        :Formal class 'RasterBrick' [package "raster"] with 12 slots
    ##   ..$ Prediction  :Formal class 'RasterBrick' [package "raster"] with 12 slots
    ##   ..$ CI Deviation:Formal class 'RasterBrick' [package "raster"] with 12 slots
    ##   ..$ PI Deviation:Formal class 'RasterBrick' [package "raster"] with 12 slots

You can quickly load a prior reconstruction by setting `force.redo =
FALSE`:

``` r
# Generate models and perform the reconstruction and error predictions.
mvnp_recon <- paleocar(predictands = mvnp_prism,
                       label = "mvnp_prism",
                       chronologies = itrdb,
                       calibration.years = 1924:1983,
                       prediction.years = 600:1299,
                       out.dir = testDir,
                       force.redo = F,
                       verbose = T)
```

    ## 
    ## Calculating all models
    ## 
    ## Generating prediction

    ## 
    ## The entire reconstruction took 0 minutes

#### Plot results

``` r
mvnp_recon$predictions$Prediction %>%
  raster::mean() %>%
  raster::plot()
```

![](README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
mvnp_recon$predictions$`CI Deviation` %>%
  raster::mean() %>%
  raster::plot()
```

![](README_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

``` r
mvnp_recon$predictions$`PI Deviation` %>%
  raster::mean() %>%
  raster::plot()
```

![](README_files/figure-gfm/unnamed-chunk-8-3.png)<!-- -->
