paleocar
========

[![Build Status](https://api.travis-ci.org/bocinsky/paleocar.png)](https://travis-ci.org/bocinsky/paleocar)

An *R* package implementing functions to perform spatio-temporal paleoclimate reconstruction from tree-rings using the CAR (Correlation Adjusted corRelation) approach of Zuber and Strimmer (2011). It is optimized for speed and memory use.

To install, use the following command in *R*:

`devtools::install_github("bocinsky/paleocar")`

This is based on the approach used in Bocinsky and Kohler (2014): 

Bocinsky, R. K. and Kohler, T. A. (2014). A 2,000-year reconstruction of the rain-fed maize agricultural niche in the US Southwest. *Nature Communications*, 5:5618. doi: [10.1038/ncomms6618](http://www.nature.com/ncomms/2014/141204/ncomms6618/full/ncomms6618.html).

The primary difference is that here model selection is performed by minimizing the corrected Akaike's Information Criterion.

A more recent reference would be Bocinsky et al. (2016):

Bocinsky, R. K., Rush, J., Kintigh, K. W., and Kohler, T. A. (2016). Exploration and exploitation in the macrohistory of the pre-Hispanic Pueblo Southwest. *Science Advances*, 2:[e1501532](http://advances.sciencemag.org/content/2/4/e1501532).