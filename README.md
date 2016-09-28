# paleocar
An *R* package implementing functions to perform spatio-temporal paleoclimate reconstruction from tree-rings using the CAR (Correlation Adjusted corRelation) approach of Zuber and Strimmer (2011). It is optimized for speed and memory use.

To install, use the following command in *R*:

`devtools::install_github("bocinsky/paleocar")`

This is based on the approach used in Bocinsky and Kohler (2014): Bocinsky, R. K. and Kohler, T. A. (2014). The primary difference is that here model selection is performed by minimizing the corrected Akaike's Information Criterion.

A 2,000-year reconstruction of the rain-fed maize agricultural niche in the US Southwest. *Nature Communications*, 5:5618. doi: [10.1038/ncomms6618](http://www.nature.com/ncomms/2014/141204/ncomms6618/full/ncomms6618.html).