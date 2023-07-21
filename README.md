<!-- badges: start -->
[![R-CMD-check](https://github.com/Mengbo-Li/protDP/workflows/R-CMD-check/badge.svg)](https://github.com/Mengbo-Li/protDP/actions)
[![Codecov](https://codecov.io/gh/Mengbo-Li/protDP/branch/main/graph/badge.svg)](https://app.codecov.io/gh/Mengbo-Li/protDP?branch=main)
<!-- badges: end -->




# Missing values are informative in label-free shotgun proteomics data

The relationship between missingness and intensity is evaluated in label-free
shotgun proteomics data. Using a series of publicly available datasets, we first 
empirically investigate the missingness-intensity relationship with observed 
data. We also propose a model for detection probability which describes
the probability of an observation being detected given its underlying intensity
for label-free shotgun proteomics data. 

## Reference
Li, M. and Smyth, G.K., 2023. Neither random nor censored: estimating intensity-dependent probabilities for missing values in label-free proteomics. Bioinformatics, 39(5), [btad200](https://doi.org/10.1093/bioinformatics/btad200). 

## Installation 

To install the package, use the following script in R:

```
# install.packages("devtools")
devtools::install_github("Mengbo-Li/protDP")
```

## Quick start

If you are interested in investigating the relationship between intensity and detection/missing values on your own proteomics dataset, then you can try the following as a quick start: 

```
dpcfit <- dpc(dat)
plotDPC(dpcfit)
```
where `dat` is the log2-intensity matrix, with rows being precursors/proteins and columns being the samples, and `dat` contains some missing values as `NA`. 

The `plotDPC()` function then visualises the detection probability curve (DPC) from which you can inspect whether missingness is dependent on the underlying intensity values on your dataset. This tells you whether missingness is missing not at random (MNAR). 


## More examples

See data examples at https://mengbo-li.github.io/protDP/articles/protDP.html. 

## Contact

Open an [**issue**](https://github.com/Mengbo-Li/protDP/issues) should there be
any questions. 
