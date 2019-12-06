MTSYS
=====

[![Travis-CI Build Status](https://travis-ci.org/okayaa/MTSYS.svg?branch=master)](https://travis-ci.org/okayaa/MTSYS) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/MTSYS)](https://CRAN.R-project.org/package=MTSYS) [![Coverage Status](https://img.shields.io/codecov/c/github/okayaa/MTSYS/master.svg)](https://codecov.io/github/okayaa/MTSYS?branch=master)

MTSYS provides a collection of multivariate analysis methods in 
Mahalanobis-Taguchi System (MTS), which was developed for the field of quality 
engineering. MTS consists of two families depending on their purpose. One is a 
family of Mahalanobis-Taguchi (MT) methods (in the broad sense) for diagnosis 
and the other is a family of Taguchi (T) methods for forecasting.

Overview
--------

The following methods are implemented.

### A family of MT methods
-   MT method
-   MTA method
-   RT method

### A family of T methods
-   T(1) method
-   Ta method
-   Tb method

For details, see the following referenses.

Installation
------------

Install the release version from CRAN:

``` r
install.packages("MTSYS")
```

Or the development version from github

``` r
# install.packages("devtools")
devtools::install_github("okayaa/MTSYS")
```

Example
-------

``` r
library(MTSYS)

# 40 data for versicolor in the iris dataset
iris_versicolor <- iris[61:100, -5]

unit_space_MT <- MT(unit_space_data = iris_versicolor)

# 10 data for each kind (setosa, versicolor, virginica) in the iris dataset
iris_test <- iris[c(1:10, 51:60, 101:110), -5]

diagnosis_MT <- diagnosis(unit_space = unit_space_MT, newdata = iris_test, 
                          threshold = 4)

(diagnosis_MT$le_threshold)
#>     1     2     3     4     5     6     7     8     9    10
#> FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#>   51    52    53    54    55    56    57    58    59    60   
#> TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
#>  101   102   103   104   105   106   107   108   109   110 
#> TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
```

References
----------

-   Woodall, W. H., Koudelik, R., Tsui, K. L., Kim, S. B., Stoumbos, Z. G., and
      Carvounis, C. P. (2003) A review and analysis of the Mahalanobis-Taguchi 
      system. *Technometrics, 45*(1), 1-15. \<[doi:10.1198/004017002188618626](http://dx.doi.org/10.1198/004017002188618626)\>
-   Kawada, H., and Nagata, Y. (2015) An application of a generalized inverse 
      regression estimator to Taguchi's T-Method. *Total Quality Science, 1*(1), 
      12-21. \<[doi:10.17929/tqs.1.12](http://dx.doi.org/10.17929/tqs.1.12)\>
