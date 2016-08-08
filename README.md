A fast algorithm for S-estimators
================
Matias Salibian
2016-08-08

A fast algorithm for S-estimators for linear regression
-------------------------------------------------------

This repository contains stand-alone `R` code implementing the fast algorithm for S-regression estimators proposed in [Salibian-Barrera, M. and Yohai, V.J. (2006).](http://dx.doi.org/10.1198/106186006X113629). The algorithm is also implemented in the function `robustbase::lmrob` that computes S- and MM-estimators. See the R package [robustbase](http://cran.r-project.org/web/packages/robustbase/index.html) for more information.

Here is a simple example on how to use this code. We first read the function:

``` r
source('fasts.R')
```

We will apply it to the well-known Boston Housing data set.

``` r
data(Boston, package='MASS')
x <- model.matrix(medv ~ ., data=Boston)
y <- Boston$medv
sest <- fast.s(x=x, y=y, int=0, N=500, k=2, best.r=5, seed=123)
```

The arguments of `fast.s` are:

-   x: design matrix, covariates for each observation in rows
-   y: vector of responses
-   int: if equal to 1 a column of ones is added to x. This corresponds to including an intercept in the fitted model.
-   N: number of random sub-samples to be generated
-   k: number of refining iterations in each subsample. k = 0 means "raw-subsampling"
-   b: right hand side of the M-scale equation. Defaults to 0.50 yielding a highly-robust estimator.
-   cc: tuning constant corresponding to the argument `b` above to ensure consistency.
-   best.r: number of partially refined candidates to retain. These are finally iterated until convergence.
-   seed: seed used in the generation of the random sub-samples.

We can compare the S- and LS-estimators:

``` r
lmest <- lm(medv ~ ., data=Boston)
round(cbind(sest[[1]], coef(lmest)), 3)
```

    ##               [,1]    [,2]
    ## (Intercept)  8.116  36.459
    ## crim        -0.142  -0.108
    ## zn           0.036   0.046
    ## indus        0.034   0.021
    ## chas         2.346   2.687
    ## nox         -3.369 -17.767
    ## rm           4.844   3.810
    ## age         -0.044   0.001
    ## dis         -0.717  -1.476
    ## rad          0.064   0.306
    ## tax         -0.006  -0.012
    ## ptratio     -0.540  -0.953
    ## black        0.013   0.009
    ## lstat       -0.247  -0.525

The estimated residual scales are

``` r
round(c(sest[[2]], summary(lmest)$sigma), 2)
```

    ## [1] 2.72 4.75

I will later add here an illustration of the better predictions obtained with the robust estimator.
