
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BMGGUM

<!-- badges: start -->

<!-- badges: end -->

The goal of BMGGUM is to estimate Multidimensional Generalized Graded
Unfolding Model (MGGUM) using Bayesian method. Specifically,the R
package **rstan** that utilizes the Hamiltonian Monte Carlo sampling
algorithm was used for estimation. Below are some important features of
the BMGGUM package:

1.  Allows users to incorporate person covariates (e.g., age, gender,
    education) into the estimation to improve estimation accuracy.
2.  Automatically deals with missing data in a way similar to how full
    information maximum likelihood handles missing data.
3.  Allows users to estimate the **multidimensional version** of three
    unfolding models that are available in the software GGUM2004
    (Roberts, Fang, Cui, & Wang, 2006).
      - UM8: The Generalized Graded Unfolding Model (GGUM).
      - UM4: The Partial Credit Unfolding Model, which is the GGUM with
        all alphas constrained to 1.
      - UM7: The Generalized Rating Scale Unfolding Model, which is the
        GGUM with equal taus across items.
4.  Five functions (i.e., BMGGUM(), Extract.BMGGUM(), Modfit.BMGGUM(),
    Bayesplot.BMGGUM(), and Itemplot.BMGGUM()) are provided for model
    estimation, results extraction, model fit examination (e.g.,WAIC,
    loo), and plottings, respectively.

## Installation

You can install the development version of BMGGUM from GitHub:

``` r
devtools::install_github("Naidantu/BMGGUM")
```

## Example

This is a basic example which shows you how to prepare data, fit the
model, extract and plot results.

``` r
library(BMGGUM)

## basic example code
# Response data
GGUM.Data <- c(1,4,4,1,1,1,1,1,1,1,4,1,1,3,1,1,NA,2,NA,3,2,2,2,1,3,2,NA,2,1,1,2,1,NA,NA,NA,1,3,NA,1,2)
GGUM.Data <- matrix(GGUM.Data,nrow = 10)

# delindex
delindex <- c(1,-1,2,1,3,-1,4,1)
delindex <- matrix(delindex,nrow = 2)

# ind
ind <- c(1,1,2,2)
ind <- t(ind)

# covariate
covariate <- c(0.70, -1.25, 0.48, -0.47, 0.86, 1.25, 1.17, -1.35, -0.84, -0.55)

# Fit the MGGUM model
mod <- BMGGUM(GGUM.Data=GGUM.Data, delindex=delindex, trait=2, ind=ind, option=4, model="UM8", covariate=covariate)
#> [1] "Case 9 was deleted because they endorse the same response option across all items"

# Extract the theta estimates 
theta <- Extract.BMGGUM(x=mod, pars='theta')
# Turn the theta estimates into p*trait matrix where p equals sample size and trait equals the number of latent traits
theta <- theta[,1]
# nrow=trait
theta <- matrix(theta, nrow=2)  
theta <- t(theta)
# theta estimates in p*trait matrix format
theta
#>              [,1]       [,2]
#>  [1,]  1.02027771  0.1112450
#>  [2,] -1.50303636 -0.4757993
#>  [3,] -1.23103651 -0.5853762
#>  [4,]  0.69315111  0.9592051
#>  [5,] -0.03894941 -0.7695685
#>  [6,]  0.15619541 -0.6873740
#>  [7,]  0.84861339  0.3318221
#>  [8,]  0.12068104  0.4478592
#>  [9,]  0.66487328  0.8901462

# Extract the tau estimates 
tau <- Extract.BMGGUM(x=mod, pars='tau')
# Turn the tau estimates into I*(option-1) matrix where I equals the number of items and option equals the number of response options
tau <- tau[,1]
# nrow=option-1
tau <- matrix(tau, nrow=3)  
tau <- t(tau)
# tau estimates in I*(option-1) matrix format
tau
#>            [,1]       [,2]       [,3]
#> [1,] -0.8559964 -1.4779004 -2.3206627
#> [2,] -1.5405565 -1.8804769 -1.0592084
#> [3,] -3.1956647 -0.8870041 -0.8631355
#> [4,] -2.2876189 -1.3062246 -1.0204106

# Extract the lambda estimates 
lambda <- Extract.BMGGUM(x=mod, pars='lambda')
# lambda[1,1] is the coefficient linking person covariate 1 to latent trait 1
# lambda[1,2] is the coefficient linking person covariate 1 to latent trait 2
lambda
#>                   mean    se_mean        sd       2.5%        50%     97.5%
#> lambda[1,1]  0.2563650 0.03887061 0.5254527 -0.6997604  0.2333346 1.4214331
#> lambda[1,2] -0.2542906 0.02739855 0.6164561 -1.5319360 -0.2272128 0.8995756
#>                n_eff     Rhat
#> lambda[1,1] 182.7362 1.023394
#> lambda[1,2] 506.2316 1.000784

# Obtain model fit statistic waic 
waic <- Modfit.BMGGUM(x=mod, index='waic')
loo <- Modfit.BMGGUM(mod)

# Obtain the density plots for alpha
Bayesplot.BMGGUM(x=mod, pars='alpha', plot='density', inc_warmup=F)
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r

# Obtain item plots with ORCs for item 1, 2, 3
Itemplot.BMGGUM(x=mod, items = 1:3)
```

<img src="man/figures/README-example-2.png" width="100%" />
