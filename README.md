
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bmggum

<!-- badges: start -->

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/bmggum)](https://cran.r-project.org/package=bmggum)
[![Total
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/bmggum)](https://cranlogs.r-pkg.org/badges/grand-total/bmggum "Total downloads")
<!-- badges: end -->

The goal of bmggum is to estimate Multidimensional Generalized Graded
Unfolding Model (MGGUM) using Bayesian method. Specifically,the R
package **rstan** that utilizes the Hamiltonian Monte Carlo sampling
algorithm was used for estimation. Below are some important features of
the bmggum package:

1.  Allows users to incorporate person covariates (e.g., age, gender,
    education) into the estimation to improve estimation accuracy.
2.  Automatically deals with missing data in a way similar to how full
    information maximum likelihood handles missing data.
3.  Allows users to estimate the **multidimensional version** of three
    unfolding models that are available in the software GGUM2004
    (Roberts, Fang, Cui, & Wang, 2006).
    -   UM8: The Generalized Graded Unfolding Model (GGUM).
    -   UM4: The Partial Credit Unfolding Model, which is the GGUM with
        all alphas constrained to 1.
    -   UM7: The Generalized Rating Scale Unfolding Model, which is the
        GGUM with equal taus across items.
4.  Five functions (i.e., bmggum( ), extract( ), modfit( ), bayesplot(
    ), and itemplot( )) are provided for model estimation, results
    extraction, model fit examination (e.g.,waic, loo, chisq/df), and
    plottings, respectively.

## Installation

You can install bmggum from CRAN:

``` r
install.packages("bmggum")
```

You can install the development version of bmggum from GitHub:

``` r
devtools::install_github("Naidantu/bmggum")
```

## Example

This is a basic example that shows you how to prepare data, fit the
model, extract and plot results.

``` r
library(bmggum)

## basic example code
## Step 1: Input data
# 1.1 Response data in wide format
GGUM.Data <- c(1,4,4,1,1,1,1,1,1,1,4,1,1,3,1,1,NA,2,NA,3,2,2,2,1,3,2,NA,2,1,1,2,1,NA,NA,NA,1,3,NA,1,2)
GGUM.Data <- matrix(GGUM.Data,nrow = 10)

# 1.2 A two-row data matrix: the first row is the item number (1,2,3,4...); the second row indicates the signs of delta for each item (-1,0,1,...). For items that have negative deltas for sure, "-1" should be assigned; for items that have positive deltas, "1" should be assigned; for items whose deltas may be either positive or negative (e.g., intermediate items), "0" should assigned. We recommend at least two positive and two negative items per trait for better estimation.
delindex <- c(1,-1,2,1,3,-1,4,1)
delindex <- matrix(delindex,nrow = 2)

# 1.3 A row vector mapping each item to each trait. For example, c(1,1,1,2,2,2) means that the first 3 items belong to trait 1 and the last 3 items belong to trait 2.
ind <- c(1,1,2,2)
ind <- t(ind)

# 1.4 An p*c person covariate matrix where p equals sample size and c equals the number of covariates. The default is NULL, meaning no person covariate.
covariate <- c(0.70, -1.25, 0.48, -0.47, 0.86, 1.25, 1.17, -1.35, -0.84, -0.55)

## Step 2: Fit the MGGUM model
mod <- bmggum(GGUM.Data=GGUM.Data, delindex=delindex, trait=2, ind=ind, option=4, model="UM8", covariate=covariate)
#> [1] "Case 9 was deleted because they endorse the same response option across all items"

## Step 3: Extract the estimated results 
# 3.1 Extract the theta estimates 
theta <- extract(x=mod, pars='theta')
# Turn the theta estimates into p*trait matrix where p equals sample size and trait equals the number of latent traits
theta <- theta[,1]
# nrow=trait
theta <- matrix(theta, nrow=2)  
theta <- t(theta)
# theta estimates in p*trait matrix format
theta
#>              [,1]        [,2]
#>  [1,]  1.01537640  0.05016495
#>  [2,] -1.53290020 -0.34579264
#>  [3,] -1.22716264 -0.55977442
#>  [4,]  0.73438219  1.01150777
#>  [5,] -0.02289508 -0.83776626
#>  [6,]  0.11009091 -0.86844999
#>  [7,]  0.85762438  0.26520432
#>  [8,]  0.11996516  0.53755372
#>  [9,]  0.63673648  0.86268826

# 3.2 Extract the tau estimates 
tau <- extract(x=mod, pars='tau')
# Turn the tau estimates into I*(option-1) matrix where I equals the number of items and option equals the number of response options
tau <- tau[,1]
# nrow=option-1
tau <- matrix(tau, nrow=3)  
tau <- t(tau)
# tau estimates in I*(option-1) matrix format
tau
#>            [,1]       [,2]       [,3]
#> [1,] -0.8412185 -1.5671934 -2.2907142
#> [2,] -1.5971053 -1.9075990 -1.1402262
#> [3,] -3.2093114 -0.8936897 -0.8808858
#> [4,] -2.3441511 -1.3209242 -1.0455207

# 3.3 Extract the lambda estimates 
lambda <- extract(x=mod, pars='lambda')
# lambda[1,1] is the coefficient linking person covariate 1 to latent trait 1
# lambda[1,2] is the coefficient linking person covariate 1 to latent trait 2
lambda
#>                   mean    se_mean        sd       2.5%        50%    97.5%
#> lambda[1,1]  0.2474065 0.02845930 0.5054661 -0.7168954  0.2467913 1.265200
#> lambda[1,2] -0.3297602 0.05024748 0.7360003 -2.1310694 -0.2765893 1.003908
#>                n_eff     Rhat
#> lambda[1,1] 315.4537 1.005384
#> lambda[1,2] 214.5494 1.006322

## Step 4: Obtain model fit statistics 
waic <- modfit(x=mod, index='waic')
loo <- modfit(mod)

## Step 5: Plottings
# 5.1 Obtain the density plots for alpha
bayesplot(x=mod, pars='alpha', plot='density', inc_warmup=FALSE)
```

<img src="man/figures/README-example-1.png" width="70%" />

``` r
## Step 6: Plotting observable response categories (ORCs) for items
# 6.1 Obtain item plots with ORCs for item 1, 2, 3
itemplot(x=mod, items = 1:3)
```

<img src="man/figures/README-example-2.png" width="70%" />
