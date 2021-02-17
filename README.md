
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bmggum

<!-- badges: start -->
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
4.  Five functions (i.e., bmggum( ), extract.bmggum( ), modfit.bmggum(
    ), bayesplot.bmggum( ), and itemplot.bmggum( )) are provided for
    model estimation, results extraction, model fit examination
    (e.g.,waic, loo, chisq/df), and plottings, respectively.

## Installation

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
theta <- extract.bmggum(x=mod, pars='theta')
# Turn the theta estimates into p*trait matrix where p equals sample size and trait equals the number of latent traits
theta <- theta[,1]
# nrow=trait
theta <- matrix(theta, nrow=2)  
theta <- t(theta)
# theta estimates in p*trait matrix format
theta
#>              [,1]        [,2]
#>  [1,]  1.01172010  0.03007375
#>  [2,] -1.49457843 -0.28825368
#>  [3,] -1.24305523 -0.56841647
#>  [4,]  0.70989404  1.00436767
#>  [5,] -0.05374441 -0.81460972
#>  [6,]  0.10233920 -0.82928930
#>  [7,]  0.80263701  0.21347193
#>  [8,]  0.11552137  0.54509935
#>  [9,]  0.67097662  0.89182374

# 3.2 Extract the tau estimates 
tau <- extract.bmggum(x=mod, pars='tau')
# Turn the tau estimates into I*(option-1) matrix where I equals the number of items and option equals the number of response options
tau <- tau[,1]
# nrow=option-1
tau <- matrix(tau, nrow=3)  
tau <- t(tau)
# tau estimates in I*(option-1) matrix format
tau
#>            [,1]       [,2]       [,3]
#> [1,] -0.8444732 -1.5242984 -2.2944532
#> [2,] -1.5878937 -1.8973009 -1.1197230
#> [3,] -3.2735301 -0.9010291 -0.8618411
#> [4,] -2.3371485 -1.3746858 -1.0575782

# 3.3 Extract the lambda estimates 
lambda <- extract.bmggum(x=mod, pars='lambda')
# lambda[1,1] is the coefficient linking person covariate 1 to latent trait 1
# lambda[1,2] is the coefficient linking person covariate 1 to latent trait 2
lambda
#>                   mean    se_mean        sd       2.5%        50%     97.5%
#> lambda[1,1]  0.2300054 0.02218539 0.4684616 -0.6407667  0.2297074 1.2042006
#> lambda[1,2] -0.3465347 0.04595152 0.6728208 -1.8855763 -0.2677931 0.8448643
#>                n_eff     Rhat
#> lambda[1,1] 445.8758 1.000741
#> lambda[1,2] 214.3873 1.003417

## Step 4: Obtain model fit statistics 
waic <- modfit.bmggum(x=mod, index='waic')
loo <- modfit.bmggum(mod)

## Step 5: Plottings
# 5.1 Obtain the density plots for alpha
bayesplot.bmggum(x=mod, pars='alpha', plot='density', inc_warmup=F)
```

<img src="man/figures/README-example-1.png" width="70%" />

``` r
## Step 6: Plotting observable response categories (ORCs) for items
# 6.1 Obtain item plots with ORCs for item 1, 2, 3
itemplot.bmggum(x=mod, items = 1:3)
```

<img src="man/figures/README-example-2.png" width="70%" />
