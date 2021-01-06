
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
mod <- BMGGUM(GGUM.Data=GGUM.Data, delindex=delindex, trait=2, ind=ind, option=4, model="UM8", covariate=covariate)
#> [1] "Case 9 was deleted because they endorse the same response option across all items"

## Step 3: Extract the estimated results 
# 3.1 Extract the theta estimates 
theta <- Extract.BMGGUM(x=mod, pars='theta')
# Turn the theta estimates into p*trait matrix where p equals sample size and trait equals the number of latent traits
theta <- theta[,1]
# nrow=trait
theta <- matrix(theta, nrow=2)  
theta <- t(theta)
# theta estimates in p*trait matrix format
theta
#>              [,1]       [,2]
#>  [1,]  1.00755202  0.1051978
#>  [2,] -1.48304663 -0.3815527
#>  [3,] -1.24215031 -0.5797086
#>  [4,]  0.66723737  0.9480621
#>  [5,]  0.09211222 -0.7688336
#>  [6,]  0.25556022 -0.7676078
#>  [7,]  0.87850501  0.3453732
#>  [8,] -0.04140219  0.4592622
#>  [9,]  0.65549867  0.8783695

# 3.2 Extract the tau estimates 
tau <- Extract.BMGGUM(x=mod, pars='tau')
# Turn the tau estimates into I*(option-1) matrix where I equals the number of items and option equals the number of response options
tau <- tau[,1]
# nrow=option-1
tau <- matrix(tau, nrow=3)  
tau <- t(tau)
# tau estimates in I*(option-1) matrix format
tau
#>            [,1]      [,2]       [,3]
#> [1,] -0.8295094 -1.498559 -2.3144571
#> [2,] -1.5646050 -1.889293 -1.0866118
#> [3,] -3.2392564 -0.884879 -0.8864353
#> [4,] -2.3187464 -1.317304 -1.0316517

# 3.3 Extract the lambda estimates 
lambda <- Extract.BMGGUM(x=mod, pars='lambda')
# lambda[1,1] is the coefficient linking person covariate 1 to latent trait 1
# lambda[1,2] is the coefficient linking person covariate 1 to latent trait 2
lambda
#>                   mean    se_mean        sd       2.5%        50%     97.5%
#> lambda[1,1]  0.3133386 0.05326958 0.5466074 -0.6521527  0.2707874 1.5867761
#> lambda[1,2] -0.2758883 0.03527324 0.6619620 -1.8982070 -0.2321969 0.8771598
#>                n_eff     Rhat
#> lambda[1,1] 105.2913 1.017840
#> lambda[1,2] 352.1887 1.005127

## Step 4: Obtain model fit statistics 
waic <- Modfit.BMGGUM(x=mod, index='waic')
loo <- Modfit.BMGGUM(mod)

## Step 5: Plottings
# 5.1 Obtain the density plots for alpha
Bayesplot.BMGGUM(x=mod, pars='alpha', plot='density', inc_warmup=F)
```

<img src="man/figures/README-example-1.png" width="70%" />

``` r

## Step 6: Plotting observable response categories (ORCs) for items
# 6.1 Obtain item plots with ORCs for item 1, 2, 3
Itemplot.BMGGUM(x=mod, items = 1:3)
```

<img src="man/figures/README-example-2.png" width="70%" />
