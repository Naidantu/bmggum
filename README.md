
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
#> Warning: There were 7 divergent transitions after warmup. See
#> http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> http://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> http://mc-stan.org/misc/warnings.html#tail-ess

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
#>  [1,]  1.02653892  0.1190288
#>  [2,] -1.49707999 -0.5484603
#>  [3,] -1.22685276 -0.5326736
#>  [4,]  0.73825429  0.8939866
#>  [5,]  0.01329808 -0.7303305
#>  [6,]  0.11594903 -0.7063715
#>  [7,]  0.87265448  0.3888945
#>  [8,]  0.08223168  0.3530410
#>  [9,]  0.64802058  0.7450873

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
#> [1,] -0.8738538 -1.5282223 -2.2845420
#> [2,] -1.5256997 -1.8672470 -1.0791201
#> [3,] -3.1769028 -0.8747523 -0.8697593
#> [4,] -2.3305834 -1.3352755 -1.0232087

# Extract the lambda estimates 
lambda <- Extract.BMGGUM(x=mod, pars='lambda')
# lambda[1,1] is the coefficient linking person covariate 1 to latent trait 1
# lambda[1,2] is the coefficient linking person covariate 1 to latent trait 2
lambda
#>                   mean    se_mean        sd       2.5%        50%    97.5%
#> lambda[1,1]  0.2715172 0.02619172 0.5021564 -0.7016386  0.2607421 1.280896
#> lambda[1,2] -0.2109062 0.07406314 0.7034274 -1.6469719 -0.2123634 1.187023
#>                 n_eff     Rhat
#> lambda[1,1] 367.57829 1.004650
#> lambda[1,2]  90.20576 1.042771

# Obtain model fit statistic waic 
waic <- Modfit.BMGGUM(x=mod, index='waic')
#> Warning: 
#> 14 (46.7%) p_waic estimates greater than 0.4. We recommend trying loo instead.
loo <- Modfit.BMGGUM(mod)
#> Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.

# Obtain the density plots for alpha
Bayesplot.BMGGUM(x=mod, pars='alpha', plot='density', inc_warmup=F)
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r

# Obtain item plots with ORCs for item 1, 2, 3
Itemplot.BMGGUM(x=mod, items = 1:3)
```

<img src="man/figures/README-example-2.png" width="100%" />
