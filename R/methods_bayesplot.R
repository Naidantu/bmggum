#' @title bayesian convergence diagnosis plotting function
#' @description This function provides plots including density plots, trace plots, and auto-correlation plots to aid model convergence diagnosis.
#' @param x returned object
#' @param pars Names of plotted parameters. They can be "theta", "alpha", "delta", "tau", "cor", "lambda", or a subset of parameters. See vignette for bmggum for more details.
#' @param plot Types of plots.They can be "density", "trace", or "autocorrelation".
#' @param inc_warmup Whether to include warmup iterations or not when plotting. The default is FALSE.
#' @return Selected plots for selected parameters
#' @examples
#' \donttest{
#' Data <- c(1,4,4,1,1,1,1,1,1,1,4,1,1,3,1,1,NA,2,NA,3,2,2,2,1,3,2,NA,2,1,1)
#' Data <- matrix(Data,nrow = 10)
#' deli <- c(1,-1,2,1,3,-1)
#' deli <- matrix(deli,nrow = 2)
#' ind <- c(1,1,2)
#' ind <- t(ind)
#' cova <- c(0.70, -1.25, 0.48, -0.47, 0.86, 1.25, 1.17, -1.35, -0.84, -0.55)
#' mod <- bmggum(GGUM.Data=Data, delindex=deli, trait=2, ind=ind, option=4, covariate=cova)
#' bayesplot(mod, 'alpha', 'density', inc_warmup=TRUE)
#' bayesplot(mod, 'delta', 'trace', inc_warmup=FALSE)
#' bayesplot(mod, 'tau', 'autocorrelation', inc_warmup=TRUE)}
#' @export
bayesplot <- function(x, pars, plot, inc_warmup=FALSE){
  UseMethod("bayesplot")
}


#' @export
#' @method bayesplot bmggum 
bayesplot.bmggum <- function(x, pars, plot, inc_warmup=FALSE){
  
  x <- extract(x, 'fit')
  if (pars=="cor"){
    pars="Cor"
  }
  if (plot=="trace"){
    
    ret <- rstan::stan_trace(x, pars, inc_warmup = inc_warmup)
    
  }
  if (plot=="density"){
    
    ret <- rstan::stan_dens(x, pars, inc_warmup = inc_warmup)
    
  }
  if (plot=="autocorrelation"){
    
    ret <- rstan::stan_ac(x, pars, inc_warmup = inc_warmup)
    
  }
  ret
}
