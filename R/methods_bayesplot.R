#' @title bayesian convergence diagnosis plotting function
#' @description This function provides plots including density plots, trace plots, and auto-correlation plots to aid model convergence diagnosis.
#' @param x returned object
#' @param pars Names of plotted parameters. They can be "theta", "alpha", "delta", "tau", "cor", "lambda", or a subset of parameters. See vignette for bmggum for more details.
#' @param plot Types of plots.They can be "density", "trace", or "autocorrelation".
#' @param inc_warmup Whether to include warmup iterations or not when plotting. The default is FALSE.
#' @return Selected plots for selected parameters
#' @examples
#' Data <- c(1,4,2,3)
#' Data <- matrix(Data,nrow = 2)
#' deli <- c(1,-1,2,1)
#' deli <- matrix(deli,nrow = 2)
#' ind <- c(1,2)
#' ind <- t(ind)
#' cova <- c(0.70, -1.25)
#' mod <- bmggum(GGUM.Data=Data,delindex=deli,trait=2,ind=ind,option=4,covariate=cova,iter=5,chains=1)
#' bayesplot(mod, 'alpha', 'density', inc_warmup=FALSE)
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
