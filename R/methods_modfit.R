#' @title Model fit
#' @description This function provides model fit statistics.
#' @param x returned object
#' @param index Model fit indices. They can be "waic", which is the widely applicable information criterion, "loo", which is the leave-one-out cross-validation, or "chisq.df", which is the adjusted chi-square degrees of freedom ratios for each trait separately that were introduced by Drasgow et al. (1995). The default is loo. Note that chisq.df can only be computed when the sample size is large. See documentation for loo and GGUM for more details.
#' @return Selected model fit statistics
#' @examples
#' Data <- c(1,4,2,3)
#' Data <- matrix(Data,nrow = 2)
#' deli <- c(1,-1,2,1)
#' deli <- matrix(deli,nrow = 2)
#' ind <- c(1,2)
#' ind <- t(ind)
#' cova <- c(0.70, -1.25)
#' mod <- bmggum(GGUM.Data=Data,delindex=deli,trait=2,ind=ind,option=4,covariate=cova,iter=5,chains=1)
#' waic <- modfit(mod, 'waic')
#' @export
modfit <- function(x, index="loo"){
  UseMethod("modfit")
}


#' @export
#' @method modfit bmggum 
modfit.bmggum <- function(x, index="loo"){
  
  x1 <- extract(x, 'fit')
  log_lik=loo::extract_log_lik(x1, merge_chains = FALSE)
  rel_n_eff=loo::relative_eff(exp(log_lik))
  
  if (index=="chisq.df"){
    dimension <- extract(x, 'dimension')
    data <- extract(x, 'data')-1
    tau <- extract(x, 'tau')
    tau <- tau[,1]
    alpha <- extract(x, pars='alpha')
    alpha <- alpha[,1]
    delta <- extract(x, pars='delta')
    delta <- delta[,1]
    C <- length(tau)/length(delta)  #response options = C+1
    I <- length(alpha)
    tau <- matrix(tau, nrow=C)
    tau <- t(tau)
    rtau <- -1*tau
    if (C==1){
      taus <- data.frame(tau, rep(0,I), rtau)
    }
    else{
      taus <- data.frame(tau, rep(0,I), rtau[,order(ncol(rtau):1)])
    }
    chisq.df <- NULL
    for (i in 1:max(dimension)){
      Data <- NULL
      Alpha <- NULL
      Delta <- NULL
      Taus <- NULL
      list <- NULL
      for (j in 1:I){
        if (dimension[j]==i){
          Data <- cbind(Data, data[,j])
          Alpha <- c(Alpha, alpha[j])
          Delta <- c(Delta, delta[j])
          Taus <- rbind(Taus, taus[j,])
        }
      }
      list <- list(data=Data, C=C, model="GGUM", alpha=Alpha, delta=Delta, taus=Taus)
      class(list) <- "GGUM"
      chisq.df <- c(chisq.df, GGUM::MODFIT(list))
    }
  }
  
  ret <- switch(index,
                waic=loo::waic(log_lik),
                loo=loo::loo(log_lik, r_eff = rel_n_eff, cores = 2),
                chisq.df=chisq.df)
  ret
}
