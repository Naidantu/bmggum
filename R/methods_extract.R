#' @title results extraction
#' @description This function extracts estimation results.
#' @param x returned object
#' @param pars Names of extracted parameters. They can be "theta" (person trait estimates), "alpha" (item discrimination parameters), "delta" (item location parameters), "tau" (item threshold parameters), "cor" (correlations among latent traits), "lambda" (regression coefficients linking person covariates to latent traits), "data" (GGUM.Data after deleting respondents who endorse the same response options across all items), "fit" (the stanfit object), and "dimension" (the input row vector mapping each item to each trait). Note that when the model is UM4 in which alpha is fixed to 1, the extracted alpha is a n*1 matrix where n equals to the number of items.
#' @return Selected results output
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
#' theta <- extract(mod, 'theta')
#' alpha <- extract(mod, 'alpha')
#' cor <- extract(mod, 'cor')
#' data <- extract(mod, 'data')}
#' @export
extract <- function(x, pars){
  UseMethod("extract")
}


#' @export
#' @method extract bmggum 
extract.bmggum <- function(x, pars){
  
  ret <- switch(pars,
                theta=x[["Theta.est"]],
                alpha=as.matrix(x[["Alpha.est"]]),
                delta=x[["Delta.est"]],
                tau=x[["Tau.est"]],
                cor=x[["Cor.est"]],
                lambda=x[["Lamda.est"]],
                data=x[["Data"]],
                fit=x[["Fit"]],
                dimension=x[["Dimension"]])
  
  ret
}
