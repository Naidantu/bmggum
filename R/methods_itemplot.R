#' @title item plotting function including observable response categories (ORCs)
#' @description This function provides item plots including observable response categories plots.
#' @param x returned object
#' @param items The items to be plotted. The default is all the items.
#' @return Selected ORC plots for selected items
#' @examples
#' Data <- c(1,4,2,3)
#' Data <- matrix(Data,nrow = 2)
#' deli <- c(1,-1,2,1)
#' deli <- matrix(deli,nrow = 2)
#' ind <- c(1,2)
#' ind <- t(ind)
#' cova <- c(0.70, -1.25)
#' mod <- bmggum(GGUM.Data=Data,delindex=deli,trait=2,ind=ind,option=4,covariate=cova,iter=5,chains=1)
#' itemplot(mod, items=1)
#' @export
itemplot <- function(x, items=NULL){
  UseMethod("itemplot")
}


#' @export
#' @method itemplot bmggum 
itemplot.bmggum <- function(x, items=NULL){
  
  Pr.GGUM <- function(theta, alpha, delta, tau) {# ATTENTION! The first tau is always fixed to zero and has to be input into the tau vector.
    
    N.Cat<-dim(tau)[2]
    N.Resp<-ncol(theta)
    N.Items<-length(alpha)
    Prob<-matrix(NA,nrow=N.Resp,ncol = N.Items*N.Cat)
    Prob.New<-matrix(NA,nrow=N.Resp,ncol = N.Items*N.Cat)
    
    for (i in 1:N.Resp){
      for (q in 1:N.Items){
        for (k in 1:N.Cat){
          Prob[i,((q-1)*N.Cat+k)]<-(exp(alpha[q]*((k-1)*(theta[q,i]-delta[q])-sum(tau[q,1:k])))+exp(alpha[q]*((2*N.Cat-k)*(theta[q,i]-delta[q])-sum(tau[q,1:k]))))
        }
        Prob.New[i,((q-1)*N.Cat+1):(q*N.Cat)]<- Prob[i,((q-1)*N.Cat+1):(q*N.Cat)]/sum(Prob[i,((q-1)*N.Cat+1):(q*N.Cat)])
      }
    }
    return(Prob.New)
  }
  
  alpha <- extract(x, pars='alpha')
  delta <- extract(x, pars='delta')
  alpha <- alpha[,1]
  delta <- delta[,1]
  tau <- extract(x, pars='tau')
  tau <- tau[,1]
  C <- length(tau)/length(alpha)  #response options = C+1
  I <- length(alpha)
  tau <- matrix(tau, nrow=C)
  tau <- t(tau)
  tau <- cbind(rep(0,I), tau)
  th.lims <- cbind(delta - 4, delta + 4)
  th.vals <- t(apply(th.lims, 1, function(vec) {seq(vec[1], vec[2], length.out = 160)}))
  
  Prob <- Pr.GGUM(th.vals, alpha, delta, tau)
  
  if (is.null(items)) {
    I.plot <- 1:I
  } else {
    I.plot <- items
  }
  
  Prob <- Prob[, ((min(I.plot)-1)*(C+1)+1):max(I.plot*(C+1))]
  Prob <- array(Prob, dim=c(160, C+1, length(I.plot)))
  
  d <- dim(Prob)
  th.vals <- t(th.vals)
  theta <- matrix(NA, nrow = 160*(C+1)*length(I.plot), ncol = 1)
  for (i in 1:length(I.plot)){
    it <- I.plot[i]
    for (j in 1:(C+1)){
      theta[((160*(j-1)+1)+(i-1)*160*(C+1)):(160*j+(i-1)*160*(C+1)), 1] <- th.vals[,it]
    }
  }
  probability <- NULL
  options <- NULL
  plotdata <- data.frame(theta = theta, probability = c(Prob), plots = gl(d[3], prod(d[-3])), options = gl(d[2], d[1]))
  label = character(length(I.plot))
  for (i in 1:length(I.plot)){
    label[i] <- paste0("Item ", I.plot[i])
  }
  plotdata$plots <- factor(plotdata$plots, labels=label)
  return(ggplot2::ggplot(plotdata, ggplot2::aes(x=theta, y=probability, col = options)) + ggplot2::geom_line() + ggplot2::facet_wrap(~plots) + ggplot2::labs(x="theta"))
  
}


