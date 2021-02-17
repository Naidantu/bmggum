#' @title Bayesian Multidimensional Generalized Graded Unfolding Model (bmggum)
#' @description This function implements full Bayesian estimation of Multidimensional Generalized Graded Unfolding Model (MGGUM) using rstan
#' @param GGUM.Data Response data in wide format
#' @param delindex A two-row data matrix: the first row is the item number (1, 2, 3, 4...); the second row indicates the signs of delta for each item (-1,0,1,...). For items that have negative deltas for sure, "-1" should be assigned; for items that have positive deltas, "1" should be assigned; for items whose deltas may be either positive or negative (e.g., intermediate items), "0" should assigned. We recommend at least two positive and two negative items per trait for better estimation.
#' @param trait The number of latent traits.
#' @param ind A row vector mapping each item to each trait. For example, c(1, 1, 1, 2, 2, 2) means that the first 3 items belong to trait 1 and the last 3 items belong to trait 2.
#' @param option The number of response options.
#' @param model Models fitted. They can be "UM8", "UM7", and "UM4". The default is UM8, which is the GGUM model. UM4 is UM8 with alpha = 1, called partial credit unfolding model. UM7 is UM8 with equal taus across items, called generalized rating scale unfolding model.
#' @param covariate An p*c person covariate matrix where p equals sample size and c equals the number of covariates. The default is NULL, meaning no person covariate.
#' @param iter The number of iterations. The default value is 1000. See documentation for rstan for more details.
#' @param chains The number of chains. The default value is 3. See documentation for rstan for more details.
#' @param warmup The number of warmups to discard. The default value is 0.5*iterations. See documentation for rstan for more details.
#' @param adapt_delta Target average proposal acceptance probability during Stan's adaptation period. The default value is 0.90. See documentation for rstan for more details.
#' @param max_treedepth Cap on the depth of the trees evaluated during each iteration. The default value is 15. See documentation for rstan for more details.
#' @param init Initial values for estimated parameters. The default is random initial values. See documentation for rstan for more details.
#' @param thin Thinning. The default value is 1. See documentation for rstan for more details.
#' @param cores The number of computer cores used for parallel computing. The default value is 2.
#' @param ma Mean of the prior distribution for alpha, which follows a lognormal distribution. The default value is 0.
#' @param va Standard deviation of the prior distribution for alpha. The default value is 0.5.
#' @param mdne Mean of the prior distribution for negative deltas, which follows a normal distribution. The default value is -1.
#' @param mdnu Mean of the prior distribution for neutral deltas, which follows a normal distribution. The default value is 0.
#' @param mdpo Mean of the prior distribution for positive deltas, which follows a normal distribution. The default value is 1.
#' @param vd Standard deviation of the prior distribution for deltas. The default value is 1.
#' @param mt Means of the prior distributions for taus, which follows a normal distribution. The default values are seq(-3, 0, 3/(options-1)). The last one has to be 0. For items with only 2 options, we recommend to use (-2, 0) as means of priors.
#' @param vt Standard deviation of the prior distribution for taus. The default value is 2.
#' @return Result object that stores information including the (1) stanfit object, (2) estimated item parameters, (3) estimated person parameters, (4) correlations among traits, (5) regression coefficients linking person covariates to each trait, (6) response data (excluding respondents who endorse a single option across all items), and (7) the input row vector mapping each item to each trait. Note that when covariates are included, output (4) represents residual correlations among the traits after controlling for the covariates.
#' @examples
#' \dontrun{
#' mod <- bmggum(GGUM.Data, delindex, trait, ind, covariate, option)}
#' @export
bmggum <- function(GGUM.Data, delindex, trait, ind, option, model="UM8", covariate=NULL, iter=1000, chains=3,
                   warmup=floor(iter/2), adapt_delta=0.90, max_treedepth=15, init="random", thin=1, cores=2,
                   ma=0, va=0.5, mdne=-1, mdnu=0, mdpo=1, vd=1, mt=seq(-3,0,3/(option-1)), vt=2){

dimension <- NULL
dimension <- ind
if (option==2){
    print(paste("For items with only 2 options, we recommend to manually input mt=c(-2,0) in the function as the means of priors for taus to replace the default (-3,0)"))
  }

  if (model=="UM8"){

    if (is.null(covariate)){

      #delete participants that endorse a single option across all items
      Valid <- matrix(NA,nrow=nrow(GGUM.Data),ncol=1)
      for (i in 1:nrow(GGUM.Data)){
        if (colSums(is.na(as.matrix(GGUM.Data[i,])))==(ncol(delindex))-1){
          Valid[i]=1
        } else{
          Valid[i]<-stats::var(GGUM.Data[i,], na.rm=T)
        }
      }

      for (i in 1:nrow(GGUM.Data)){
        if (Valid[i]==0){
          print(paste("Case", i, "was deleted because they endorse the same response option across all items"))
        }
      }
      GGUM.Data<-GGUM.Data[which(Valid!=0),]

      #sample size
      N1 <- nrow(GGUM.Data)

      ###### Data preparation
      #add delta and ind to response data
      GGUM.Data <- rbind(GGUM.Data, delindex)
      GGUM.Data <- rbind(GGUM.Data, ind)
      #reorder response data based on delta
      GGUM.Data <- GGUM.Data[, order(GGUM.Data[N1+2,],decreasing=FALSE)]
      #reordered dimension index
      ind <- GGUM.Data[N1+3,]
      #reordered delindex
      delindex <- GGUM.Data[(N1+1):(N1+2),]
      #reordered response data
      GGUM.Data <- GGUM.Data[-((N1+1):(N1+3)),]
      #### indicate which dimension does each item belong to while taking into consideration of missing data
      #0=missing
      #1=non-missing
      Missing <- matrix(NA,nrow=(ncol(delindex))*N1,ncol=1)
      MissPattern<-data.frame(Missing=((as.numeric(is.na(t(GGUM.Data)))*-1)+1),ID=seq(1,((ncol(delindex))*N1),1))
      NoMiss<-subset(MissPattern,Missing==1)
      #table(MissPattern$Missing)
      ind<-rep(ind,N1)[NoMiss$ID]

      #####long format response data
      Data<-suppressWarnings(edstan::irt_data(response_matrix =GGUM.Data))
      I<-dim(GGUM.Data)[1]
      J<-dim(GGUM.Data)[2]
      K=max(Data$y)

      #####number of ne, nu, po delta items
      ne=0
      nu=0
      po=0
      for (i in 1:(ncol(delindex))){
        if (delindex[2,i]<0){
          ne=ne+1
        } else if (delindex[2,i]==0){
          nu=nu+1
        } else if (delindex[2,i]>0){
          po=po+1
        }
      }

      data_list<-list(
        I=Data$I,
        J=Data$J,
        N=Data$N,
        II=Data$ii,
        JJ=Data$jj,
        K=max(Data$y),
        I_ne=ne,
        I_nu=nu,
        I_po=po,
        I_NN=ne+nu,
        N_mis=(Data$I*Data$J-Data$N),
        y=Data$y,
        M=2*K,
        trait=trait,
        ind=ind,
        theta_mu=as.array(c(rep(0,trait))),
        ma=ma,
        va=va,
        mdne=mdne,
        mdnu=mdnu,
        mdpo=mdpo,
        vd=vd,
        mt=mt,
        vt=vt)

      #######################################################
      #######################################################
      #                                                     #
      #          Bayesian estimation of GGUM                #
      #                                                     #
      #######################################################
      #######################################################
      BZ_bmggum<-'

functions {

real GGUM(int y, real theta, real alpha, real delta, vector tau, int K, int M) {

     vector[K] prob;

     vector[K] nominator;

     for (a in 1:K) {

     nominator[a]=exp(alpha*((a-1)*(theta-delta)-sum(tau[1:a])))+exp(alpha*((M-a)*(theta-delta)-sum(tau[1:a])));

     }

     prob=nominator/sum(nominator);

     return categorical_lpmf(y|prob);
}
}

data {
  int<lower=1> K;
  int<lower=2> M;  // M=2K
  int<lower=1> I;
  int<lower=1> J;
  int<lower=1> N;
  int<lower=0> I_ne;
  int<lower=0> I_nu;
  int<lower=0> I_po;
  int<lower=1> I_NN;
  int<lower=0> N_mis;
  int<lower=1, upper=I> II[N];
  int<lower=1, upper=J> JJ[N];
  int<lower=0, upper=6> y[N];

  // user-defined priors
  real ma;
  real va;
  real mdne;
  real mdnu;
  real mdpo;
  real vd;
  vector<lower=-5,upper=0>[K] mt;
  real vt;

  int<lower=1,upper=I>trait;
  int<lower=1,upper=I>ind[N];
  vector[trait] theta_mu;
}

parameters {

  vector<lower=0,upper=4>[I] alpha;
  vector<lower=0,upper=5>[I_po]   delta_po;
  vector<lower=-5,upper=0>[I_ne]  delta_ng;
  vector<lower=-5,upper=5>[I_nu]  delta_nu;
  vector<lower=-5,upper=0>[K-1]   tau_raw[I];

  vector[trait] theta[J];
  cholesky_factor_corr[trait] L_Omega;
}

transformed parameters{

  vector[I_NN] delta1;
  vector[I] delta;
  vector[K] tau[I];

  delta1=append_row(delta_ng,delta_nu);
  delta=append_row(delta1,delta_po);

  for (i in 1:I){

    tau[i,1]=0;
  }

  for (d in 2:K){

    tau[,d]=tau_raw[,(d-1)];
}
}


model {

    alpha ~    lognormal(ma,va);
    delta_ng ~ normal(mdne,vd);
    delta_nu ~ normal(mdnu,vd);
    delta_po ~ normal(mdpo,vd);
    L_Omega  ~ lkj_corr_cholesky(1);

    for (t in 1:(K-1)){
    tau_raw[,t] ~ normal(mt[t],vt);
    }


	  theta~ multi_normal_cholesky(theta_mu,L_Omega);

for (n in 1:N){

  target += GGUM(y[n],theta[JJ[n],ind[n]],alpha[II[n]],delta[II[n]],tau[II[n],],K,M);

}
}

generated quantities{
matrix[trait,trait] Cor;
vector[N] log_lik;

Cor=multiply_lower_tri_self_transpose(L_Omega);

  for (n in 1:N) {
    log_lik[n] = GGUM(y[n],theta[JJ[n],ind[n]],alpha[II[n]],delta[II[n]],tau[II[n],],K,M);
  }
}
'

##################################################
#       Input response data estimation        #
##################################################

rstan::rstan_options(auto_write = TRUE,javascript = FALSE)

BZ_bmggum<-rstan::stan(model_code=BZ_bmggum,data=data_list,
                       iter=iter, chains=chains,cores=cores, warmup=warmup,
                       init=init, thin=thin,
                       control=list(adapt_delta=adapt_delta,max_treedepth=max_treedepth))

#####Extract some parameters
THETA<-rstan::summary(BZ_bmggum, pars = c("theta"), probs = c(0.025,0.5,0.975))$summary
Alpha_ES<-rstan::summary(BZ_bmggum, pars = c("alpha"), probs = c(0.025,0.5,0.975))$summary
Delta_ES<-rstan::summary(BZ_bmggum, pars = c("delta"), probs = c(0.025, 0.5,0.975))$summary
Tau_ES<-rstan::summary(BZ_bmggum, pars = c("tau_raw"), probs = c(0.025, 0.5,0.975))$summary
Cor_ES<-rstan::summary(BZ_bmggum, pars = c("Cor"), probs = c(0.025, 0.5,0.975))$summary

#####reorder item parameters and data to the original input data order
delindex <- t(delindex)
estimates <- cbind(Alpha_ES, Delta_ES, delindex)
estimates <- estimates[order(estimates[,8+8+1],decreasing=FALSE),]
Alpha_ES <- estimates[,1:8]
Delta_ES <- estimates[,9:16]
rownames(Alpha_ES) <- c()
rownames(Delta_ES) <- c()
delindex1 <- delindex[rep(seq_len(nrow(delindex)), each = (option-1)), ]
estimates1 <- cbind(Tau_ES, delindex1)
estimates1 <- estimates1[order(estimates1[,8+1],decreasing=FALSE),]
Tau_ES <- estimates1[,1:8]
delindex <- t(delindex)
GGUM.Data <- rbind(GGUM.Data, delindex)
GGUM.Data <- GGUM.Data[, order(GGUM.Data[N1+1,],decreasing=FALSE)]
GGUM.Data <- GGUM.Data[-((N1+1):(N1+2)),]

#####save estimated parameters to an R object
MGGUM.summary<-list(Theta.est=THETA,
                    Alpha.est=Alpha_ES,
                    Delta.est=Delta_ES,
                    Tau.est=Tau_ES,
                    Cor.est=Cor_ES,
                    Data=GGUM.Data,
                    Fit=BZ_bmggum,
                    Dimension=dimension)
    } else{

  #delete participants that endorse a single option across all items
      Valid <- matrix(NA,nrow=nrow(GGUM.Data),ncol=1)
      for (i in 1:nrow(GGUM.Data)){
        if (colSums(is.na(as.matrix(GGUM.Data[i,])))==(ncol(delindex))-1){
          Valid[i]=1
        } else{
          Valid[i]<-stats::var(GGUM.Data[i,], na.rm=T)
        }
      }

  for (i in 1:nrow(GGUM.Data)){
    if (Valid[i]==0){
      print(paste("Case", i, "was deleted because they endorse the same response option across all items"))
    }
  }
  GGUM.Data<-GGUM.Data[which(Valid!=0),]
  covariate <- as.matrix(covariate)
  covariate<-covariate[which(Valid!=0),]
  covariate <- as.matrix(covariate)

  #number of covariates
  P <- ncol(covariate)
  #sample size
  N1 <- nrow(GGUM.Data)

  ###### Data preparation
  #add delta and ind to response data
  GGUM.Data <- rbind(GGUM.Data, delindex)
  GGUM.Data <- rbind(GGUM.Data, ind)
  #reorder response data based on delta
  GGUM.Data <- GGUM.Data[, order(GGUM.Data[N1+2,],decreasing=FALSE)]
  #reordered dimension index
  ind <- GGUM.Data[N1+3,]
  #reordered delindex
  delindex <- GGUM.Data[(N1+1):(N1+2),]
  #reordered response data
  GGUM.Data <- GGUM.Data[-((N1+1):(N1+3)),]
  #### indicate which dimension does each item belong to while taking into consideration of missing data
  #0=missing
  #1=non-missing
  Missing <- matrix(NA,nrow=(ncol(delindex))*N1,ncol=1)
  MissPattern<-data.frame(Missing=((as.numeric(is.na(t(GGUM.Data)))*-1)+1),ID=seq(1,((ncol(delindex))*N1),1))
  NoMiss<-subset(MissPattern,Missing==1)
  #table(MissPattern$Missing)
  ind<-rep(ind,N1)[NoMiss$ID]

  #####long format response data
  Data<-suppressWarnings(edstan::irt_data(response_matrix =GGUM.Data))
  I<-dim(GGUM.Data)[1]
  J<-dim(GGUM.Data)[2]
  K=max(Data$y)

  #####number of ne, nu, po delta items
  ne=0
  nu=0
  po=0
  for (i in 1:(ncol(delindex))){
    if (delindex[2,i]<0){
      ne=ne+1
    } else if (delindex[2,i]==0){
      nu=nu+1
    } else if (delindex[2,i]>0){
      po=po+1
    }
  }

  #####long format person covariates
  PC1<-covariate[rep(seq_len(nrow(covariate)), each=trait), ]
  PC1 <- as.matrix(PC1)

  data_list<-list(
    I=Data$I,
    J=Data$J,
    N=Data$N,
    II=Data$ii,
    JJ=Data$jj,
    K=max(Data$y),
    I_ne=ne,
    I_nu=nu,
    I_po=po,
    I_NN=ne+nu,
    N_mis=(Data$I*Data$J-Data$N),
    y=Data$y,
    M=2*K,
    trait=trait,
    ind=ind,
    theta_mu=as.array(c(rep(0,trait))),
    P=P,
    PC1=PC1,
    ma=ma,
    va=va,
    mdne=mdne,
    mdnu=mdnu,
    mdpo=mdpo,
    vd=vd,
    mt=mt,
    vt=vt)

  #######################################################
  #######################################################
  #                                                     #
  #          Bayesian estimation of GGUM                #
  #                                                     #
  #######################################################
  #######################################################
  BZ_bmggum<-'

functions {

real GGUM(int y, real theta, real alpha, real delta, vector tau, int K, int M) {

     vector[K] prob;

     vector[K] nominator;

     for (a in 1:K) {

     nominator[a]=exp(alpha*((a-1)*(theta-delta)-sum(tau[1:a])))+exp(alpha*((M-a)*(theta-delta)-sum(tau[1:a])));

     }

     prob=nominator/sum(nominator);

     return categorical_lpmf(y|prob);
}
}

data {
  int<lower=1> K;
  int<lower=2> M;  // M=2K
  int<lower=1> I;
  int<lower=1> J;
  int<lower=1> N;
  int<lower=1> I_ne;
  int<lower=0> I_nu;
  int<lower=1> I_po;
  int<lower=1> I_NN;
  int<lower=0> N_mis;
  int<lower=1, upper=I> II[N];
  int<lower=1, upper=J> JJ[N];
  int<lower=0, upper=6> y[N];

  // user-defined priors
  real ma;
  real va;
  real mdne;
  real mdnu;
  real mdpo;
  real vd;
  vector<lower=-5,upper=0>[K] mt;
  real vt;

  int<lower=1,upper=I>trait;
  int<lower=1,upper=I>ind[N];
  vector[trait] theta_mu;
  int<lower=1> P;
  matrix[trait*J,P] PC1;
}

parameters {

  vector<lower=0,upper=4>[I] alpha;
  vector<lower=0,upper=5>[I_po]   delta_po;
  vector<lower=-5,upper=0>[I_ne]  delta_ng;
  vector<lower=-5,upper=5>[I_nu]  delta_nu;
  vector<lower=-5,upper=0>[K-1]   tau_raw[I];

  vector[trait] theta[J];
  cholesky_factor_corr[trait] L_Omega;
  matrix[P,trait] lambda;    // person covariate regression coefficient
}

transformed parameters{

  vector[I_NN] delta1;
  vector[I] delta;
  vector[K] tau[I];

  delta1=append_row(delta_ng,delta_nu);
  delta=append_row(delta1,delta_po);

  for (i in 1:I){

    tau[i,1]=0;
  }

  for (d in 2:K){

    tau[,d]=tau_raw[,(d-1)];
}
}


model {

    matrix[trait, trait] center_mu;  //PC1*lambda matrix
    row_vector[trait] center_mu0;    //first row of PC1*lambda matrix
    vector[trait] center_mu1;        //change center_mu0 row to column

    alpha ~    lognormal(ma,va);
    delta_ng ~ normal(mdne,vd);
    delta_nu ~ normal(mdnu,vd);
    delta_po ~ normal(mdpo,vd);
    L_Omega  ~ lkj_corr_cholesky(1);

    for (t in 1:(K-1)){
    tau_raw[,t] ~ normal(mt[t],vt);
    }


	//theta~ multi_normal_cholesky(theta_mu,L_Omega);
    for (j in 1:J){
      center_mu=PC1[(trait*j-(trait-1)):trait*j,]*lambda;
      center_mu0=center_mu[1,];
      center_mu1=to_vector(center_mu0);
      theta[j]~ multi_normal_cholesky(theta_mu + center_mu1, L_Omega);
    }
     to_vector(lambda) ~ normal (0, 1);

for (n in 1:N){

  target += GGUM(y[n],theta[JJ[n],ind[n]],alpha[II[n]],delta[II[n]],tau[II[n],],K,M);

}
}

generated quantities{
matrix[trait,trait] Cor;
vector[N] log_lik;

Cor=multiply_lower_tri_self_transpose(L_Omega);

  for (n in 1:N) {
    log_lik[n] = GGUM(y[n],theta[JJ[n],ind[n]],alpha[II[n]],delta[II[n]],tau[II[n],],K,M);
  }
}
'

##################################################
#       Input response data estimation        #
##################################################

rstan::rstan_options(auto_write = TRUE,javascript = FALSE)

BZ_bmggum<-rstan::stan(model_code=BZ_bmggum,data=data_list,
                       iter=iter, chains=chains,cores=cores, warmup=warmup,
                       init=init, thin=thin,
                       control=list(adapt_delta=adapt_delta,max_treedepth=max_treedepth))

#####Extract some parameters
THETA<-rstan::summary(BZ_bmggum, pars = c("theta"), probs = c(0.025,0.5,0.975))$summary
Alpha_ES<-rstan::summary(BZ_bmggum, pars = c("alpha"), probs = c(0.025,0.5,0.975))$summary
Delta_ES<-rstan::summary(BZ_bmggum, pars = c("delta"), probs = c(0.025, 0.5,0.975))$summary
Tau_ES<-rstan::summary(BZ_bmggum, pars = c("tau_raw"), probs = c(0.025, 0.5,0.975))$summary
Cor_ES<-rstan::summary(BZ_bmggum, pars = c("Cor"), probs = c(0.025, 0.5,0.975))$summary
Lamda_ES<-rstan::summary(BZ_bmggum, pars = c("lambda"), probs = c(0.025, 0.5,0.975))$summary

#####reorder item parameters and data to the original input data order
delindex <- t(delindex)
estimates <- cbind(Alpha_ES, Delta_ES, delindex)
estimates <- estimates[order(estimates[,8+8+1],decreasing=FALSE),]
Alpha_ES <- estimates[,1:8]
Delta_ES <- estimates[,9:16]
rownames(Alpha_ES) <- c()
rownames(Delta_ES) <- c()
delindex1 <- delindex[rep(seq_len(nrow(delindex)), each = (option-1)), ]
estimates1 <- cbind(Tau_ES, delindex1)
estimates1 <- estimates1[order(estimates1[,8+1],decreasing=FALSE),]
Tau_ES <- estimates1[,1:8]
delindex <- t(delindex)
GGUM.Data <- rbind(GGUM.Data, delindex)
GGUM.Data <- GGUM.Data[, order(GGUM.Data[N1+1,],decreasing=FALSE)]
GGUM.Data <- GGUM.Data[-((N1+1):(N1+2)),]

#####save estimated parameters to an R object
MGGUM.summary<-list(Theta.est=THETA,
                    Alpha.est=Alpha_ES,
                    Delta.est=Delta_ES,
                    Tau.est=Tau_ES,
                    Cor.est=Cor_ES,
                    Lamda.est=Lamda_ES,
                    Data=GGUM.Data,
                    Fit=BZ_bmggum,
                    Dimension=dimension)
}
  } else if (model=="UM4"){
    if (is.null(covariate)){

      #delete participants that endorse a single option across all items
      Valid <- matrix(NA,nrow=nrow(GGUM.Data),ncol=1)
      for (i in 1:nrow(GGUM.Data)){
        if (colSums(is.na(as.matrix(GGUM.Data[i,])))==(ncol(delindex))-1){
          Valid[i]=1
        } else{
          Valid[i]<-stats::var(GGUM.Data[i,], na.rm=T)
        }
      }

      for (i in 1:nrow(GGUM.Data)){
        if (Valid[i]==0){
          print(paste("Case", i, "was deleted because they endorse the same response option across all items"))
        }
      }
      GGUM.Data<-GGUM.Data[which(Valid!=0),]

      #sample size
      N1 <- nrow(GGUM.Data)

      ###### Data preparation
      #add delta and ind to response data
      GGUM.Data <- rbind(GGUM.Data, delindex)
      GGUM.Data <- rbind(GGUM.Data, ind)
      #reorder response data based on delta
      GGUM.Data <- GGUM.Data[, order(GGUM.Data[N1+2,],decreasing=FALSE)]
      #reordered dimension index
      ind <- GGUM.Data[N1+3,]
      #reordered delindex
      delindex <- GGUM.Data[(N1+1):(N1+2),]
      #reordered response data
      GGUM.Data <- GGUM.Data[-((N1+1):(N1+3)),]
      #### indicate which dimension does each item belong to while taking into consideration of missing data
      #0=missing
      #1=non-missing
      Missing <- matrix(NA,nrow=(ncol(delindex))*N1,ncol=1)
      MissPattern<-data.frame(Missing=((as.numeric(is.na(t(GGUM.Data)))*-1)+1),ID=seq(1,((ncol(delindex))*N1),1))
      NoMiss<-subset(MissPattern,Missing==1)
      #table(MissPattern$Missing)
      ind<-rep(ind,N1)[NoMiss$ID]

      #####long format response data
      Data<-suppressWarnings(edstan::irt_data(response_matrix =GGUM.Data))
      I<-dim(GGUM.Data)[1]
      J<-dim(GGUM.Data)[2]
      K=max(Data$y)

      #####number of ne, nu, po delta items
      ne=0
      nu=0
      po=0
      for (i in 1:(ncol(delindex))){
        if (delindex[2,i]<0){
          ne=ne+1
        } else if (delindex[2,i]==0){
          nu=nu+1
        } else if (delindex[2,i]>0){
          po=po+1
        }
      }

      data_list<-list(
        I=Data$I,
        J=Data$J,
        N=Data$N,
        II=Data$ii,
        JJ=Data$jj,
        K=max(Data$y),
        I_ne=ne,
        I_nu=nu,
        I_po=po,
        I_NN=ne+nu,
        N_mis=(Data$I*Data$J-Data$N),
        y=Data$y,
        M=2*K,
        trait=trait,
        ind=ind,
        theta_mu=as.array(c(rep(0,trait))),
        ma=ma,
        va=va,
        mdne=mdne,
        mdnu=mdnu,
        mdpo=mdpo,
        vd=vd,
        mt=mt,
        vt=vt)

      #######################################################
      #######################################################
      #                                                     #
      #          Bayesian estimation of GGUM                #
      #                                                     #
      #######################################################
      #######################################################
      BZ_bmggum<-'

functions {

real GGUM(int y, real theta, real delta, vector tau, int K, int M) {

     vector[K] prob;

     vector[K] nominator;

     for (a in 1:K) {

     nominator[a]=exp(1*((a-1)*(theta-delta)-sum(tau[1:a])))+exp(1*((M-a)*(theta-delta)-sum(tau[1:a])));

     }

     prob=nominator/sum(nominator);

     return categorical_lpmf(y|prob);
}
}

data {
  int<lower=1> K;
  int<lower=2> M;  // M=2K
  int<lower=1> I;
  int<lower=1> J;
  int<lower=1> N;
  int<lower=1> I_ne;
  int<lower=0> I_nu;
  int<lower=1> I_po;
  int<lower=1> I_NN;
  int<lower=0> N_mis;
  int<lower=1, upper=I> II[N];
  int<lower=1, upper=J> JJ[N];
  int<lower=0, upper=6> y[N];

  // user-defined priors
  real ma;
  real va;
  real mdne;
  real mdnu;
  real mdpo;
  real vd;
  vector<lower=-5,upper=0>[K] mt;
  real vt;

  int<lower=1,upper=I>trait;
  int<lower=1,upper=I>ind[N];
  vector[trait] theta_mu;
}

parameters {

  //vector<lower=0,upper=4>[I] alpha;
  vector<lower=0,upper=5>[I_po]   delta_po;
  vector<lower=-5,upper=0>[I_ne]  delta_ng;
  vector<lower=-5,upper=5>[I_nu]  delta_nu;
  vector<lower=-5,upper=0>[K-1]   tau_raw[I];

  vector[trait] theta[J];
  cholesky_factor_corr[trait] L_Omega;
}

transformed parameters{

  vector[I_NN] delta1;
  vector[I] delta;
  vector[K] tau[I];

  delta1=append_row(delta_ng,delta_nu);
  delta=append_row(delta1,delta_po);

  for (i in 1:I){

    tau[i,1]=0;
  }

  for (d in 2:K){

    tau[,d]=tau_raw[,(d-1)];
}
}


model {

    //alpha ~    lognormal(ma,va);
    delta_ng ~ normal(mdne,vd);
    delta_nu ~ normal(mdnu,vd);
    delta_po ~ normal(mdpo,vd);
    L_Omega  ~ lkj_corr_cholesky(1);

    for (t in 1:(K-1)){
    tau_raw[,t] ~ normal(mt[t],vt);
    }


	  theta~ multi_normal_cholesky(theta_mu,L_Omega);

for (n in 1:N){

  target += GGUM(y[n],theta[JJ[n],ind[n]],delta[II[n]],tau[II[n],],K,M);

}
}

generated quantities{
matrix[trait,trait] Cor;
vector[N] log_lik;

Cor=multiply_lower_tri_self_transpose(L_Omega);

  for (n in 1:N) {
    log_lik[n] = GGUM(y[n],theta[JJ[n],ind[n]],delta[II[n]],tau[II[n],],K,M);
  }
}
'

##################################################
#       Input response data estimation        #
##################################################

rstan::rstan_options(auto_write = TRUE,javascript = FALSE)

BZ_bmggum<-rstan::stan(model_code=BZ_bmggum,data=data_list,
                       iter=iter, chains=chains,cores=cores, warmup=warmup,
                       init=init, thin=thin,
                       control=list(adapt_delta=adapt_delta,max_treedepth=max_treedepth))

#####Extract some parameters
THETA<-rstan::summary(BZ_bmggum, pars = c("theta"), probs = c(0.025,0.5,0.975))$summary
#Alpha_ES<-rstan::summary(BZ_bmggum, pars = c("alpha"), probs = c(0.025,0.5,0.975))$summary
Delta_ES<-rstan::summary(BZ_bmggum, pars = c("delta"), probs = c(0.025, 0.5,0.975))$summary
Tau_ES<-rstan::summary(BZ_bmggum, pars = c("tau_raw"), probs = c(0.025, 0.5,0.975))$summary
Cor_ES<-rstan::summary(BZ_bmggum, pars = c("Cor"), probs = c(0.025, 0.5,0.975))$summary
Alpha_ES<-rep(1, ncol(delindex))
Alpha_ES <- as.matrix(Alpha_ES)

#####reorder item parameters and data to the original input data order
delindex <- t(delindex)
estimates <- cbind(Delta_ES, delindex)
estimates <- estimates[order(estimates[,8+1],decreasing=FALSE),]
#Alpha_ES <- estimates[,1:8]
Delta_ES <- estimates[,1:8]
#rownames(Delta_ES) <- c()
delindex1 <- delindex[rep(seq_len(nrow(delindex)), each = (option-1)), ]
estimates1 <- cbind(Tau_ES, delindex1)
estimates1 <- estimates1[order(estimates1[,8+1],decreasing=FALSE),]
Tau_ES <- estimates1[,1:8]
delindex <- t(delindex)
GGUM.Data <- rbind(GGUM.Data, delindex)
GGUM.Data <- GGUM.Data[, order(GGUM.Data[N1+1,],decreasing=FALSE)]
GGUM.Data <- GGUM.Data[-((N1+1):(N1+2)),]

#####save estimated parameters to an R object
MGGUM.summary<-list(Theta.est=THETA,
                    Alpha.est=Alpha_ES,
                    Delta.est=Delta_ES,
                    Tau.est=Tau_ES,
                    Cor.est=Cor_ES,
                    Data=GGUM.Data,
                    Fit=BZ_bmggum,
                    Dimension=dimension)
    } else{

  #delete participants that endorse a single option across all items
      Valid <- matrix(NA,nrow=nrow(GGUM.Data),ncol=1)
      for (i in 1:nrow(GGUM.Data)){
        if (colSums(is.na(as.matrix(GGUM.Data[i,])))==(ncol(delindex))-1){
          Valid[i]=1
        } else{
          Valid[i]<-stats::var(GGUM.Data[i,], na.rm=T)
        }
      }

  for (i in 1:nrow(GGUM.Data)){
    if (Valid[i]==0){
      print(paste("Case", i, "was deleted because they endorse the same response option across all items"))
    }
  }
  GGUM.Data<-GGUM.Data[which(Valid!=0),]
  covariate <- as.matrix(covariate)
  covariate<-covariate[which(Valid!=0),]
  covariate <- as.matrix(covariate)

  #number of covariates
  P <- ncol(covariate)
  #sample size
  N1 <- nrow(GGUM.Data)

  ###### Data preparation
  #add delta and ind to response data
  GGUM.Data <- rbind(GGUM.Data, delindex)
  GGUM.Data <- rbind(GGUM.Data, ind)
  #reorder response data based on delta
  GGUM.Data <- GGUM.Data[, order(GGUM.Data[N1+2,],decreasing=FALSE)]
  #reordered dimension index
  ind <- GGUM.Data[N1+3,]
  #reordered delindex
  delindex <- GGUM.Data[(N1+1):(N1+2),]
  #reordered response data
  GGUM.Data <- GGUM.Data[-((N1+1):(N1+3)),]
  #### indicate which dimension does each item belong to while taking into consideration of missing data
  #0=missing
  #1=non-missing
  Missing <- matrix(NA,nrow=(ncol(delindex))*N1,ncol=1)
  MissPattern<-data.frame(Missing=((as.numeric(is.na(t(GGUM.Data)))*-1)+1),ID=seq(1,((ncol(delindex))*N1),1))
  NoMiss<-subset(MissPattern,Missing==1)
  #table(MissPattern$Missing)
  ind<-rep(ind,N1)[NoMiss$ID]

  #####long format response data
  Data<-suppressWarnings(edstan::irt_data(response_matrix =GGUM.Data))
  I<-dim(GGUM.Data)[1]
  J<-dim(GGUM.Data)[2]
  K=max(Data$y)

  #####number of ne, nu, po delta items
  ne=0
  nu=0
  po=0
  for (i in 1:(ncol(delindex))){
    if (delindex[2,i]<0){
      ne=ne+1
    } else if (delindex[2,i]==0){
      nu=nu+1
    } else if (delindex[2,i]>0){
      po=po+1
    }
  }

  #####long format person covariates
  PC1<-covariate[rep(seq_len(nrow(covariate)), each=trait), ]
  PC1 <- as.matrix(PC1)

  data_list<-list(
    I=Data$I,
    J=Data$J,
    N=Data$N,
    II=Data$ii,
    JJ=Data$jj,
    K=max(Data$y),
    I_ne=ne,
    I_nu=nu,
    I_po=po,
    I_NN=ne+nu,
    N_mis=(Data$I*Data$J-Data$N),
    y=Data$y,
    M=2*K,
    trait=trait,
    ind=ind,
    theta_mu=as.array(c(rep(0,trait))),
    P=P,
    PC1=PC1,
    ma=ma,
    va=va,
    mdne=mdne,
    mdnu=mdnu,
    mdpo=mdpo,
    vd=vd,
    mt=mt,
    vt=vt)

  #######################################################
  #######################################################
  #                                                     #
  #          Bayesian estimation of GGUM                #
  #                                                     #
  #######################################################
  #######################################################
  BZ_bmggum<-'

functions {

real GGUM(int y, real theta, real delta, vector tau, int K, int M) {

     vector[K] prob;

     vector[K] nominator;

     for (a in 1:K) {

     nominator[a]=exp(1*((a-1)*(theta-delta)-sum(tau[1:a])))+exp(1*((M-a)*(theta-delta)-sum(tau[1:a])));

     }

     prob=nominator/sum(nominator);

     return categorical_lpmf(y|prob);
}
}

data {
  int<lower=1> K;
  int<lower=2> M;  // M=2K
  int<lower=1> I;
  int<lower=1> J;
  int<lower=1> N;
  int<lower=1> I_ne;
  int<lower=0> I_nu;
  int<lower=1> I_po;
  int<lower=1> I_NN;
  int<lower=0> N_mis;
  int<lower=1, upper=I> II[N];
  int<lower=1, upper=J> JJ[N];
  int<lower=0, upper=6> y[N];

  // user-defined priors
  real ma;
  real va;
  real mdne;
  real mdnu;
  real mdpo;
  real vd;
  vector<lower=-5,upper=0>[K] mt;
  real vt;

  int<lower=1,upper=I>trait;
  int<lower=1,upper=I>ind[N];
  vector[trait] theta_mu;
  int<lower=1> P;
  matrix[trait*J,P] PC1;
}

parameters {

  //vector<lower=0,upper=4>[I] alpha;
  vector<lower=0,upper=5>[I_po]   delta_po;
  vector<lower=-5,upper=0>[I_ne]  delta_ng;
  vector<lower=-5,upper=5>[I_nu]  delta_nu;
  vector<lower=-5,upper=0>[K-1]   tau_raw[I];

  vector[trait] theta[J];
  cholesky_factor_corr[trait] L_Omega;
  matrix[P,trait] lambda;    // person covariate regression coefficient
}

transformed parameters{

  vector[I_NN] delta1;
  vector[I] delta;
  vector[K] tau[I];

  delta1=append_row(delta_ng,delta_nu);
  delta=append_row(delta1,delta_po);

  for (i in 1:I){

    tau[i,1]=0;
  }

  for (d in 2:K){

    tau[,d]=tau_raw[,(d-1)];
}
}


model {

    matrix[trait, trait] center_mu;  //PC1*lambda matrix
    row_vector[trait] center_mu0;    //first row of PC1*lambda matrix
    vector[trait] center_mu1;        //change center_mu0 row to column

    //alpha ~    lognormal(ma,va);
    delta_ng ~ normal(mdne,vd);
    delta_nu ~ normal(mdnu,vd);
    delta_po ~ normal(mdpo,vd);
    L_Omega  ~ lkj_corr_cholesky(1);

    for (t in 1:(K-1)){
    tau_raw[,t] ~ normal(mt[t],vt);
    }


	//theta~ multi_normal_cholesky(theta_mu,L_Omega);
    for (j in 1:J){
      center_mu=PC1[(trait*j-(trait-1)):trait*j,]*lambda;
      center_mu0=center_mu[1,];
      center_mu1=to_vector(center_mu0);
      theta[j]~ multi_normal_cholesky(theta_mu + center_mu1, L_Omega);
    }
     to_vector(lambda) ~ normal (0, 1);

for (n in 1:N){

  target += GGUM(y[n],theta[JJ[n],ind[n]],delta[II[n]],tau[II[n],],K,M);

}
}

generated quantities{
matrix[trait,trait] Cor;
vector[N] log_lik;

Cor=multiply_lower_tri_self_transpose(L_Omega);

  for (n in 1:N) {
    log_lik[n] = GGUM(y[n],theta[JJ[n],ind[n]],delta[II[n]],tau[II[n],],K,M);
  }
}
'

##################################################
#       Input response data estimation        #
##################################################

rstan::rstan_options(auto_write = TRUE,javascript = FALSE)

BZ_bmggum<-rstan::stan(model_code=BZ_bmggum,data=data_list,
                       iter=iter, chains=chains,cores=cores, warmup=warmup,
                       init=init, thin=thin,
                       control=list(adapt_delta=adapt_delta,max_treedepth=max_treedepth))

#####Extract some parameters
THETA<-rstan::summary(BZ_bmggum, pars = c("theta"), probs = c(0.025,0.5,0.975))$summary
#Alpha_ES<-rstan::summary(BZ_bmggum, pars = c("alpha"), probs = c(0.025,0.5,0.975))$summary
Delta_ES<-rstan::summary(BZ_bmggum, pars = c("delta"), probs = c(0.025, 0.5,0.975))$summary
Tau_ES<-rstan::summary(BZ_bmggum, pars = c("tau_raw"), probs = c(0.025, 0.5,0.975))$summary
Cor_ES<-rstan::summary(BZ_bmggum, pars = c("Cor"), probs = c(0.025, 0.5,0.975))$summary
Lamda_ES<-rstan::summary(BZ_bmggum, pars = c("lambda"), probs = c(0.025, 0.5,0.975))$summary
Alpha_ES<-rep(1, ncol(delindex))
Alpha_ES <- as.matrix(Alpha_ES)

#####reorder item parameters and data to the original input data order
delindex <- t(delindex)
estimates <- cbind(Delta_ES, delindex)
estimates <- estimates[order(estimates[,8+1],decreasing=FALSE),]
#Alpha_ES <- estimates[,1:8]
Delta_ES <- estimates[,1:8]
#rownames(Delta_ES) <- c()
delindex1 <- delindex[rep(seq_len(nrow(delindex)), each = (option-1)), ]
estimates1 <- cbind(Tau_ES, delindex1)
estimates1 <- estimates1[order(estimates1[,8+1],decreasing=FALSE),]
Tau_ES <- estimates1[,1:8]
delindex <- t(delindex)
GGUM.Data <- rbind(GGUM.Data, delindex)
GGUM.Data <- GGUM.Data[, order(GGUM.Data[N1+1,],decreasing=FALSE)]
GGUM.Data <- GGUM.Data[-((N1+1):(N1+2)),]

#####save estimated parameters to an R object
MGGUM.summary<-list(Theta.est=THETA,
                    Alpha.est=Alpha_ES,
                    Delta.est=Delta_ES,
                    Tau.est=Tau_ES,
                    Cor.est=Cor_ES,
                    Lamda.est=Lamda_ES,
                    Data=GGUM.Data,
                    Fit=BZ_bmggum,
                    Dimension=dimension)
}
  } else if (model=="UM7"){
    if (is.null(covariate)){

      #delete participants that endorse a single option across all items
      Valid <- matrix(NA,nrow=nrow(GGUM.Data),ncol=1)
      for (i in 1:nrow(GGUM.Data)){
        if (colSums(is.na(as.matrix(GGUM.Data[i,])))==(ncol(delindex))-1){
          Valid[i]=1
        } else{
          Valid[i]<-stats::var(GGUM.Data[i,], na.rm=T)
        }
      }

      for (i in 1:nrow(GGUM.Data)){
        if (Valid[i]==0){
          print(paste("Case", i, "was deleted because they endorse the same response option across all items"))
        }
      }
      GGUM.Data<-GGUM.Data[which(Valid!=0),]

      #sample size
      N1 <- nrow(GGUM.Data)

      ###### Data preparation
      #add delta and ind to response data
      GGUM.Data <- rbind(GGUM.Data, delindex)
      GGUM.Data <- rbind(GGUM.Data, ind)
      #reorder response data based on delta
      GGUM.Data <- GGUM.Data[, order(GGUM.Data[N1+2,],decreasing=FALSE)]
      #reordered dimension index
      ind <- GGUM.Data[N1+3,]
      #reordered delindex
      delindex <- GGUM.Data[(N1+1):(N1+2),]
      #reordered response data
      GGUM.Data <- GGUM.Data[-((N1+1):(N1+3)),]
      #### indicate which dimension does each item belong to while taking into consideration of missing data
      #0=missing
      #1=non-missing
      Missing <- matrix(NA,nrow=(ncol(delindex))*N1,ncol=1)
      MissPattern<-data.frame(Missing=((as.numeric(is.na(t(GGUM.Data)))*-1)+1),ID=seq(1,((ncol(delindex))*N1),1))
      NoMiss<-subset(MissPattern,Missing==1)
      #table(MissPattern$Missing)
      ind<-rep(ind,N1)[NoMiss$ID]

      #####long format response data
      Data<-suppressWarnings(edstan::irt_data(response_matrix =GGUM.Data))
      I<-dim(GGUM.Data)[1]
      J<-dim(GGUM.Data)[2]
      K=max(Data$y)

      #####number of ne, nu, po delta items
      ne=0
      nu=0
      po=0
      for (i in 1:(ncol(delindex))){
        if (delindex[2,i]<0){
          ne=ne+1
        } else if (delindex[2,i]==0){
          nu=nu+1
        } else if (delindex[2,i]>0){
          po=po+1
        }
      }

      data_list<-list(
        I=Data$I,
        J=Data$J,
        N=Data$N,
        II=Data$ii,
        JJ=Data$jj,
        K=max(Data$y),
        I_ne=ne,
        I_nu=nu,
        I_po=po,
        I_NN=ne+nu,
        N_mis=(Data$I*Data$J-Data$N),
        y=Data$y,
        M=2*K,
        trait=trait,
        ind=ind,
        theta_mu=as.array(c(rep(0,trait))),
        ma=ma,
        va=va,
        mdne=mdne,
        mdnu=mdnu,
        mdpo=mdpo,
        vd=vd,
        mt=mt,
        vt=vt)

      #######################################################
      #######################################################
      #                                                     #
      #          Bayesian estimation of GGUM                #
      #                                                     #
      #######################################################
      #######################################################
      BZ_bmggum<-'

functions {

real GGUM(int y, real theta, real alpha, real delta, vector tau, int K, int M) {

     vector[K] prob;

     vector[K] nominator;

     for (a in 1:K) {

     nominator[a]=exp(alpha*((a-1)*(theta-delta)-sum(tau[1:a])))+exp(alpha*((M-a)*(theta-delta)-sum(tau[1:a])));

     }

     prob=nominator/sum(nominator);

     return categorical_lpmf(y|prob);
}
}

data {
  int<lower=1> K;
  int<lower=2> M;  // M=2K
  int<lower=1> I;
  int<lower=1> J;
  int<lower=1> N;
  int<lower=1> I_ne;
  int<lower=0> I_nu;
  int<lower=1> I_po;
  int<lower=1> I_NN;
  int<lower=0> N_mis;
  int<lower=1, upper=I> II[N];
  int<lower=1, upper=J> JJ[N];
  int<lower=0, upper=6> y[N];

  // user-defined priors
  real ma;
  real va;
  real mdne;
  real mdnu;
  real mdpo;
  real vd;
  vector<lower=-5,upper=0>[K] mt;
  real vt;

  int<lower=1,upper=I>trait;
  int<lower=1,upper=I>ind[N];
  vector[trait] theta_mu;
}

parameters {

  vector<lower=0,upper=4>[I] alpha;
  vector<lower=0,upper=5>[I_po]   delta_po;
  vector<lower=-5,upper=0>[I_ne]  delta_ng;
  vector<lower=-5,upper=5>[I_nu]  delta_nu;
  vector<lower=-5,upper=0>[K-1]   tau_raw;

  vector[trait] theta[J];
  cholesky_factor_corr[trait] L_Omega;
}

transformed parameters{

  vector[I_NN] delta1;
  vector[I] delta;
  vector[K] tau;

  delta1=append_row(delta_ng,delta_nu);
  delta=append_row(delta1,delta_po);

  tau[1]=0;
  for (d in 2:K){

    tau[d]=tau_raw[d-1];
}
}


model {

    alpha ~    lognormal(ma,va);
    delta_ng ~ normal(mdne,vd);
    delta_nu ~ normal(mdnu,vd);
    delta_po ~ normal(mdpo,vd);
    L_Omega  ~ lkj_corr_cholesky(1);

    for (t in 1:(K-1)){
    tau_raw[t] ~ normal(mt[t],vt);
    }


	  theta~ multi_normal_cholesky(theta_mu,L_Omega);

for (n in 1:N){

  target += GGUM(y[n],theta[JJ[n],ind[n]],alpha[II[n]],delta[II[n]],tau[],K,M);

}
}

generated quantities{
matrix[trait,trait] Cor;
vector[N] log_lik;

Cor=multiply_lower_tri_self_transpose(L_Omega);

  for (n in 1:N) {
    log_lik[n] = GGUM(y[n],theta[JJ[n],ind[n]],alpha[II[n]],delta[II[n]],tau[],K,M);
  }
}
'

##################################################
#       Input response data estimation        #
##################################################

rstan::rstan_options(auto_write = TRUE,javascript = FALSE)

BZ_bmggum<-rstan::stan(model_code=BZ_bmggum,data=data_list,
                       iter=iter, chains=chains,cores=cores, warmup=warmup,
                       init=init, thin=thin,
                       control=list(adapt_delta=adapt_delta,max_treedepth=max_treedepth))

#####Extract some parameters
THETA<-rstan::summary(BZ_bmggum, pars = c("theta"), probs = c(0.025,0.5,0.975))$summary
Alpha_ES<-rstan::summary(BZ_bmggum, pars = c("alpha"), probs = c(0.025,0.5,0.975))$summary
Delta_ES<-rstan::summary(BZ_bmggum, pars = c("delta"), probs = c(0.025, 0.5,0.975))$summary
Tau_ES<-rstan::summary(BZ_bmggum, pars = c("tau_raw"), probs = c(0.025, 0.5,0.975))$summary
Cor_ES<-rstan::summary(BZ_bmggum, pars = c("Cor"), probs = c(0.025, 0.5,0.975))$summary

#####reorder item parameters and data to the original input data order
delindex <- t(delindex)
estimates <- cbind(Alpha_ES, Delta_ES, delindex)
estimates <- estimates[order(estimates[,8+8+1],decreasing=FALSE),]
Alpha_ES <- estimates[,1:8]
Delta_ES <- estimates[,9:16]
rownames(Alpha_ES) <- c()
rownames(Delta_ES) <- c()
#change tau estimates to the same number of rows as alpha UM8
Tau_ES <- Tau_ES[rep(seq_len(nrow(Tau_ES)), nrow(delindex)), ]
# delindex1 <- delindex[rep(seq_len(nrow(delindex)), each = (option-1)), ]
# estimates1 <- cbind(Tau_ES, delindex1)
# estimates1 <- estimates1[order(estimates1[,8+1],decreasing=FALSE),]
# Tau_ES <- estimates1[,1:8]
delindex <- t(delindex)
GGUM.Data <- rbind(GGUM.Data, delindex)
GGUM.Data <- GGUM.Data[, order(GGUM.Data[N1+1,],decreasing=FALSE)]
GGUM.Data <- GGUM.Data[-((N1+1):(N1+2)),]

#####save estimated parameters to an R object
MGGUM.summary<-list(Theta.est=THETA,
                    Alpha.est=Alpha_ES,
                    Delta.est=Delta_ES,
                    Tau.est=Tau_ES,
                    Cor.est=Cor_ES,
                    Data=GGUM.Data,
                    Fit=BZ_bmggum,
                    Dimension=dimension)
    } else{

  #delete participants that endorse a single option across all items
      Valid <- matrix(NA,nrow=nrow(GGUM.Data),ncol=1)
      for (i in 1:nrow(GGUM.Data)){
        if (colSums(is.na(as.matrix(GGUM.Data[i,])))==(ncol(delindex))-1){
          Valid[i]=1
        } else{
          Valid[i]<-stats::var(GGUM.Data[i,], na.rm=T)
        }
      }

  for (i in 1:nrow(GGUM.Data)){
    if (Valid[i]==0){
      print(paste("Case", i, "was deleted because they endorse the same response option across all items"))
    }
  }
  GGUM.Data<-GGUM.Data[which(Valid!=0),]
  covariate <- as.matrix(covariate)
  covariate<-covariate[which(Valid!=0),]
  covariate <- as.matrix(covariate)

  #number of covariates
  P <- ncol(covariate)
  #sample size
  N1 <- nrow(GGUM.Data)

  ###### Data preparation
  #add delta and ind to response data
  GGUM.Data <- rbind(GGUM.Data, delindex)
  GGUM.Data <- rbind(GGUM.Data, ind)
  #reorder response data based on delta
  GGUM.Data <- GGUM.Data[, order(GGUM.Data[N1+2,],decreasing=FALSE)]
  #reordered dimension index
  ind <- GGUM.Data[N1+3,]
  #reordered delindex
  delindex <- GGUM.Data[(N1+1):(N1+2),]
  #reordered response data
  GGUM.Data <- GGUM.Data[-((N1+1):(N1+3)),]
  #### indicate which dimension does each item belong to while taking into consideration of missing data
  #0=missing
  #1=non-missing
  Missing <- matrix(NA,nrow=(ncol(delindex))*N1,ncol=1)
  MissPattern<-data.frame(Missing=((as.numeric(is.na(t(GGUM.Data)))*-1)+1),ID=seq(1,((ncol(delindex))*N1),1))
  NoMiss<-subset(MissPattern,Missing==1)
  #table(MissPattern$Missing)
  ind<-rep(ind,N1)[NoMiss$ID]

  #####long format response data
  Data<-suppressWarnings(edstan::irt_data(response_matrix =GGUM.Data))
  I<-dim(GGUM.Data)[1]
  J<-dim(GGUM.Data)[2]
  K=max(Data$y)

  #####number of ne, nu, po delta items
  ne=0
  nu=0
  po=0
  for (i in 1:(ncol(delindex))){
    if (delindex[2,i]<0){
      ne=ne+1
    } else if (delindex[2,i]==0){
      nu=nu+1
    } else if (delindex[2,i]>0){
      po=po+1
    }
  }

  #####long format person covariates
  PC1<-covariate[rep(seq_len(nrow(covariate)), each=trait), ]
  PC1 <- as.matrix(PC1)

  data_list<-list(
    I=Data$I,
    J=Data$J,
    N=Data$N,
    II=Data$ii,
    JJ=Data$jj,
    K=max(Data$y),
    I_ne=ne,
    I_nu=nu,
    I_po=po,
    I_NN=ne+nu,
    N_mis=(Data$I*Data$J-Data$N),
    y=Data$y,
    M=2*K,
    trait=trait,
    ind=ind,
    theta_mu=as.array(c(rep(0,trait))),
    P=P,
    PC1=PC1,
    ma=ma,
    va=va,
    mdne=mdne,
    mdnu=mdnu,
    mdpo=mdpo,
    vd=vd,
    mt=mt,
    vt=vt)

  #######################################################
  #######################################################
  #                                                     #
  #          Bayesian estimation of GGUM                #
  #                                                     #
  #######################################################
  #######################################################
  BZ_bmggum<-'

functions {

real GGUM(int y, real theta, real alpha, real delta, vector tau, int K, int M) {

     vector[K] prob;

     vector[K] nominator;

     for (a in 1:K) {

     nominator[a]=exp(alpha*((a-1)*(theta-delta)-sum(tau[1:a])))+exp(alpha*((M-a)*(theta-delta)-sum(tau[1:a])));

     }

     prob=nominator/sum(nominator);

     return categorical_lpmf(y|prob);
}
}

data {
  int<lower=1> K;
  int<lower=2> M;  // M=2K
  int<lower=1> I;
  int<lower=1> J;
  int<lower=1> N;
  int<lower=1> I_ne;
  int<lower=0> I_nu;
  int<lower=1> I_po;
  int<lower=1> I_NN;
  int<lower=0> N_mis;
  int<lower=1, upper=I> II[N];
  int<lower=1, upper=J> JJ[N];
  int<lower=0, upper=6> y[N];

  // user-defined priors
  real ma;
  real va;
  real mdne;
  real mdnu;
  real mdpo;
  real vd;
  vector<lower=-5,upper=0>[K] mt;
  real vt;

  int<lower=1,upper=I>trait;
  int<lower=1,upper=I>ind[N];
  vector[trait] theta_mu;
  int<lower=1> P;
  matrix[trait*J,P] PC1;
}

parameters {

  vector<lower=0,upper=4>[I] alpha;
  vector<lower=0,upper=5>[I_po]   delta_po;
  vector<lower=-5,upper=0>[I_ne]  delta_ng;
  vector<lower=-5,upper=5>[I_nu]  delta_nu;
  vector<lower=-5,upper=0>[K-1]   tau_raw;

  vector[trait] theta[J];
  cholesky_factor_corr[trait] L_Omega;
  matrix[P,trait] lambda;    // person covariate regression coefficient
}

transformed parameters{

  vector[I_NN] delta1;
  vector[I] delta;
  vector[K] tau;

  delta1=append_row(delta_ng,delta_nu);
  delta=append_row(delta1,delta_po);

  //for (i in 1:I){

    tau[1]=0;
  //}

  for (d in 2:K){

    tau[d]=tau_raw[d-1];
}
}


model {

    matrix[trait, trait] center_mu;  //PC1*lambda matrix
    row_vector[trait] center_mu0;    //first row of PC1*lambda matrix
    vector[trait] center_mu1;        //change center_mu0 row to column

    alpha ~    lognormal(ma,va);
    delta_ng ~ normal(mdne,vd);
    delta_nu ~ normal(mdnu,vd);
    delta_po ~ normal(mdpo,vd);
    L_Omega  ~ lkj_corr_cholesky(1);

    for (t in 1:(K-1)){
    tau_raw[t] ~ normal(mt[t],vt);
    }


	//theta~ multi_normal_cholesky(theta_mu,L_Omega);
    for (j in 1:J){
      center_mu=PC1[(trait*j-(trait-1)):trait*j,]*lambda;
      center_mu0=center_mu[1,];
      center_mu1=to_vector(center_mu0);
      theta[j]~ multi_normal_cholesky(theta_mu + center_mu1, L_Omega);
    }
     to_vector(lambda) ~ normal (0, 1);

for (n in 1:N){

  target += GGUM(y[n],theta[JJ[n],ind[n]],alpha[II[n]],delta[II[n]],tau[],K,M);

}
}

generated quantities{
matrix[trait,trait] Cor;
vector[N] log_lik;

Cor=multiply_lower_tri_self_transpose(L_Omega);

  for (n in 1:N) {
    log_lik[n] = GGUM(y[n],theta[JJ[n],ind[n]],alpha[II[n]],delta[II[n]],tau[],K,M);
  }
}
'

##################################################
#       Input response data estimation        #
##################################################

rstan::rstan_options(auto_write = TRUE,javascript = FALSE)

BZ_bmggum<-rstan::stan(model_code=BZ_bmggum,data=data_list,
                       iter=iter, chains=chains,cores=cores, warmup=warmup,
                       init=init, thin=thin,
                       control=list(adapt_delta=adapt_delta,max_treedepth=max_treedepth))

#####Extract some parameters
THETA<-rstan::summary(BZ_bmggum, pars = c("theta"), probs = c(0.025,0.5,0.975))$summary
Alpha_ES<-rstan::summary(BZ_bmggum, pars = c("alpha"), probs = c(0.025,0.5,0.975))$summary
Delta_ES<-rstan::summary(BZ_bmggum, pars = c("delta"), probs = c(0.025, 0.5,0.975))$summary
Tau_ES<-rstan::summary(BZ_bmggum, pars = c("tau_raw"), probs = c(0.025, 0.5,0.975))$summary
Cor_ES<-rstan::summary(BZ_bmggum, pars = c("Cor"), probs = c(0.025, 0.5,0.975))$summary
Lamda_ES<-rstan::summary(BZ_bmggum, pars = c("lambda"), probs = c(0.025, 0.5,0.975))$summary

#####reorder item parameters and data to the original input data order
delindex <- t(delindex)
estimates <- cbind(Alpha_ES, Delta_ES, delindex)
estimates <- estimates[order(estimates[,8+8+1],decreasing=FALSE),]
Alpha_ES <- estimates[,1:8]
Delta_ES <- estimates[,9:16]
rownames(Alpha_ES) <- c()
rownames(Delta_ES) <- c()
#change tau estimates to the same number of rows as alpha UM8
Tau_ES <- Tau_ES[rep(seq_len(nrow(Tau_ES)), nrow(delindex)), ]
# delindex1 <- delindex[rep(seq_len(nrow(delindex)), each = (option-1)), ]
# estimates1 <- cbind(Tau_ES, delindex1)
# estimates1 <- estimates1[order(estimates1[,8+1],decreasing=FALSE),]
# Tau_ES <- estimates1[,1:8]
delindex <- t(delindex)
GGUM.Data <- rbind(GGUM.Data, delindex)
GGUM.Data <- GGUM.Data[, order(GGUM.Data[N1+1,],decreasing=FALSE)]
GGUM.Data <- GGUM.Data[-((N1+1):(N1+2)),]

#####save estimated parameters to an R object
MGGUM.summary<-list(Theta.est=THETA,
                    Alpha.est=Alpha_ES,
                    Delta.est=Delta_ES,
                    Tau.est=Tau_ES,
                    Cor.est=Cor_ES,
                    Lamda.est=Lamda_ES,
                    Data=GGUM.Data,
                    Fit=BZ_bmggum,
                    Dimension=dimension)
}
  }

return(MGGUM.summary)
}


#' @title bmggum results extraction
#' @description This function extracts bmggum estimation results.
#' @param x bmggum returned object
#' @param pars Names of extracted parameters. They can be "theta" (person trait estimates), "alpha" (ttem discrimination parameters), "delta" (ttem location parameters), "tau" (ttem threshold parameters), "cor" (correlations among latent traits), "lambda" (regression coefficients linking person covariates to latent traits), "data" (GGUM.Data after deleting respondents who endorse the same response options across all items), "fit" (the stanfit object), and "dimension" (the input row vector mapping each item to each trait). Note that when the model is UM4 in which alpha is fixed to 1, the extracted alpha is a n*1 matrix where n equals to the number of items.
#' @return Selected results output
#' @examples
#' \dontrun{
#' theta <- extract.bmggum(mod, 'theta')
#' alpha <- extract.bmggum(mod, 'alpha')
#' cor <- extract.bmggum(mod, 'cor')
#' data <- extract.bmggum(mod, 'data')
#' mod <- bmggum(GGUM.Data, delindex, trait, ind, covariate, option)}
#' @export
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


#' @title bmggum Model fit
#' @description This function provides model fit statistics.
#' @param x bmggum returned object
#' @param index Model fit indices. They can be "waic", which is the widely applicable information criterion, "loo", which is the leave-one-out cross-validation, or "chisq.df", which is the adjusted chi-square degrees of freedom ratios for each trait separately that were introduced by Drasgow et al. (1995). The default is loo. Note that chisq.df can only be computed when the sample size is large. See documentation for loo and GGUM for more details.
#' @return Selected model fit statistics
#' @examples
#' \dontrun{
#' waic <- modfit.bmggum(mod, 'waic')
#' loocv <- modfit.bmggum(mod, 'loo')
#' chisq.df <- modfit.bmggum(mod, 'chisq.df')
#' mod <- bmggum(GGUM.Data, delindex, trait, ind, covariate, option)}
#' @export
modfit.bmggum <- function(x, index="loo"){

  x1 <- extract.bmggum(x, 'fit')
  log_lik=loo::extract_log_lik(x1, merge_chains = FALSE)
  rel_n_eff=loo::relative_eff(exp(log_lik))

  if (index=="chisq.df"){
    dimension <- extract.bmggum(x, 'dimension')
    data <- extract.bmggum(x, 'data')-1
    tau <- extract.bmggum(x, 'tau')
    tau <- tau[,1]
    alpha <- extract.bmggum(x, pars='alpha')
    alpha <- alpha[,1]
    delta <- extract.bmggum(x, pars='delta')
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


#' @title bmggum bayesian convergence diagnosis plotting function
#' @description This function provides plots including density plots, trace plots, and auto-correlation plots to aid model convergence diagnosis.
#' @param x bmggum returned object
#' @param pars Names of plotted parameters. They can be "theta", "alpha", "delta", "tau", "cor", "lambda", or a subset of parameters. See vignette for bmggum for more details.
#' @param plot Types of plots.They can be "density", "trace", or "autocorrelation".
#' @param inc_warmup Whether to include warmup iterations or not when plotting. The default is FALSE.
#' @return Selected plots for selected parameters
#' @examples
#' \dontrun{
#' bayesplot.bmggum(mod, 'alpha', 'density', inc_warmup=T)
#' bayesplot.bmggum(mod, 'delta', 'trace', inc_warmup=F)
#' bayesplot.bmggum(mod, 'tau', 'autocorrelation', inc_warmup=T)
#' mod <- bmggum(GGUM.Data, delindex, trait, ind, covariate, option)}
#' @export
bayesplot.bmggum <- function(x, pars, plot, inc_warmup=F){

  x <- extract.bmggum(x, 'fit')
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


#' @title bmggum item plotting function including observable response categories (ORCs)
#' @description This function provides item plots including observable response categories plots.
#' @param x bmggum returned object
#' @param items The items to be plotted. The default is all the items.
#' @return Selected ORC plots for selected items
#' @examples
#' \dontrun{
#' itemplot.bmggum(mod, items=1:5)
#' mod <- bmggum(GGUM.Data, delindex, trait, ind, covariate, option)}
#' @export
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

  alpha <- extract.bmggum(x, pars='alpha')
  delta <- extract.bmggum(x, pars='delta')
  alpha <- alpha[,1]
  delta <- delta[,1]
  tau <- extract.bmggum(x, pars='tau')
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


