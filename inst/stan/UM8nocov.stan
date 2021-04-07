//
// This Stan program defines a UM8 model with no covariates

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
