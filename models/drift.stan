functions {
  real partial_sum(real[] slice_beta,
                   int start, int end,
                   vector l_snpsd,
                   vector f,
                   vector l_fp,
                   vector l_mfp) {
    return normal_lpdf(slice_beta |
                               0,
                               exp(l_snpsd[start:end])) + 
                               beta_lpdf(f[start:end]|exp(l_fp[start:end]),exp(l_mfp[start:end]));
  }
}
data {
  int<lower=0> J; // number of snps
   real beta[J]; // effect sizes
   vector<lower=0,upper=1>[J] f; // minor allele frequences
   real<lower=0> Fst; // Known value of Fst
}
transformed data{
  real log_Fst = log(Fst);
  real lm1_Fst = log(1-Fst);
}
parameters {
  real<lower=-2,upper=2> S; 
  real<lower=0, upper=2> sigma_b;
  vector<lower=0,upper=0.5>[J] p; // true minor allele frequences
}
transformed parameters {
  vector[J] l_snpsd;
  vector[J] log_Fp;
  vector[J] log_F1mp;
  vector[J] l_p;
  vector[J] l_mp;
  l_p=log(p);
  l_mp=log1m(p);
  log_F1mp = lm1_Fst -log_Fst + l_mp;
  log_Fp = lm1_Fst -log_Fst + l_p;
  
  l_snpsd=(0.5*S)*(l_p+l_mp) +log(sigma_b);
}
model {
  // vectorised likelihood
   //f ~ beta( exp(log_Fp), exp(log_F1mp));
   //beta ~ normal(0,exp(l_snpsd));
  
   target += reduce_sum(partial_sum, beta, 1,
                        l_snpsd, f, log_Fp, log_F1mp);
}