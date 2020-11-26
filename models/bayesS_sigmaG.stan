functions {
  real partial_sum(real[] slice_beta,
                   int start, int end,
                   vector l_sd) {
    return normal_lpdf(slice_beta | 0, exp(l_sd[start:end]));
  }
}
data {
  int<lower=0> J; // number of snps
  real beta[J]; // effect sizes
  vector<lower=0,upper=1>[J] f; // minor allele frequences
  real<lower=0> sigma_g;
}
transformed data{
  vector[J] l_f;
  vector[J] l_mf;
  vector[J] l_term;
  real A;
  real l_sigma;
 A=log(sum(2*f.*(1-f)));
  l_f = log(f);
  l_mf = log1m(f);
  l_term =log(2) + l_f+l_mf;
  l_sigma = log(sigma_g);
}
parameters {
  real S_coef; // no need to specify S uniform
}
model {
  vector[J] l_snpsd;
  S_coef ~ normal(0,1);
  l_snpsd= 0.5*S_coef*(l_term) + l_sigma - 0.5*(1+S_coef)*A ; 
  target += reduce_sum(partial_sum,beta, 425000,
                     l_snpsd);

}
