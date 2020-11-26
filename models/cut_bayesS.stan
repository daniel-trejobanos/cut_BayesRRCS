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
}
transformed data{
  vector[J] l_f;
  vector[J] l_mf;
  vector[J] l_term;
  real A;
  real l_a;
  l_f = log(f);
  l_mf = log1m(f);
  l_term =log(2) + l_f+l_mf;
  A=sum(2*f.*(1-f));
  l_a = log(A);
}
parameters {
  real<lower=0> sigma;
  real<lower=-2,upper=2> S_coef; // no need to specify S uniform
}
transformed parameters{
  real var_g;
  var_g = exp((1+S_coef)*l_a + 2*log(sigma));
}
model {
  vector[J] l_snpsd;
  real l_sigma;
   S_coef ~ normal(0,sqrt(1));
   sigma ~ normal(0,1);
//  l_sigma=log(sqrt(0.57));
   l_sigma=log(sigma);
  l_snpsd= 0.5*S_coef*(l_term) + l_sigma ; 
    //beta ~ normal(0,exp(l_sigma-0.5*log(J))); // vectorised likelihood
  target += reduce_sum(partial_sum,beta, 425000,
                     l_snpsd);

}
