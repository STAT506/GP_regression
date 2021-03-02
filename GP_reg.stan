data {
  int<lower=0> N; // number of data points
  vector[N] y; // response
  matrix[N,N] dist; // distance matrix
  vector[N] x; // covariate
  int<lower=0> N_preds;
  matrix[N_preds, N_preds] dist_preds;
  matrix[N, N_preds] dist_12; 
  vector[N_preds] x_preds; // covariate
}

parameters {
  real<lower = 0.25, upper = 9> phi;
  real<lower = 0> sigmasq;
  real<lower = 0> tausq;
  real beta;
}

transformed parameters{
  vector[N] mu_vec;
  vector[N] tausq_vec;
  corr_matrix[N] Sigma;
  
  for(i in 1:N) mu_vec[i] = x[i] * beta;
  for(i in 1:N) tausq_vec[i] = tausq;
  
  for(i in 1:(N-1)){
   for(j in (i+1):N){
     Sigma[i,j] = exp((-1)*dist[i,j]/ phi);
     Sigma[j,i] = Sigma[i,j];
   }
 }
 for(i in 1:N) Sigma[i,i] = 1;

}

model {
  y ~ multi_normal(mu_vec ,sigmasq * Sigma + diag_matrix(tausq_vec));
  phi ~ inv_gamma(10, 10);
  sigmasq ~ inv_gamma(10, 10);
  tausq ~ inv_gamma(10, 2);
  beta ~ normal(0, 10);
}
 

generated quantities {
  vector[N_preds] y_preds;
  vector[N] y_diff;
  vector[N_preds] mu_preds;
  corr_matrix[N_preds] Sigma_preds;
  vector[N_preds] tausq_preds;
  matrix[N, N_preds] Sigma_12;

  for(i in 1:N_preds) tausq_preds[i] = tausq;
  for(i in 1:N_preds) mu_preds[i] = x_preds[i] * beta;
  for(i in 1:N) y_diff[i] = y[i] - x[i] * beta;
  

  for(i in 1:(N_preds-1)){
   for(j in (i+1):N_preds){
     Sigma_preds[i,j] = exp((-1)*dist_preds[i,j]/ phi);
     Sigma_preds[j,i] = Sigma_preds[i,j];
   }
 }
 for(i in 1:N_preds) Sigma_preds[i,i] = 1;
 
   for(i in 1:(N)){
   for(j in (1):N_preds){
     Sigma_12[i,j] = exp((-1)*dist_12[i,j]/ phi);
   }
 }

 y_preds = multi_normal_rng(mu_preds + (sigmasq * Sigma_12)' * inverse(sigmasq * Sigma) * (y_diff), sigmasq * Sigma_preds + diag_matrix(tausq_preds) - (sigmasq * Sigma_12)' * inverse(sigmasq * Sigma + diag_matrix(tausq_vec)) * (sigmasq * Sigma_12) );
}
