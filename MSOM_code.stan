data {
  int<lower=0> J;  // Number of sites
  int<lower=0> K[J];  // Number of occasions for each site
  int<lower=0> n1;  // Number of observed species
  int<lower=0> n0;  // Number of unobserved species
  int<lower=0> n;  // Total number of species (n0+n1)
  int<lower=0> PA;  // Number of parks
  int<lower=0> park[J];  // Index for each park
  
  real pasize[PA];  // Park size
  real punish[PA];  // Punishment factor
  real reach[PA];  // Reach factor
  real elev[J];  // Elevation
  real pop[J];  // Population
  real dist[J];  // Distance
  real hour[J,max(K)];  // Camera hour
  real angle[J,max(K)];  // Camera angle
  
  int<lower=0> Y[J,max(K),n];  // Observation data with augmentation
}

parameters {
  // a0, b0
  real a0[n,PA];
  real b0[n];
  real mu_b0;
  real<lower=0> tau_b0;
  
  // a1-a6
  real a1[n];
  real a2[n];
  real a3[n];
  real a4[n];
  real a5[n];
  real a6[n];  

  real mu_a1;
  real mu_a2;
  real mu_a3;
  real mu_a4;
  real mu_a5;
  real mu_a6;
  
  real<lower=0> tau_a1;
  real<lower=0> tau_a2;
  real<lower=0> tau_a3;
  real<lower=0> tau_a4;
  real<lower=0> tau_a5;
  real<lower=0> tau_a6;
  
  // b1-b2
  real b1[n];
  real b2[n];
  
  real mu_b1;
  real mu_b2;
  
  real<lower=0> tau_b1;
  real<lower=0> tau_b2;
  // park
  real a_pa[n, PA];
  real<lower=0> tau_pa[n];
  
  // omega
  real<lower=0, upper=1> omega;
  
  // p_fit & p_fit_new
  //real<lower=0> p_fit;
  //real<lower=0> p_fit_new;
}

transformed parameters {
  real mu_pa[n, PA];
  for (i in 1:n) {
    for (f in 1:PA) {
      mu_pa[i, f] = a0[i, f] + a4[i] * pasize[f] + a5[i] * punish[f] + a6[i] * reach[f];
    }
  }
  
  real<lower=0, upper=1> psi[J,n];
  real<lower=0, upper=1> p[J,max(K),n];
  for (j in 1:J) {
    for (i in 1:n) {
      psi[j, i] = inv_logit(a_pa[i, park[j]] + a1[i] * elev[j] + a2[i] * pop[j] + a3[i] * dist[j]);

      for (k in 1:K[j]) {
        p[j, k, i] = inv_logit(b0[i] + b1[i] * hour[j, k] + b2[i] * angle[j, k]);
      }
    }
  }
}

model {
  int w[n];
  int Z[J,n];
  int Ynew[J,max(K),n];
  real Nsite[J];
  int N0;
  int N;
  // a0, b0
  inv_logit(mu_b0) ~ uniform(0, 1);
  tau_b0 ~ gamma(0.1, 0.1);
  // a1-a6
  mu_a1 ~ normal(0, 0.01);
  mu_a2 ~ normal(0, 0.01);
  mu_a3 ~ normal(0, 0.01);
  mu_a4 ~ normal(0, 0.01);
  mu_a5 ~ normal(0, 0.01);
  mu_a6 ~ normal(0, 0.01);
  tau_a1 ~ gamma(0.1, 0.1);
  tau_a2 ~ gamma(0.1, 0.1);
  tau_a3 ~ gamma(0.1, 0.1);
  tau_a4 ~ gamma(0.1, 0.1);
  tau_a5 ~ gamma(0.1, 0.1);
  tau_a6 ~ gamma(0.1, 0.1);
  // b1-b2
  mu_b1 ~ normal(0, 0.01);
  mu_b2 ~ normal(0, 0.01);
  tau_b1 ~ gamma(0.1, 0.1);
  tau_b2 ~ gamma(0.1, 0.1);
  // omega
  omega ~ uniform(0, 1);
  
  for (i in 1:n) {
    b0[i] ~ normal(mu_b0, tau_b0);
    
    a1[i] ~ normal(mu_a1, tau_a1);
    a2[i] ~ normal(mu_a2, tau_a2);
    a3[i] ~ normal(mu_a3, tau_a3);
    a4[i] ~ normal(mu_a4, tau_a4);
    a5[i] ~ normal(mu_a5, tau_a5);
    a6[i] ~ normal(mu_a6, tau_a6);
    
    b1[i] ~ normal(mu_b1, tau_b1);
    b2[i] ~ normal(mu_b2, tau_b2);
    
    w[i] ~ bernoulli(omega);
    
    sqrt(1/tau_pa[i]) ~ uniform(0, 10);

    for (f in 1:PA) {
      a0[i, f] ~ uniform(-10, 10);
      a_pa[i, f] ~ normal(mu_pa[i, f], tau_pa[i]);
      //mu_pa[i, f] = a0[i, f] + a4[i] * pasize[f] + a5[i] * punish[f] + a6[i] * reach[f];
    }
  }

  for (j in 1:J) {
    for (i in 1:n) {
      //psi[j, i] = inv_logit(a_pa[i, park[j]] + a1[i] * elev[j] + a2[i] * pop[j] + a3[i] * dist[j]);
      Z[j, i] ~ bernoulli(psi[j, i] * w[i]);

      for (k in 1:K[j]) {
        //p[j, k, i] = inv_logit(b0[i] + b1[i] * hour[j, k] + b2[i] * angle[j, k]);
        Y[j, k, i] ~ bernoulli(p[j, k, i] * Z[j, i]);
        Ynew[j, k, i] ~ bernoulli(p[j, k, i] * Z[j, i]);
        // Calculate the discrepancy measure
        
        //Create simulated dataset to calculate the Bayesian p-value
        //d[j,k,i] = abs(X[j,k,i] - mu.p[j,k,i]) 
        //dnew[j,k,i] = abs(Xnew[j,k,i] - mu.p[j,k,i]) 
        //d2[j,k,i] = pow(d[j,k,i],2)  
        //dnew2[j,k,i] = pow(dnew[j,k,i],2) 
        
      //dsum[j,i] = sum(d2[j,1:K[j],i]) 
      //dnewsum[j,i] = sum(dnew2[j,1:K[j],i])
      
        //p_fit = sum(pow(Y[j, 1:K[j], 1:n] - p[j, k, i] * Z[j, i], 2));
        //p_fit_new = sum(pow(Ynew[j, 1:K[j], 1:n] - p[j, k, i] * Z[j, i], 2));
      }
    }
  }

  // Calculate site richness
  for (j in 1:J) {
    Nsite[j] = dot_product(Z[j,1:n], w[1:n]);
  }

  // Calculate total estimated richness
  N0 = sum(w[(n1 + 1):n]);
  N = n1 + N0;
}
