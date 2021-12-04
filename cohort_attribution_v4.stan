
functions {

  /**
   * Return the Poisson-binomial log probability mass for the specified
   * count y and vector of probabilities theta.  The result is the log
   * probability of y successes in N = size(theta) independent
   * trials with probabilities of success theta[1], ..., theta[N].
   *
   * See:  https://en.wikipedia.org/wiki/Poisson_binomial_distribution
   *
   * @param y number of successes
   * @param theta vector of probabilities
   * @return Poisson-binomial log probability mass
   */

  real poisson_binomial_lpmf(int y, vector theta) {

    int N = rows(theta);
    matrix[N + 1, N + 1] alpha;
    vector[N] log_theta = log(theta);
    vector[N] log1m_theta = log1m(theta);

    if (y < 0 || y > N)
      reject("poisson binomial variate y out of range, y = ", y,
             " N = ", N);
    for (n in 1:N)
      if (theta[n] < 0 || theta[n] > 1)
        reject("poisson binomial parameter out of range,",
               " theta[", n, "] =", theta[n]);

    //if (N == 0)
    //  return y == 0 ? 0 : negative_infinity();

    // dynamic programming with cells
    // alpha[n + 1, tot + 1] = log prob of tot successes in first n trials
    alpha[1, 1] = 0;
    for (n in 1:N) {
      // tot = 0
      alpha[n + 1, 1] = alpha[n, 1] + log1m_theta[n];

      // 0 < tot < n
      for (tot in 1:min(y, n - 1))
        alpha[n + 1, tot + 1]
            = log_sum_exp(alpha[n, tot] + log_theta[n],
                          alpha[n, tot  + 1] + log1m_theta[n]);

      // tot = n
      if (n > y) continue;
      alpha[n + 1, n + 1] = alpha[n, n] + log_theta[n];
    }
    return alpha[N + 1, y + 1];

  }

}

data {

  // the total number of installs in EA server on day A and B
  real<lower=0> I_ea[2];
  //number of campaigns with impressions for day A
  int<lower=0> num_cmps;
  // vector of total number of installs in the postback (a scalar for each campaign)
  // note: the total number of installs can reflect installs from both day a and day b and
  // so is a relatively weak constraint
  int<lower=0> K_c_skan[num_cmps];
  
  vector[K_c_skan[1]] bernoulli_probs_c1;
  vector[K_c_skan[2]] bernoulli_probs_c2;
  vector[K_c_skan[3]] bernoulli_probs_c3;
  
  //negative binomial prior means for inst_c_a

  // int<lower=0> prior_mu_c_a[num_cmps];
  // int<lower=0> prior_mu_c_a[num_cmps];
  real<lower=0> prior_mu_c_a[num_cmps];
  real<lower=0> prior_var_c_a[num_cmps];
  
  real<lower=0> org_var;
}

transformed data {
  int <lower=0> K_skan_sum = sum(K_c_skan);
  

}

parameters {
  // number of installs for campaigns in C on day A
  real<lower=0, upper=K_c_skan[1]> inst_c1_a;
  real<lower=0, upper=K_c_skan[2]> inst_c2_a;
  real<lower=0, upper=K_c_skan[3]> inst_c3_a;

  // residual variance
  //real<lower=0> sigma;

}

 transformed parameters {
  //number of organic installs on day a
  real<lower=0> ORG[2];
  
  //total number of paid installs across campaigns for day a 
 // real<lower=0,upper=K_skan_sum> inst_c_sum[2];

  
  
  //inst_c_a_sum = sum(inst_c_a);
  //inst_c_sum[1] = inst_c1_a+inst_c2_a+inst_c3_a;
  //inst_c_sum[2] = (K_c_skan[1]-inst_c1_a)+(K_c_skan[2]-inst_c2_a)+(K_c_skan[3]-inst_c3_a);

  //inst_c_sum[1] = inst_c1_a+inst_c2_a+inst_c3_a;

  // ORG[1] = I_ea[1]-inst_c_sum[1];
  // ORG[2] = I_ea[2]-inst_c_sum[2];
  ORG[1] = I_ea[1]-(inst_c1_a+inst_c2_a+inst_c3_a);
  ORG[2] = I_ea[2]-((K_c_skan[1]-inst_c1_a)+(K_c_skan[2]-inst_c2_a)+(K_c_skan[3]-inst_c3_a));
  

}

model {
  
  //Likelihood for total installs in EA server on day a
  //I_ea[1] ~ normal(ORG[1], sigma)T[0,I_ea[1]];
  //I_ea[2] ~ normal(ORG[2], sigma)T[0,I_ea[2]];
  I_ea[1] ~ normal(ORG[1], prior_mu_c_a[1]+prior_mu_c_a[2]+prior_mu_c_a[3])T[0,I_ea[1]];
  I_ea[2] ~ normal(ORG[2], (K_c_skan[1]-prior_mu_c_a[1])+(K_c_skan[2]-prior_mu_c_a[2])+(K_c_skan[3]-prior_mu_c_a[3]))T[0,I_ea[2]];
  
  //target +=  normal_lpdf(I_ea[1] | ORG[1], inst_c1_a+inst_c2_a+inst_c3_a);//T[0,I_ea_a];
  //target +=  normal_lpdf(I_ea[2] | ORG[2], (K_c_skan[1]-inst_c1_a)+(K_c_skan[2]-inst_c2_a)+(K_c_skan[3]-inst_c3_a));//T[0,I_ea_a];

  //prior for variance
  //sigma ~ normal(prior_mu_c_a[1]+prior_mu_c_a[2]+prior_mu_c_a[3], mean([prior_var_c_a[1],prior_var_c_a[2],prior_var_c_a[3]]) );
  
  //prior for total installs in a day
  //inst_c_a ~ normal(prior_mu_c_a, 10);
  //inst_c1_a ~ normal(prior_mu_c_a[1], prior_var_c_a[1]);
  //inst_c2_a ~ normal(prior_mu_c_a[2], prior_var_c_a[2]);
  //inst_c3_a ~ normal(prior_mu_c_a[3], prior_var_c_a[3]);
 // inst_c1_a ~ poisson_binomial_lpmf(2, bernoulli_probs);
  //inst_c1_a ~ poisson_binomial_lpmf(inst_c1_a, bernoulli_probs);
  //inst_c1_a ~ poisson_binomial_lpmf(inst_c1_a, theta);
  //inst_c1_a ~ test_lpmf(inst_c1_a, theta);

  inst_c1_a ~ normal(prior_mu_c_a[1], prior_var_c_a[1]);
  inst_c2_a ~ normal(prior_mu_c_a[2], prior_var_c_a[2]);
  inst_c3_a ~ normal(prior_mu_c_a[3], prior_var_c_a[3]);

  //target += poisson_binomial_lpmf( prior_mu_c_a[1] | bernoulli_probs_c1);
  //target += poisson_binomial_lpmf( prior_mu_c_a[2] | bernoulli_probs_c2);
 // target += poisson_binomial_lpmf( prior_mu_c_a[3] | bernoulli_probs_c3);

  //target += test_lpmf( 2 | [3,2]);

 // target += normal_lpdf( inst_c1_a | 3, bernoulli_probs);

  //inst_c2_a ~ normal(prior_mu_c_a[2], prior_var_c_a[2]);
  //inst_c3_a ~ normal(prior_mu_c_a[3], prior_var_c_a[3]);
  
  //prior for organic installs 
  ORG ~ normal(110, org_var);
  //target += normal_lpdf(ORG | 110, org_var);

}
