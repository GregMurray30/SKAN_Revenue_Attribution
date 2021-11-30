

data {

  // the total number of installs in EA server on day A
  real<lower=0> I_ea_a;
  //number of campaigns with impressions for day A
  int<lower=0> num_cmps;
  // vector of total number of installs in the postback (a scalar for each campaign)
  // note: the total number of installs can reflect installs from both day a and day b and
  // so is a relatively weak constraint
  int<lower=0> K_c_skan[num_cmps];
  
  //negative binomial prior means for inst_c_a
  real<lower=0> prior_mu_c_a[num_cmps];
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
  real<lower=0> sigma;

}

 transformed parameters {
  //number of organic installs on day a
  real<lower=0> ORG_a;
  
  //total number of paid installs across campaigns for day a 
  real<lower=0,upper=K_skan_sum> inst_c_a_sum;

  ORG_a = I_ea_a-(inst_c1_a+inst_c2_a+inst_c3_a);
  
  //inst_c_a_sum = sum(inst_c_a);
  inst_c_a_sum = inst_c1_a+inst_c2_a+inst_c3_a;

}

model {
  
  //Likelihood for total installs in EA server on day a
  I_ea_a ~ normal(ORG_a, sigma)T[0,I_ea_a];
  
  //prior for variance
  sigma ~ normal(10, 10);
  //prior for total installs in a day
  //inst_c_a ~ normal(prior_mu_c_a, 10);
  inst_c1_a ~ normal(prior_mu_c_a[1], 15);
  inst_c2_a ~ normal(prior_mu_c_a[2], 15);
  inst_c3_a ~ normal(prior_mu_c_a[3], 15);

  //prior for organic installs 
  ORG_a ~ normal(100, 10);
  
}
