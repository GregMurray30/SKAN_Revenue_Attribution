

data {

  // the total number of installs in EA server on day A
  real<lower=0> I_ea_a;
  //number of campaigns with impressions for day A
  int<lower=0> num_cmps;
  // vector of total number of installs in the postback (a scalar for each conversion value)
  // note: the total number of installs can reflect installs from both day a and day b and
  // so is a relatively weak constraint
  int<lower=0> K_skan;
  
  //bernoulli prior means for inst_c_a
  real<lower=0> prior_mu_c_a[num_cmps];
}


parameters {
  // number of installs for campaign C on day A
  real<lower=0, upper=K_skan> inst_c_a[num_cmps];


  // residual variance
  real<lower=0> sigma;

}

 transformed parameters {
  //number of organic installs on day a
  real<lower=0> ORG_a;
  real<lower=0,upper=K_skan> inst_c_a_sum;

  ORG_a = I_ea_a-sum(inst_c_a);
  
  inst_c_a_sum = sum(inst_c_a);
  
}

model {
  
  //Likelihood for total installs in EA server on day a
  I_ea_a ~ normal(ORG_a, sigma);
  
  //prior for variance
  sigma ~ normal(10, 10);
  //prior for total installs in a day
  inst_c_a ~ normal(prior_mu_c_a, 5);
  //prior for organic installs 
  ORG_a ~ normal(100, 10);
  
}
