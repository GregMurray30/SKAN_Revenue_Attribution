# SKAN_Revenue_Attribution
Cohort Install Attribution Model

Author: Greg Murray, Senior Data Scientist EA Mobile

INTRODUCTION

The goal for this work is to attribute the installs and revenue from the SKAdNetwork (SKAN) postbacks to a cohort date(s). This is needed in order to train the elasticity model which will be used to inform User Acquisition spend decisions. It is vital that the estimate of the installs and revenue for each day and cohort, respectively, be as accurate as possible for the elasticity model to provide useful spend signals. 

This attribution will be accomplished in two parts. First, the installs must be attributed to a date. Second, once the installs are attributed, daily cohorts’ DX revenue can be calculated by combining the revenue estimates associated with the corresponding conversion value. 

 This model will leverage the knowledge of SKAN along with the postback conversion values and batch send times to infer the marginal probabilities of the install distributions for each campaign-conversion_value (Part I). Once calculated, we can use those marginal probabilities as inputs to a second model that uses maximum likelihood estimation techniques to find the most probable install ratio (Part II). More specifically, we will utilize the knowledge that:
	the postback chain “breaks” exactly 24 hours after the last conversion event if no other conversion value increase occurs to reset the 24-hour timer
	the conversion value indicates the days since install date in the first two bits (up to 3 days after install date)
	the batches are sent every 24 hours and include all postbacks that occurred since the last batch send and none of the postbacks that occurred before the previous batch send*
	the postbacks contain the total number of installs for every campaign-conversion_value pair
	the conversion values indicate a range of revenue (spender buckets)

It is preferrable that the model in Part II be used in production, but in case any of the model’s assumptions are egregiously violated, or time/computational constraints render that model impractical in the short term, the model in Part I can be used on its own. It should be noted, however, that the model in Part I will likely not provide powerful signals because it is only able to leverage data with low information content (SKAN send times).

Both model derivations deal only with installs as cohorts’ DX revenue can be inferred once installs for each conversion value are attributed to a date.

*This is not entirely certain but is the consensus among most familiar with SKAN at EA, Singular and Facebook

PART I
Marginal Cohort-Install Distributions

In this section we will derive the Bernoulli probability distribution over possible install dates from the batch send times based on a number of elements and axioms defined below.

Definitions
	D_install is a discrete random variable representing day of install 
	D_last is a discrete random variable representing day of last event 
	D_chain is a discrete random variable representing day of chain break 
	d_skan is the day of the SKAN postback receipt (deterministic)
	T_install is a continuous random variable representing time (seconds) of install [0, 86400) 
	T_last is a continuous random variable representing time (seconds)  of last event  [0, 86400)
	T_chain is a continuous random variable representing time (seconds) of chain break [0, 86400)
	t_skan is the time (seconds) of SKAN postback batch receipt [0, 86400) (deterministic)
	r is the number of days from retention bits in the conversion value [0-3] (deterministic)
	h is the threshold in seconds that marks the start of a new date [0, 86400); assumed to be 86400 (24:00) (deterministic)
	n is the conversion value [0-63] (deterministic)


Axioms
	D_install= D_last-r
	D_last= D_chain-1
	d_skan-1≤D_chain≤ d_skan
	t_x=□( t_x  mod 86400)
	t_skan= T_chain+X, X ~ U[0,86400)*
	p(D_install= d_skan-r-2)+p(D_install=d_skan-r-1 )=1
	T_chain,t_skan< h
	h=86400

*Assumption



D_j □(∶=) D_install+r+1
(i-viii.)
 D_j=d_skan-1= D_chain  □(⇔┬  )  〖 T〗_chain>t_skan
□(⇒┬  ) p(D_j= d_skan-1)= 〖p( T〗_chain>t_skan)
(v., vii., viii.)
=〖f_unif ( t_skan<T〗_chain<h)
=|h-t_skan |/86400

where f(x) is the probability density function for a uniform distribution. So,
 p( D_install= d_skan-r-2┤|  t_skan) = |h-t_skan |/86400.


PART II
Maximum Posterior Estimation

Having inferred a Bernoulli probability distribution over possible install dates from the batch send times, we can go one step further and combine our estimated distribution with the data in the EA server by simulating the conversion values for our users there. To combine them in this way, we have to be cognizant of the co-dependencies of all the parameters and random variables involved. The distribution of installs for any given campaign is conditionally dependent on the probability distribution of installs from all other campaigns (for that game and conversion value) during that date range, and directly dependent on:
	 the installs in the EA server on that day with the corresponding simulated conversion value (constraint)
	 the total count of installs in that campaign’s postback (constraint) 
Because of the enormous number of co-dependencies, estimating the probability directly is intractable. 

Fortunately, there is a way of estimating the ratio of installs between day_a and day_b without directly estimating the probability distribution. To do so, we use a discrete R.V. (ie: negative binomial or Poisson) to model the conditional probability that: install_count on day_a = x, given the installs in the EA server. This avoids the more consistent (for small sample sizes) but nearly impossible estimation of the marginal probability that install_date = day_a with a Bernoulli random variable as was done in Part I. 

Because we are mainly interested in the most likely value of the parameters (campaigns’ installs on day a & b campaign) given the install counts in the EA server on both days*, and they are iterable quanta (computationally feasible), we can simply maximize the product of the likelihood and prior (maximum posterior estimation) to turn the problem into one of constrained optimization. The only stipulation is that in order to compute the likelihood term, we must calculate the probability distribution of organic install counts ahead of time. Luckily, we have estimates of those distributions as we are currently able to distinguish paid and organic installs and can supplement those estimates with k-factor analysis currently underway.

*posterior

Definitions
	d_a □(∶=)  d_skan-r-2
	d_b □(∶=)  d_skan-r-1
	i_(x_(t_skan ) n)^EA is the count of installs in EA server on day x, where x∈{a,b}, and 〖time〗_EA>t_skan if day=a and 〖time〗_EA≤t_skan if day=b, for conversion value n. 

	I_(x_(t_skan ) n)^c is a negative binomial random variable of count of installs on day x where 〖time〗_EA>t_skan if day =a, and 〖time〗_EA≤t_skan if day =b, for campaign c and conversion value n

	〖ORGANIC〗_xn is a discrete random variable of count of organic installs on day x for conversion value n
	K_n^c    is the total count of installs in the SKAN postback for campaign c and conversion value n
	C_n is the total count of campaigns currently active for the game being modelled and conversion value n

*Subscript notation will omit 〖t_skan〗_ and n as the time operator is implied by day=a or day=b (>t_skan or <=t_skan), and the conversion value n is constant for each model but is important to note in the definitions for context.

Constraints
	〖organic〗_a+∑_(c=1)^C▒i_a^c = i_a^EA
	〖organic〗_b+∑_(c=1)^C▒i_b^c = i_b^EA
	i_a^c+i_b^c=k^c   ,∀c 

Axioms
	p(I_x^1= i_x^1)⊥〖p(I〗_x^2= i_x^2) ⊥⋯⊥ 〖p(I〗_x^c= i_x^c),∀c ,∀x,c ∈C,x∈(a,b)*
	p(I_x^c= i_x^c  ┤|  i_a^EA,i_b^EA  )≠ p(I_x^c= i_x^c),∀c ,∀x, c ∈C, x∈(a,b)
	 I_x^EA 〖=ORGANIC〗_x+∑_(c=1)^C▒〖 I_x^c 〗 ,∀c ,∀x, c ∈C,x∈(a,b)
	 i_x^EA 〖=organic〗_x+∑_(c=1)^C▒〖 i_x^c 〗 ,∀c ,∀x, c ∈C,x∈(a,b)
	〖p(organic〗_a) ⊥〖p(organic〗_b)
	∑_(i=0)^(K^c)▒〖p(I_a^c= i_a^c  )+p(I_b^c= K^c-i_a^c ) 〗=1,∀c ,c ∈C
	〖p(organic〗_x)⊥ p(I_x^c= i_x^c  ),∀c ,∀x, c ∈C, x∈(a,b)*
	p(I_x^c= i_x^c  ,I_x^c'= i_x^c'  ┤|  i_x^EA  ) ≠ p(I_x^c= i_x^c  ┤|  i_x^EA  )*p(I_x^c'= i_x^c'  ┤|  i_x^EA  ),∀c ,∀x, c ∈C,x∈(a,b)

*Assumption


We first define our posterior as seen below in purple. This arises naturally from the business problem and formulation of the model. Because of the conditional dependence of all campaigns’ install counts (for each day) on the install counts in the EA server, estimating this probability distribution directly is practically impossible as it requires estimating the conditional probabilities of all possible combinations of campaigns! Consequently, even if we had access to a way of directly estimating those conditional probabilities (we don’t), the computational complexity alone renders a direct estimation of the posterior infeasible. 
	Instead, we can leverage the clever simplicity of Bayes theorem so need only estimate the right side of the proportionality below,


Posterior ∝ L*π.

Posterior□(∶=)  p( I_a^1= i_a^1,〖I_b^1= i_b^1,I〗_a^2= i_a^2,〖 I〗_b^2= i_b^2,…,I_a^c= i_a^c,I_b^c= i_b^c   |  〖I_a^EA=i〗_a^EA,I_b^EA=i_b^EA )


π□(∶=)  p( I_a^1= i_a^1,〖I_b^1= i_b^1,I〗_a^2= i_a^2,〖 I〗_b^2= i_b^2,…,I_a^c= i_a^c,I_b^c= i_b^c )

(i., Negative Binomial PMF)
□(=) ∏_(c=1)^C▒〖C_(〖(k〗^c-i_a^c))^(k^c )*p(I_a^c= i_a^c )^(i_a^c )*p(I_b^c= i_b^c )^(i_b^c ) 〗  

Finally, we end up with a calculable prior since we can plug in the marginal probabilities calculated in Part I. However, it is preferable to avoid potential rounding errors with arbitrarily small products of fractions – especially in the case of many campaigns or large values of k^c, so we take the negative log prior instead,
(vi., monotonicity of natural log)

□(π_log=-) ∑_(c=1)^C▒〖log(C_(〖(k〗^c-i_a^c))^(k^c ))+log⁡(i_a^c*p(I_a^c= i_a^c )^  )+log⁡(〖〖(k〗^c-i〗_a^c )*(1-p〖(I_a^c= i_a^c ))〗^  〗.


The likelihood derivation is slightly more involved,

L□(∶=)  p(〖I_a^EA=i〗_a^EA,I_b^EA=i_b^EA   |  I_a^1= i_a^1,〖I_b^1= i_b^1,I〗_a^2= i_a^2,〖 I〗_b^2= i_b^2,…,I_a^c= i_a^c,I_b^c= i_b^c )

(ii., iii.)

= p(〖ORGANIC〗_a+〖∑_(c=1)^C▒〖 I_a^c 〗  =organic〗_a+∑_(c=1)^C▒〖 i_a^c  〗,〖 ORGANIC〗_b+〖∑_(c=1)^C▒〖 I_b^c 〗  =organic〗_b+∑_(c=1)^C▒〖 i_b^c  〗   ┤|  〖 I〗_a^1= i_a^1,I_b^1= i_b^1,I_a^2= i_a^2  〖,I〗_b^2= i_b^2,〖…,I_a^c= i_a^c,I〗_b^c= i_b^c  )

 ≡p(〖ORGANIC〗_a 〖  =organic〗_a,∑_(c=1)^C▒〖 I_a^c 〗=∑_(c=1)^C▒〖 i_a^c 〗  ,〖  ORGANIC〗_b 〖=organic〗_b,∑_(c=1)^C▒〖 I_b^c 〗=∑_(c=1)^C▒〖 i_b^c  〗   ┤|  〖 I〗_a^1= i_a^1,I_b^1= i_b^1,I_a^2= i_a^2  〖,I〗_b^2= i_b^2,〖…,I_a^c= i_a^c,I〗_b^c= i_b^c  )
(vii.)
= p(〖ORGANIC〗_a 〖  =organic〗_a,〖 ORGANIC〗_b 〖=organic〗_b  ┤|  〖 I〗_a^1= i_a^1,I_b^1= i_b^1,I_a^2= i_a^2  〖,I〗_b^2= i_b^2,〖…,I_a^p= i_a^p,I〗_b^p= i_b^p  ) * p(∑_(c=1)^C▒〖 I_a^c 〗=∑_(c=1)^C▒〖 i_a^c 〗,∑_(c=1)^C▒〖 I_b^c 〗=∑_(c=1)^C▒〖 i_b^c  〗  ┤|  〖 I〗_a^1= i_a^1,I_b^1= i_b^1,I_a^2= i_a^2  〖,I〗_b^2= i_b^2,〖…,I_a^c= i_a^c,I〗_b^c= i_b^c  )

(p(X=x|X=x)=1)
= p(〖ORGANIC〗_a 〖  =organic〗_a,〖 ORGANIC〗_b 〖=organic〗_b  ┤|  〖 I〗_a^1= i_a^1,I_b^1= i_b^1,I_a^2= i_a^2  〖,I〗_b^2= i_b^2,〖…,I_a^c= i_a^c,I〗_b^c= i_b^c  )

(iii.)
= p(〖ORGANIC〗_a=i_a^EA-∑_(c=1)^C▒〖 i_a^c 〗,〖ORGANIC〗_b=i_b^EA-∑_(c=1)^C▒〖 i_b^c 〗   ┤|   I_a^1= i_a^1,I_b^1= i_b^1,I_a^2= i_a^2  〖,I〗_b^2= i_b^2,〖…,I_a^c= i_a^c,I〗_b^c= i_b^c  )

(v.)

L = p(〖ORGANIC〗_a=i_a^EA-∑_(c=1)^C▒〖 i_a^c 〗   |  I_a^1= i_a^1,I_a^2= i_a^2  ,…,I_a^c= i_a^c )*p(〖ORGANIC〗_b=i_b^EA-∑_(c=1)^C▒〖 i_b^c 〗   ┤|  I_b^1= i_b^1,I_b^2= i_b^2  ,…,I_b^c= i_b^c).

So to get the parameters that maximize the posterior probability for our negative binomial distribution of installs we take,

argmax┬(i_(a,)^1 i_b^1,i_a^2,i_b^2…,i_a^c,i_b^c  )⁡〖L_log 〗+π_log

subject to 
                      〖organic〗_a+∑_(c=1)^C▒i_a^c = i_a^EA,
          〖organic〗_b+∑_(c=1)^C▒i_b^c = i_b^EA,
                     i_a^c+i_b^c=k^c   ,∀c .


PART III 
Quantum (Rounding) Error

While we transformed an intractable problem into a manageable one by discretizing our posterior in a form we can calculate more easily, it did not come without cost. Ideally, we would want the parameters of the Bernoulli probability distribution itself because that is the maximum likelihood estimator and ensures our distribution of installs will account for the uncertainty in our estimations more accurately than the discrete count solution. This is especially true when the SKAN install counts are low for a given conversion value but luckily the difference between the two distributions’ means converges to zero as sample size goes to infinity. For example, if there are 1000 installs in a given conversion value and the true probability is 25.7% an install occurred on day_a from that conversion value, then the maximum likelihood estimate (MLE) of our model will be 257 installs on day_a and there will be zero “quantizing” error. However, if there are only 10 installs for that conversion value in the SKAN postback, then the model’s MLE will be 3 installs on day_a which is off by 16.7% of the true MLE.

PART IV
Limitations and Future Iterations

The potentially significant limitation of the maximum posterior solution derived in Part II is that axioms v. and vii. do not hold in most situations. However, since problem complexity increases when incorporating the conditional dependence of organics on campaign installs (which is really endogeneity, both caused by UA impressions to some degree), and the codependence between organic and paid installs is believed to be marginal, this assumption allows for a tractable solution that is believed to not overly bias the results. Fortunately, there has been work done that measures the relationship between organics and UA spend which would serve as a useful starting point. Future iterations of this work should try to account for those co-dependencies between organic and paid installs in the likelihood term.

The biggest limitation of this approach is undoubtedly the large number of corner solutions that result from a large number of parameters and weak signal. In other words, there are many ways you can distribute the installs across the various campaigns that result in the same posterior density since the only distinguishing factor between them is the receipt times (priors). The more granular the send times the more effective this model will be in identifying unique optimal parameter values.

Lastly, another assumption that can be improved is that of the uniform random variable for time between chain break and postback send time. One possible improvement is to model the rate of install with a Poisson or time series trained on different panels. 



![image](https://user-images.githubusercontent.com/16949600/144113082-846189a9-13ab-43db-b137-f09fd5462ab8.png)

