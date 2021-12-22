# POC SKAN_Revenue_Attribution
Cohort Install Attribution Model



This work is a quick/dirty prototype that provides a means to derive priors for the media mix model for campaign revenue attribution for SKAN.
Files cohort_attribution.RMD and cohort_attribution_v4.stan serve as a toy simulation for demonstration purposes. See the former's library calls for dependencies.


See SKAN_Revenue_Cohort_Attribution_V3.pdf for POC description and model derivations.

NEXT STEPS:
Setting the constraints for each campaign to the total installs in the postback for that campaign rather than the sum of all installs for all campaigns combined yields better performance in almost all cases, most notably with the log likelihood and smaller SD. The only issue is an engineering one as I had to hardcode the number of campaigns with separate campaign-install parameters  as I'm not sure how to set dynamic constraints for a vector - clearly this presents scaling issues.

The assumptions around independence of the organics (both from previous days' organics and the spend on previous days) is known to be violated to some extent, the latter more severely than the former, but as a simplifying assumption it helps make the math simpler in these early iterations. Future versions should at least include the organics on day "a" in the conditional probability of organics on day "b". In addition, the probability of organics given the sum of the paid installs on that day should be updated to account for the halo effect estimated in previous analyses, which is a function of the (very) short-term time series of total marketing spend.

Lastly, the Poisson Binomial (PBD) PMF is not working correctly as it results in very odd shaped posteriors with massive variance and posterior estimates worse than the priors. Currently, I am using a normal with the mean and variance calculated from the Poisson Binomial distributions to model the priors for each campaign and the likelihood, which seems to work relatively well as an appoximation of the PBD.
