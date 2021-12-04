# POC SKAN_Revenue_Attribution
Cohort Install Attribution Model



This work is a quick/dirty prototype that provides a means to derive priors for the media mix model for campaign revenue attribution for SKAN.
Files cohort_attribution.RMD and cohort_attribution_v4.stan serve as a toy simulation for demonstration purposes. See the former's library calls for dependencies.


See SKAN_Revenue_Cohort_Attribution_V3.pdf for POC description and model derivations.

NEXT STEPS:
Setting the constraints for each campaign to the total installs in the postback for that campaign rather than the sum of all installs for all campaigns combined yields better performance in almost all cases, most notably with the log likelihood and smaller SD. The only issue is an engineering one as I had to hardcode the number of campaigns with separate campaign-install parameters  as I'm not sure how to set dynamic constraints for a vector - clearly this presents scaling issues.

In addition, while it was not done in this simulation, the priors must be derived from the bernoulli distributions obtained via the SKAN send times, and subsequently, the variance of the priors for each inst_c_a (each campaign's installs on day a) should necessarily reflect the uncertainty in the original Bernoulli distribution.
