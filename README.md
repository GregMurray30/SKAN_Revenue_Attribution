# POC SKAN_Revenue_Attribution
Cohort Install Attribution Model



This work is a "quick and dirty" prototype that provides a means to derive priors for the media mix model for campaign revenue attribution for SKAN.
Files cohort_attribution.RMD and cohort_attribution_v2.stan serve as a toy simulation for demonstration purposes. See the former's library calls for dependencies.


See SKAN_Revenue_Cohort_Attribution_V2.docx for POC description and model derivations.

Next steps:
Setting the constraints for each campaign to the total installs in the postback for that campaign rather than the sum of all installs for all campaigns combined yields better performance in almost all cases, most notably with the log likelihood and smaller SD. The only issue is an engineering one as I had to hardcode the number of campaigns with separate campaign-install parameters  as I'm not sure how to set dynamic constraints for a vector - clearly this presents scaling issues.
