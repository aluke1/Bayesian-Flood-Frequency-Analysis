bayes_LPIII.m fits the stationary (ST) and non-stationary (NS) LPIII distributions to peak streamflow data using Bayesian inference.  
Post-processing following parameter estimation calculates and plots selected return periods vs return levels, 
input time series vs. the mean of distribution, and the density of return level estimates for a selected return period. 

Required input data is an annual maximum discharge record.  Note:  This program is not intended to inform design 
or insurance products.  Please refer to federal and state guidelines for such procedures. For methodological details 
please see sections 4.1.1 - 4.1.3 of "Predicting nonstationary flood frequencies: 
 Evidence supports an updated stationarity thesis in the United States". Luke et al 2017

In the NS model, the mean of the LPII distribution changes as a function of time.  Return periods are not calculated in a true 
"non-stationary" sense, but rather for a fixed value of the mean. The updated Stationaty return periods are obtained by calculating the 
return periods associated with the NS mean at the end of the fitting period, or t = t(end).  The updated ST return periods 
are calculated and displayed by default when applying bayes_LPIII.m for estimation of the NS LPIII model parameters. 

Required function files (must be in same directory) 
bayes_LPIII.m -> main program 
dream_zs.m    -> MCMC algorithm 
lp3inv.m      -> inverse of LPII CDF based on Wilson - Hilferty transformation 
lp3pdf.m      -> PIII probability density function calculated in log-space (returns actual density) 
prior_pdf.m   -> computes prior density of parameter combination 
prior_rnd.m   -> random draw from prior (for initialization) 
post_pdf.m    -> computes unnormalized posterior density for at proposal theta and data X.  Likelihood function in this file
record.txt    -> example record (United States Geological Survey site number 08074500, available at: URL:http://nwis.waterdata.usgs.gov/nwis/peak?) 

Please see comments in program files and "Predicting nonstationary flood frequencies: 
 Evidence supports an updated stationarity thesis in the United States" for further details.  

