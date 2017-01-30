
"""
Implements the likelihood function of Scherbaum et al. (2004)

Based on code from gmpe-smtk (Weatherill, 2014)
https://github.com/GEMScienceTools/gmpe-smtk
"""
def get_likelihood_values(ln_res):
	# res = log (ln) residuals normalised by the GMPE total sigma

	from numpy import array, fabs, isnan, sqrt, where, nan
	from scipy.special import erf
	from scipy.stats import scoreatpercentile
	
	# remove nans  
	tmpres = array(ln_res)
	nonnan = where(~isnan(tmpres))
	tmpres = tmpres[nonnan]
	
	if len(tmpres) != 0:
		zvals = fabs(tmpres)
		
		l_h_vals = 1.0 - erf(zvals / sqrt(2.))
		
		l_h_score = scoreatpercentile(l_h_vals, 50.0)
	
	else:
		l_h_vals = ln_res * nan
		l_h_score = nan
		
	return l_h_score, l_h_vals

"""
Implements of average sample log-likelihood estimator from
Scherbaum et al (2009)
"""
def get_loglikelihood_values(ln_res, periods):
    '''
    ln_res = nested array of residuals for each period in periods
    e.g. ln_res = vstack([ln_res['T1'], ln_res[T2], ln_res[Tn]])
    '''
    
    from numpy import array, log2, hstack, sum, nan, isnan, where
    from scipy.stats import norm
    
    log_residuals = array([])
    tllh = []
    per = []
    
    # enumerate across periods
    for i, period in enumerate(periods):
         if period >= 0.1 and period <= 2.0:
             # remove nans  
             tmpres = array(ln_res[:,i])
             nonnan = where(~isnan(tmpres))
             tmpres = tmpres[nonnan]
             
             if len(tmpres) != 0:
                 asll = log2(norm.pdf(tmpres, 0., 1.0))
                 
             else:
                 asll = ln_res[:,i] * nan
             
             log_residuals = hstack([log_residuals, asll])
             
             # check /per
             tllh.append(-(1. / float(len(asll))) * sum(asll))
             per.append(period)
    
    llh = -(1. / float(len(log_residuals))) * sum(log_residuals)
        
    return llh, tllh, per

"""
    Implements the Euclidean Distance-Based Ranking Method for GMPE selection
    by Kale & Akkar (2013)
    Kale, O., and Akkar, S. (2013) A New Procedure for Selecting and Ranking
    Ground Motion Predicion Equations (GMPEs): The Euclidean Distance-Based
    Ranking Method

Modified from Weatherill, 2014: https://github.com/GEMScienceTools/gmpe-smtk
"""
from numpy import fabs, mean, sum, ones, sqrt

def get_edr_values(obs, expected, stddev, periods, bandwidth=0.01, multiplier=3.0):
    """
    Calculates the EDR values for each GMPE
    :param float bandwidth:
        Discretisation width
    :param float multiplier:
        "Multiplier of standard deviation (equation 8 of Kale and Akkar)
    """
    edr_values = {}

    results = _get_edr(obs,
                       expected,
                       stddev,
                       bandwidth,
                       multiplier)
    edr_values["MDE Norm"] = results[0]
    edr_values["sqrt Kappa"] = results[1]
    edr_values["EDR"] = results[2]
    
    return edr_values

'''
def _get_gmpe_information(self, gmpe):
    """
    Extract the observed ground motions, expected and total standard
    deviation for the GMPE (aggregating over all IMS)
    """
    obs = array([], dtype=float)
    expected = array([], dtype=float)
    stddev = array([], dtype=float)
    for i, period in enumerate(periods):
        for context in self.contexts:
            obs = hstack([obs, log(context["Observations"][imtx])])
            expected = hstack([expected,
                                  context["Expected"][gmpe][imtx]["Mean"]])
            stddev = hstack([stddev,
                                context["Expected"][gmpe][imtx]["Total"]])
    return obs, expected, stddev
'''

def _get_edr(obs, expected, stddev, bandwidth=0.01, multiplier=3.0):
    """
    Calculated the Euclidean Distanced-Based Rank for a set of
    observed and expected values from a particular GMPE
    """
    nvals = len(obs)
    min_d = bandwidth / 2.
    kappa = _get_kappa(obs, expected)
    mu_d = obs - expected
    d1c = fabs(obs - (expected - (multiplier * stddev)))
    d2c = fabs(obs - (expected + (multiplier * stddev)))
    dc_max = ceil(max(array([max(d1c), max(d2c)])))
    num_d = len(arange(min_d, dc_max, bandwidth))
    mde = zeros(nvals)
    for iloc in range(0, num_d):
        d_val = (min_d + (float(iloc) * bandwidth)) * ones(nvals)
        d_1 = d_val - min_d
        d_2 = d_val + min_d
        p_1 = norm.cdf((d_1 - mu_d) / stddev) -\
            norm.cdf((-d_1 - mu_d) / stddev)
        p_2 = norm.cdf((d_2 - mu_d) / stddev) -\
            norm.cdf((-d_2 - mu_d) / stddev)
        mde += (p_2 - p_1) * d_val
    inv_n = 1.0 / float(nvals)
    mde_norm = sqrt(inv_n * sum(mde ** 2.))
    edr = sqrt(kappa * inv_n * sum(mde ** 2.))
    
    return mde_norm, sqrt(kappa), edr


def _get_kappa(obs, expected):
    """
    Returns the correction factor kappa
    """
    mu_a = mean(obs)
    mu_y = mean(expected)
    b_1 = sum((obs - mu_a) * (expected - mu_y)) /\
        sum((obs - mu_a) ** 2.)
    b_0 = mu_y - b_1 * mu_a
    y_c =  expected - ((b_0 + b_1 * obs) - obs)
    de_orig = sum((obs - expected) ** 2.)
    de_corr = sum((obs - y_c) ** 2.)
    
    return de_orig / de_corr
