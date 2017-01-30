# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 15:24:21 2016

@author: tallen
"""
def beta2bval(beta):
    from numpy import log10, exp
    return log10(exp(beta))

def bval2beta(bval):
    from numpy import log
    return log(10**bval)

# code simplified from OpenQuake HMTK
def aki_maximum_likelihood(mrng, number_obs, mc):
    '''
    mrng = range of mags
    number_obs = observations per mag bin
    mc = magnitude of completenss
    '''    
    from numpy import array, diff, exp, log, log10, sqrt, sum, min, nan, where
    
    # get indexes of mags >= Mc
    mvals = array(mrng)
    number_obs = array(number_obs)
    idx = where(mvals >= mc)[0]
    dmag = diff(mrng)[0]
    
    mvals = mvals[idx]
    number_obs = number_obs[idx]
    
    # Get Number of events, minimum magnitude and mean magnitude
    neq = sum(number_obs)
    if neq <= 1:
        # Cannot determine b-value (too few event) return NaNs
        print 'Too few events (<= 1) to calculate b-value'
        return nan, nan
    
    else:
        m_min = min(mvals)
        m_ave = sum(mvals * number_obs) / neq
        
        # Calculate b-value
        bval = log10(exp(1.0)) / (m_ave - m_min + (dmag / 2.))
        
        # Calculate sigma b from Bender estimator
        sigma_b = sum(number_obs * ((mvals - m_ave) ** 2.0)) / (neq * (neq - 1))
        sigma_b = log(10.) * (bval ** 2.0) * sqrt(sigma_b)
        return bval, sigma_b
        
def fit_a_value(cum_mags, b_val, mrng, mc, m_upper):
    '''
    cum_mags = cummulative number of events per mag bin
    mrng = magnitude range
    mc = magnitude of completenss
    '''
    from numpy import array, log10, mean, where
    
    mrng = array(mrng)
    cum_mags = array(cum_mags)
    
    idx = where((mrng >= mc) & (mrng <= m_upper))[0]
    
    a_val = mean(log10(cum_mags[idx]) + b_val*mrng[idx])
    '''
    data = Data(array(mvals[idx]),log10(array(cum_mags[idx])))
    
    model = odrpack.Model(fit_intercept)
    odr = odrpack.ODR(data, model, beta0=[2.])
    odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=odr
    out = odr.run()
    out.pprint()
    c = out.beta[0]
    '''
        
    return a_val
    
# adapted from OQ-HMTK
def weichert_algorithm(tper, fmag, nobs, mrate=0.0, bval=1.0,
                       itstab=1E-5, maxiter=1000):
    import warnings
    """
    Weichert algorithm

    :param tper: length of observation period corresponding to magnitude
    :type tper: numpy.ndarray (float)
    :param fmag: central magnitude
    :type fmag: numpy.ndarray (float)
    :param nobs: number of events in magnitude increment
    :type nobs: numpy.ndarray (int)
    :keyword mrate: reference magnitude
    :type mrate: float
    :keyword bval: initial value for b-value
    :type beta: float
    :keyword itstab: stabilisation tolerance
    :type itstab: float
    :keyword maxiter: Maximum number of iterations
    :type maxiter: Int
    :returns: b-value, sigma_b, a-value, sigma_a
    :rtype: float
    """
    import numpy as np
    
    beta = bval * np.log(10.)
    d_m = fmag[1] - fmag[0]
    itbreak = 0
    snm = np.sum(nobs * fmag)
    nkount = np.sum(nobs)
    iteration = 1
    while (itbreak != 1):
        beta_exp = np.exp(-beta * fmag)
        tjexp = tper * beta_exp
        tmexp = tjexp * fmag
        sumexp = np.sum(beta_exp)
        stmex = np.sum(tmexp)
        sumtex = np.sum(tjexp)
        stm2x = np.sum(fmag * tmexp)
        dldb = stmex / sumtex
        if np.isnan(stmex) or np.isnan(sumtex):
            warnings.warn('NaN occurs in Weichert iteration')
            return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
            #raise ValueError('NaN occers in Weichert iteration')

        d2ldb2 = nkount * ((dldb ** 2.0) - (stm2x / sumtex))
        dldb = (dldb * nkount) - snm
        betl = np.copy(beta)
        beta = beta - (dldb / d2ldb2)
        sigbeta = np.sqrt(-1. / d2ldb2)

        if np.abs(beta - betl) <= itstab:
            # Iteration has reached convergence
            bval = beta / np.log(10.0)
            sigb = sigbeta / np.log(10.)
            fngtm0 = nkount * (sumexp / sumtex)
            fn0 = fngtm0 * np.exp((beta) * (fmag[0] - (d_m / 2.0)))
            stdfn0 = fn0 / np.sqrt(nkount)
            a_m = fngtm0 * np.exp((-beta) * (mrate -
                                            (fmag[0] - (d_m / 2.0))))
            siga_m = a_m / np.sqrt(nkount)
            itbreak = 1
        else:
            iteration += 1
            if iteration > maxiter:
                warnings.warn('Maximum Number of Iterations reached')
                return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    return bval, sigb, a_m, siga_m, fn0, stdfn0
