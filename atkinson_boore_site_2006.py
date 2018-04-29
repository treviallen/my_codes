# -*- coding: utf-8 -*-
"""
Created on Tue Aug 05 16:16:38 2014

@author: tallen
"""

def atkinson_boore_site_getcoefs(T):
    from openquake.hazardlib.gsim.base import CoeffsTable
    from openquake.hazardlib import imt
    # coeficients from Atkinson & Boore etal (2006)
    
    # T = 0.3 value rounded from original

    COEFFS = CoeffsTable(sa_damping=5, table="""\
            imt blin b1 b2
            pgv -0.600 -0.495 -0.060
            pga -0.361 -0.641 -0.144
            0.025 -0.330 -0.624 -0.115
            0.031 -0.322 -0.618 -0.108
            0.04 -0.314 -0.609 -0.105
            0.05 -0.286 -0.643 -0.105
            0.063 -0.249 -0.642 -0.105
            0.079 -0.232 -0.637 -0.117
            0.1 -0.250 -0.595 -0.132
            0.125 -0.260 -0.560 -0.140
            0.159 -0.280 -0.528 -0.185
            0.2 -0.306 -0.521 -0.185
            0.25 -0.390 -0.518 -0.160
            0.3 -0.445 -0.513 -0.130
            0.4 -0.500 -0.508 -0.095
            0.5 -0.600 -0.495 -0.060
            0.625 -0.670 -0.480 -0.031
            0.769 -0.690 -0.465 -0.002
            1.0 -0.700 -0.440 0
            1.587 -0.726 -0.395 0
            2.0 -0.730 -0.375 0
            3.125 -0.740 -0.330 0
            4.0 -0.745 -0.310 0
            5.0 -0.752 -0.300 0
            """)
            
    if T == 0.0:
        ct = COEFFS[imt.PGA()]
    elif T == -1.0:
        ct = COEFFS[imt.PGV()]
    else:
        ct = COEFFS[imt.SA(damping=5, period=T)]

    return ct

# returns amp factor for given input Vs30, period and PGA
def atkinson_boore_siteamp(vs30, T, pgaBC):
    '''
    Reference PGA (pgaBC) is in g
    
    Site response factor is returned in non-log units
    '''
    from numpy import exp, log

    # limit amp factors to periods within AB06 range
    if T > 5.0:
        T = 5.0
        
    elif T < 0.025:
        T = 0.025
    
    # get coefs for period T
    coeffs = atkinson_boore_site_getcoefs(T)
    
    # set constants
    vref = 760. # m/s
    v1 = 180.
    v2 = 300.
    a1 = 0.03 # g
    a2 = 0.09 # g
    
    pgalow = 0.06 # g
        
    # first get linear
    lnFlin = coeffs['blin'] * log(vs30 / vref)

    # get non-linear
    b1 = coeffs['b1']
    b2 = coeffs['b2']
    
    # set bnl
    bnl = 0.
    if vs30 <= v1:
        bnl = b1
    elif vs30 > v1 and vs30 <= v2:
        bnl = (b1 - b2) * log(vs30 / v2) / log(v1 / v2) + b2
    elif vs30 > v2 and vs30 < vref:
        bnl = b2 * log(vs30 / vref) / log(v2 / vref)
            
    # convert pgaBC from g to cm/s**2
    pgaBC_cms2 = pgaBC * 981.
    
    # calculate site factor
    #print vs30, pgaBC_cms2
    if pgaBC_cms2 <= 60.:
        s = exp(lnFlin + bnl*log(60. / 100.))
    else:
        s = exp(lnFlin + bnl*log(pgaBC_cms2 / 100.))
       
    # return total site response factor
    return s

