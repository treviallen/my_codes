# -*- coding: utf-8 -*-
"""
Created on Tue Aug 05 16:16:38 2014

@author: tallen
"""

def boore_atkinson_site_getcoefs(T):
    from openquake.hazardlib.gsim.base import CoeffsTable
    from openquake.hazardlib import imt

    COEFFS = CoeffsTable(sa_damping=5, table="""\
            imt blin b1 b2
            pgv -0.600 -0.500 -0.06
            pga -0.360 -0.640 -0.14
            0.010 -0.360 -0.640 -0.14
            0.020 -0.340 -0.630 -0.12
            0.030 -0.330 -0.620 -0.11
            0.050 -0.290 -0.640 -0.11
            0.075 -0.230 -0.640 -0.11
            0.100 -0.250 -0.600 -0.13
            0.150 -0.280 -0.530 -0.18
            0.200 -0.310 -0.520 -0.19
            0.250 -0.390 -0.520 -0.16
            0.300 -0.440 -0.520 -0.14
            0.400 -0.500 -0.510 -0.10
            0.500 -0.600 -0.500 -0.06
            0.750 -0.690 -0.470 0.00
            1.000 -0.700 -0.440 0.00
            1.500 -0.720 -0.400 0.00
            2.000 -0.730 -0.380 0.00
            3.000 -0.740 -0.340 0.00
            4.000 -0.750 -0.310 0.00
            5.000 -0.750 -0.291 0.00
            7.500 -0.692 -0.247 0.00
            10.000 -0.650 -0.215 0.00  
            """)
            
    if T == 0.0:
        ct = COEFFS[imt.PGA()]
    elif T == -1.0:
        ct = COEFFS[imt.PGV()]
    else:
        ct = COEFFS[imt.SA(damping=5, period=T)]

    return ct

# returns amp factor for given input Vs30, period and PGA    
def boore_atkinson_siteamp(vs30, T, pga4nl):
    '''
    Reference PGA (pga4nl) is in g
    
    Site response factor is returned in non-log units
    '''
    from numpy import exp, log

    # get coefs for period T
    coeffs = boore_atkinson_site_getcoefs(T)
    
    # set constants
    vref = 760. # m/s
    a1 = 0.03 # g
    a2 = 0.09 # g
    pgalow = 0.06 # g
        
    # first get linear
    lnFlin = coeffs['blin'] * log(vs30 / vref)

    # get non-linear
    b1 = coeffs['b1']
    b2 = coeffs['b2']
    
    # set bnl
    bnl = 0
    if vs30 <= 180.:
        bnl = b1
    elif vs30 > 180. and vs30 <= 300.:
        bnl = (b1 - b2) * log(vs30 / 300.) / log(180. / 300.) + b2
    elif vs30 > 300. and vs30 < vref:
        bnl = b2 * log(vs30 / vref) / log(300. / vref)
            
    # find coefs c & d
    dx = log(a2 / a1)
    dy = bnl * log(a2 / pgalow)
    c = (3. * dy - bnl * dx) / dx**2
    d = -(2. * dy - bnl * dx) / dx**3
    
    # get lnFnl
    if pga4nl <= a1:
        lnFnl = bnl * log(pgalow / 0.1)
        
    elif pga4nl > a1 and pga4nl <= a2:
        lnFnl = bnl * log(pgalow / 0.1) + c * log(pga4nl / a1)**2 \
              + d * log(pga4nl / a1)**3
              
    elif pga4nl > a2:
        lnFnl = bnl * log(pga4nl / 0.1)
        
    # return total site response factor
    return exp(lnFlin + lnFnl), exp(lnFlin)

