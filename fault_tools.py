# -*- coding: utf-8 -*-
"""
Created on Tue May 21 16:10:47 2013

Functions for working with faults

    - mm2m: Converts units of mm to m
    - mag2disp: takes mw, ftype & scaling model
    - char_rec: Estimates characteristic return periods from slip rate, area, and M_char

@author: tallen
"""
def print_functions():
    txt = open('U:/Code/pycode/fault_tools.py').readlines()
    for line in txt:
        if line.find('def') >= 0:
            print line.strip('def ').strip('\n')


def beta2bval(beta):
    from numpy import log10, exp
    return log10(exp(beta))

def bval2beta(bval):
    from numpy import log
    return log(10**bval)

def mm2m(mm):
    return mm / 1000.

def m2mm(m):
    return m * 1000.

def sqkm2sqm(sqkm):
    return sqkm * 1.0E6

def mmpyr2mps(mmpyr):
    # first convert to m
    sr1 = mm2m(mmpyr)
    # get N secs per yr
    nsecs = 365 * 24 * 60 * 60
    return sr1 / nsecs

# assumes "sliprate" untis of mm/yr
def year_per_m(sliprate):
    mpy = mm2m(sliprate)
    return 1. / mpy

def set_kwarg_defaults(kwargs):
    model = 'WC94'
    ftype = 'all'

    # set kwarg variables
    for key in ('model', 'ftype'):
        if key in kwargs:
            # set fault type
            if key == 'ftype':
                ftype = kwargs[key]

            # set scaling relation
            if key == 'model':
                model = kwargs[key]

    return model, ftype

'''
Fault scaling tools
    mw: moment magnitude
    ftype: either 'ss', 'rs', 'ns', or 'all' - ftype='rs'
    model: scaling model name (e.g. WC94, L10) - model='WC94'
'''

'''do Wyss 1979'''
# calc mag from area in km**2
def area2mag_W79(area):
    from numpy import log10
    return log10(area) + 4.15

'''do Wells & Coppersmith 1994'''
# uses Average Displacement coefficients - units of m
def mag2avs_WC94(mw, ftype):
    if ftype == 'all':
        a = -4.80
        b = 0.69
        sig = 0.36
    elif ftype == 'ss':
        a = -6.32
        b = 0.90
        sig = 0.28
    elif ftype == 'rs':
        a = -0.74
        b = 0.08
        sig = 0.38
        print '\nNote: WC94 poorly constrained for reverse faulting events\n'
    elif ftype == 'ns':
        a = -4.45
        b = 0.63
        sig = 0.33
    return 10**(a + b * mw), sig # in m

def mag2area_WC94(mw, ftype):
    if ftype == 'all':
        a = -3.49
        b = 0.91
    elif ftype == 'ss':
        a = -3.42
        b = 0.90
    elif ftype == 'rs':
        a = -3.99
        b = 0.98
    elif ftype == 'ns':
        a = -2.87
        b = 0.82
    return 10**(a + b * mw) # in km**2

# assumes subsurface rupture
def mag2ruplen_WC94(mw, ftype):
    if ftype == 'all':
        a = -2.44
        b = 0.59
    elif ftype == 'ss':
        a = -2.57
        b = 0.62
    elif ftype == 'rs':
        a = -2.42
        b = 0.58
    elif ftype == 'ns':
        a = -1.88
        b = 0.50
    return 10**(a + b * mw) # in km

def mag2rupwid_WC94(mw, ftype):
    if ftype == 'all':
        a = -1.01
        b = 0.32
    elif ftype == 'ss':
        a = -0.76
        b = 0.27
    elif ftype == 'rs':
        a = -1.61
        b = 0.41
    elif ftype == 'ns':
        a = -1.14
        b = 0.35
    return 10**(a + b * mw) # in km

# MW to max slip
def mag2maxs_WC94(mw, ftype):
    if ftype == 'all':
        a = -5.46
        b = 0.82
        sig = 0.42
    elif ftype == 'ss':
        a = -7.03
        b = 1.03
        sig = 0.43
    elif ftype == 'rs':
        a = -1.84
        b = 0.29
        sig = 0.42
        print '\nNote: WC94 poorly constrained for reverse faulting events\n'
    elif ftype == 'ns':
        a = -5.90
        b = 0.89
        sig = 0.38
    return 10**(a + b * mw), sig # in m
    
def mag2disp_WC94(mw, ftype):
    if ftype == 'all':
        a = -4.80
        b = 0.69
        sig = 0.36
    elif ftype == 'ss':
        a = -6.32
        b = 0.90
        sig = 0.28
    elif ftype == 'rs':
        a = -0.74
        b = 0.08
        sig = 0.38
        print '\nNote: WC94 poorly constrained for reverse faulting events\n'
    elif ftype == 'ns':
        a = -4.45
        b = 0.63
        sig = 0.33
    return 10**(a + b * mw), sig # in m

# srl is in km
def srl2mag_WC94(srl, ftype):
    from numpy import log10
    if ftype == 'all':
        a = 5.08
        b = 1.16
    elif ftype == 'ss':
        a = 5.16
        b = 1.12
    elif ftype == 'rs':
        a = 5.0
        b = 1.22
    elif ftype == 'ns':
        a = 4.86
        b = 1.32
    return a + b * log10(srl) # in mw

# area is in km**2
def area2mag_WC94(area, ftype):
    from numpy import log10
    if ftype == 'all':
        a = 4.07
        b = 0.98
        sig = 0.24
    elif ftype == 'ss':
        a = 3.98
        b = 1.02
        sig = 0.23
    elif ftype == 'rs':
        a = 4.33
        b = 0.90
        sig = 0.25
    elif ftype == 'ns':
        a = 3.93
        b = 1.02
        sig = 0.25
    mw = a + b * log10(area)
    return mw, sig # in Mw units

'''do Hanks & Bakun, 2002 - Only for SS'''
def area2mag_HB02(area):
    print '\nNote: Hanks & Bakun, 2002 only appropriate for strike-slip events\n'
    from numpy import log10
    if area <= 537.:
        mw = log10(area) + 3.98
        sig = 0.03
    else:
        mw = (4./3.) * log10(area) + 3.07
        sig = 0.04

    return mw, sig

'''do Ellsworth-B, 2003 - Only for SS'''
def area2mag_E03(area):
    print '\nNote: Ellsworth-B, 2003 only appropriate for strike-slip events\n'
    from numpy import log10
    mw = log10(area) + 4.20
    sig = 0.10 # after OpenSHA: http://www.opensha.org/glossary-magScalingRelation

    return mw, sig

'''do Shaw 2009 - Only for SS'''
def area2mag_Sh09(area):
    print '\nNote: Shaw, 2009 only appropriate for strike-slip events\n'
    from numpy import log10, sqrt
    b = 6.9
    c = 3.98
    h = 15.6

    return log10(area) + (2./3.) * log10(max(1, sqrt(area / h**2))) \
           / ((1 + max(1, sqrt(area / (h**2 * b)))) / 2.0) + c

''' from Strasser etal 2010'''
''' area is in km**2'''
def area2mag_St10inter(area):
    from numpy import log10
    a   = 4.441
    b   = 0.846
    sig = 0.286
    return a + b * log10(area), sig # in mw

def area2mag_St10intra(area):
    from numpy import log10
    a   = 4.054
    b   = 0.981
    sig = 0.193
    return a + b * log10(area), sig # in mw

def mag2area_St10inter(mag):
    a   = -3.476
    b   = 0.952
    sig = 0.304
    return 10**(a + b * mag), sig # km**2

def mag2area_St10intra(mag):
    a   = -3.225
    b   = 0.890
    sig = 0.184
    return 10**(a + b * mag), sig # km**2

def srl2mag_St10inter(srl):
    from numpy import log10
    a = 4.868
    b = 1.392
    sig = 0.277
    return a + b * log10(srl), sig # in mw

def srl2mag_St10intra(srl):
    from numpy import log10
    a = 4.725
    b = 1.445
    sig = 0.234
#    print a, b
    return a + b * log10(srl), sig # in mw

def mag2srl_St10inter(mag):
    a = -2.447
    b = 0.585
    sig = 0.180
    return 10**(a + b * mag), sig # in km
    
def mag2srl_St10intra(mag):
    a = -2.350
    b = 0.562
    sig = 0.146
    return 10**(a + b * mag), sig # in km

def mag2wid_St10inter(mag):
    a =  -0.882
    b = 0.351
    sig = 0.173
    return 10**(a + b * mag), sig # in km
    
def mag2wid_St10intra(mag):
    a =  -1.058
    b = 0.356
    sig = 0.067
    return 10**(a + b * mag), sig # in km

'''do Blaser etal 2010'''
def mag2wid_Bl10rev(mag):
    a =  -1.86
    b = 0.46
    sig = 0.17
    return 10**(a + b * mag), sig # in km

def mag2len_Bl10rev(mag):
    a =  -2.37
    b = 0.57
    sig = 0.18
    return 10**(a + b * mag), sig # in km

'''do Somerville etal 2015 non-self-similar for interface'''
def mag2area_So15inter(mag):
    from numpy import log10
    from mag_tools import mw2m0
    
    m0 = mw2m0(mag) # in N-m
    a = 1.72E-09
    b = log10(4.17)
    sig = 1.481
    
    return 10**(log10(a) + b * log10(m0)), sig # in km^2
    
# do linear width
def mag2wid1_So15inter(mag):
    from numpy import log10
    from mag_tools import mw2m0
    
    m0 = mw2m0(mag)
    a = 6.75E-03
    b = log10(1.59)
    sig = 1.264
    
    return 10**(log10(a) + b * log10(m0)), sig # in km

# do bilinear width
def mag2wid2_So15inter(mag):
    from numpy import log10
    from mag_tools import mw2m0
    
    m0 = mw2m0(mag)
    a = 1.66E-04
    b = log10(1.90)
    sig = 1.259
    
    wid = 10**(log10(a) + b * log10(m0))
    
    # check if wid GT 200 km
    for i in range(0, len(wid)):
        if wid[i] > 200.:
            wid[i] = 200.
    
    return wid, sig # in km
    
# do average slip
def mag2avs_So15inter(mag):
    from numpy import log10
    from mag_tools import mw2m0
    
    m0 = mw2m0(mag)
    a = 3.39E-8
    b = 2.29
    sig = 1.522
    
    return 10**(log10(a) + b * log10(m0)), sig # in m
    
# do max slip
def mag2maxs_So15inter(mag):
    from numpy import log10
    from mag_tools import mw2m0
    
    m0 = mw2m0(mag)
    a = 2.70e-6
    b = 1.99
    sig = 1.501
    
    return 10**(log10(a) + b * log10(m0)), sig # in m
    
'''do Ichinose etal 2006 intraslab earthquakes'''
# do rupture area
def mag2area_I06intra(mag):
    from numpy import log10
    from mag_tools import mw2m0
    
    m0 = mw2m0(mag) # in N-m

    m0_dynecm = m0 * 10E7
    
    a = 0.57
    b = -13.5 # correct for dyne-cm to N-m
    sig = 0.9
    
    return 10**(a * log10(m0_dynecm) + b), sig # in km**2

'''do Somerville etal 2015 self-similar (SS) for interface'''
# do max slip
def mag2maxs_So15inter_SS(mag):
    from numpy import log10
    from mag_tools import mw2m0
    
    m0 = mw2m0(mag) # in N-m
    c = 5.00E-07
    sig = 1.508
    
    return 10**(log10(c) + (1./3.) * log10(m0)), sig # in m
    
# do average slip
def mag2avs_So15inter_SS(mag):
    from numpy import log10
    from mag_tools import mw2m0
    
    m0 = mw2m0(mag) # in N-m
    c = 1.23E-07
    sig = 1.527
    
    return 10**(log10(c) + (1./3.) * log10(m0)), sig # in m 
    
# do rupture area
def mag2area_So15inter_SS(mag):
    from numpy import log10
    from mag_tools import mw2m0
    
    m0 = mw2m0(mag) # in N-m
    c = 1.77E-10
    sig = 1.527
    
    return 10**(log10(c) + (2./3.) * log10(m0)), sig # in m     
    
'''do MUROTANI et al 2013 interface'''
def mag2area_Mu13inter(mag):
    from mag_tools import mw2m0
    
    m0 = mw2m0(mag)
    sig = 1.54
    return 1.34e-10 * m0**(2./3.), sig

# Average slip    
def mag2avs_Mu13inter(mag):
    from mag_tools import mw2m0
    
    m0 = mw2m0(mag)
    sig = 1.64
    return 1.66E-7 * m0**(1./3.), sig

'''do Leonard 2010 - active crust'''
# average slip from T6 coefs
def mag2avs_L10(mw, ftype):
    from numpy import log10
    # assume ss == all
    a = 1.67
    c = 0.833
    if ftype == 'rs':
        # first get subsurface L
        b = 4.24
        d = -1.30
    elif ftype == 'all' or ftype == 'ss':
        b = 4.17
        d = -1.34
    elif ftype == 'scr':
        b = 4.32
        d = -1.07
    flen = 10**((mw - b) / a) # in km
    return 10**(c * log10(flen) + d) # in m
    
def area2avs_L10(area, ftype):
    from numpy import log10
    # assume ss == all
    a = 0.5
    if ftype == 'ds':
        # first get subsurface L
        b = -4.42
    elif ftype == 'all' or ftype == 'ss':
        b = -4.43
    elif ftype == 'scr':
        b = -4.14
    return 10**(a * log10(area) + b) # in m

def area2mag_L10(area, ftype): # in km**2
    from numpy import log10
    a = 1.0
    if ftype == 'ds':
        b = 4.00
    elif ftype == 'ss':
        b = 3.99
    elif ftype == 'scr':
        b = 4.19
    return a * log10(area) + b # in MW
    
def len2mag_L10(length, ftype): # in km
    from numpy import log10
    a = 1.67
    if ftype == 'ds':
        b = 4.24
    elif ftype == 'ss':
        b = 4.17
    elif ftype == 'scr':
        b = 4.32
    return a * log10(length) + b # in MW
    
def mag2area_L10(mw, ftype): # in MW
    from numpy import log10
    a = 1.0
    if ftype == 'ds':
        b = 4.00
    elif ftype == 'ss':
        b = 3.99
    elif ftype == 'scr':
        b = 4.19
    return 10**((mw - b) / a) # in km**2
    
def mag2len_L10(mw, ftype): # in MW
    from numpy import log10
    a = 1.67
    if ftype == 'ds':
        b = 4.24
    elif ftype == 'ss':
        b = 4.17
    elif ftype == 'scr':
        b = 4.32
    return 10**((mw - b) / a) # in km

# uses mag2len and T5 coefs
def mag2wid_L10(mw, ftype): # in MW
    from numpy import log10
    flen = mag2len_L10(mw, ftype) * 1000. # in m
    
    a = 0.667
    if ftype == 'rs':
        b = 1.24
        print '\n only valid for RS if GT 5.5 km\n'
    elif ftype == 'ss':
        b = 1.18
        print '\n only valid for SS if 3.4-45 km\n'
    elif ftype == 'scr':
        b = 1.13
        print '\n only valid for SCR if GT 2.5 km\n'
    return 10**(a * log10(flen) + b) / 1000. # in km
    
def len2wid_L10(flen, ftype): # in MW
    from numpy import log10
    
    flen *= 1000 # convert from km to m
    a = 0.667
    if ftype == 'rs':
        b = 1.24
        print '\n only valid for RS if GT 5.5 km\n'
    elif ftype == 'ss':
        b = 1.18
        print '\n only valid for SS if 3.4-45 km\n'
    elif ftype == 'scr':
        b = 1.13
        print '\n only valid for SCR if GT 2.5 km\n'
    return 10**(a * log10(flen) + b) / 1000. # in km


'''do Leonard 2014 - SCR'''
def mag2len_L14(mw, ftype): # in MW - not sure if correct
    # assume mw is a list
    if ftype == 'scrrs':
        b = 0.953
        a = 5.12
        lsr = 10**((mw - a) / b)
        
        idx = lsr > 15.
        b = 1.667
        a = 4.32
        lsr[idx] = 10**((mw[idx] - a) / b)
    
    return lsr # in km
    
def len2mag_L14(lsr, ftype): # in km
    from numpy import log10
    from mag_tools import m02mw
    
    lsr *= 1000. # km2m
    
    if ftype == 'scrrs':
        if lsr <= 15.:
            a = 12.55
            b = 1.43
        else:
            a = 8.08
            b = 2.50
    
    logM0 = a + log10(lsr) * b # in m0
    
    return m02mw(10**logM0)

'''do Clark et al (2014) - non-extended crust'''
def srl2mag_Cea14(srl, reg):
    from numpy import log10
    from mag_tools import m02mw
    if reg == 'all':
        a = 16.23
        b = 2.03
        sig = 0.27
    elif reg == 'au': # Aust & Ungava
        a = 16.21
        b = 1.83
        sig = 0.18
        
    logM0 = a + b * log10(srl)
    return m02mw(10**logM0), sig # Mw

def mag2srl_Cea14(mag, reg):
    from numpy import log10
    from mag_tools import mw2m0
    if reg == 'all':
        a = 16.23
        b = 2.03
        sig = 0.27
    elif reg == 'au': # Aust & Ungava
        a = 16.21
        b = 1.83
        sig = 0.18
    logM0 = log10(mw2m0(mag))
    return 10**((logM0 - a) / b), sig # srl in km
    
''' do Johnston 1994 SCR '''
def mag2len_J94_SCR_lin(mw):
    return 10**((mw - 4.67) / 1.36) # in km
    
def len2mag_J94_SCR_quad(flen):
    from numpy import log10
    return 5.22 + 0.38 * log10(flen) + 0.37 * log10(flen)**2
    
''' do Allen & Hayes 2017 intraslab '''

# get mag 2 len in km**2
def mag2len_AH17_other(mw, rtype):
    b = 0.63
    if rtype == 'intra':
        a = -3.03
        sig = 0.14
    elif rtype == 'outer':
        a = -2.87
        sig = 0.08
    elif rtype == 'ss':
        a = -2.81
        sig = 0.15
        
    return 10**(a + b * mw), sig
    
# get mag 2 area in km**2
def mag2area_AH17_other(mw, rtype):
    b = 0.96
    if rtype == 'intra':
        a = -3.89
        sig = 0.19
    elif rtype == 'outer':
        a = -3.89
        sig = 0.11
    elif rtype == 'ss':
        a = -4.04
        sig = 0.20
        
    return 10**(a + b * mw), sig
    

''' STOP SCALING RELATIONS '''

# get displacement (in m) for earthquake of magnitude mw
def mag2disp(mw, **kwargs):
    # set default values
    model, ftype = set_kwarg_defaults(kwargs)

    # set models
    if model == 'WC94':
        disp = mag2disp_WC94(mw, ftype)
    elif model == 'L10AC':
        disp = mag2disp_L10(mw, ftype)

    return disp

# get subsurface rupture length (in km) for earthquake of magnitude mw
def mag2ruplen(mw, **kwargs):
    # set default values
    model, ftype = set_kwarg_defaults(kwargs)

    # set models
    if model == 'WC94':
        rup_len = mag2ruplen_WC94(mw, ftype)

    return rup_len

# get displacement (in m) for earthquake of magnitude mw
def srl2mag(srl, **kwargs):
    # set default values
    model, ftype = set_kwarg_defaults(kwargs)

    # set models
    if model == 'WC94':
        mw = srl2mag_WC94(srl, ftype)
    elif model == 'S10inter':
        mw, sig = srl2mag_St10inter(srl)
    elif model == 'S10intra':
        mw, sig = srl2mag_St10intra(srl)

    return mw

def return_area2mag_for_ss(area):
    from numpy import mean, std

    ftype = 'ss'
    # do Wells & Coppersmith 1994
    WC94mw, WC94sig = area2mag_WC94(area, ftype)
    # do Hanks & Bakun 2002
    HB02mw, HB02sig = area2mag_HB02(area)
    # do Ellsworth-B 2003
    E03mw, E03sig = area2mag_E03(area)
    # do Shaw 2009
    Sh09 = area2mag_Sh09(area)
    # do Leonard 2010
    L10mw = area2mag_L10(area, ftype)

    Mest = [WC94mw, HB02mw, E03mw, Sh09, L10mw]

    return Mest, mean(Mest), std(Mest)

def return_area2mag_for_rs(area):
    from numpy import mean, std

    ftype = 'rs'
    # do Wells & Coppersmith 1994
    WC94mw, WC94sig = area2mag_WC94(area, ftype)
    # do Leonard 2010
    L10mw = area2mag_L10(area, ftype)

    Mest = [WC94mw, L10mw]

    return Mest, mean(Mest), std(Mest)

'''
 estimate characteristic eq model return periods
    Assumes: "sliprate" untis of mm/yr
             "area" is in km**2
             "m_char" is characteristic earthquake
'''
def characteristic_mag_rec(m_char, sliprate, **kwargs):
    # set kwargs
    model, ftype = set_kwarg_defaults(kwargs)

    # get average displacement for m_char
    disp = mag2disp(m_char, model=model, ftype=ftype)
    print '\nAverage displacement for characteristic earthquake: ', \
          str("%0.1f" % disp),'m'

    # get years per characteristic earthquake
    srm = mm2m(sliprate) # in m / year
    avrec = disp / srm
    print 'Average recurrence: ',str("%0.0f" % avrec),'years\n'
    return avrec

'''
get G-R mag-recurrence from slip rate
    - sliprate in mm/yr
    - area in km**2
    - crust_thikness in km
    - mmax in Mw
needs kwargs area or crust_thikness
    - optional: model, ftype
'''
def bounded_mag_rec_from_slip(sliprate, mu, bval, mchar, mmin, **kwargs):
    from numpy import log10, nan, isnan

    '''
    Modified from: http://docs.openquake.org/oq-hazardlib/mfd.html
    '''

    # set kwarg variables
    area = nan
    usgs = False
    for key in ('crust_thikness', 'area', 'usgs'):
        if key in kwargs:
            # set fault type
            if key == 'crust_thikness':
                crust_thikness = kwargs[key]

            # set fault area
            if key == 'area':
                area = kwargs[key]

            # set boolean USGS model
            if key == 'usgs':
                usgs = kwargs[key]

    model, ftype = set_kwarg_defaults(kwargs)

    try:
        # get mmax
        if usgs == True:
            mmax = mchar
        else:
            mmax = mchar + 0.25 # following YC85

        if isnan(area):
            rup_len = mag2ruplen(mmax, model=model, ftype=ftype)
            area = crust_thikness * rup_len

        # set constants
#        phi = 1.27 # from Hyndman and Weichert [1983]
        c = 1.5 # Hanks & Kanamori - logM0 = 1.5 Mw + 9.05 (M0 in N-m)
        d = 9.05

        # get moment rate in N-m / yr
        M0r = mu * sqkm2sqm(area) * mm2m(sliprate)

        cmb = c - bval

        # get A0, where A0 = log(N0)
        A0 = log10(M0r * cmb) - log10(10**(cmb * mmax) - 10**(cmb * mmin)) \
             - log10(bval) - d

        return A0, mmax

    except:
        print "\nRequires kwargs 'area' or 'crust_thikness\n"

def characteristic_mag_rec_from_slip(sliprate, mu, bval, mmin, mchar, **kwargs):
    from numpy import log10, nan, isnan, exp
    '''
    modified from OpenQuake:
    http://docs.openquake.org/oq-hazardlib/_modules/openquake/hazardlib/mfd/youngs_coppersmith_1985.html#YoungsCoppersmith1985MFD
    '''
    # set kwarg variables
    area = nan
    usgs = False
    for key in ('crust_thikness', 'area', 'usgs'):
        if key in kwargs:
            # set fault type
            if key == 'crust_thikness':
                crust_thikness = kwargs[key]

            # set fault area
            if key == 'area':
                area = kwargs[key]

            # set boolean USGS model
            if key == 'usgs':
                usgs = kwargs[key]

    model, ftype = set_kwarg_defaults(kwargs)

    try:
        # boxcar function width
        DELTA_CHAR = 0.5
        if usgs == True:
            mmax = mchar + 0.1
        else:
            mmax = mchar + DELTA_CHAR / 2 # following YC85

        if isnan(area):
            rup_len = mag2ruplen(mchar, model=model, ftype=ftype)
            area = crust_thikness * rup_len

        beta = bval2beta(bval)

        # seismic moment (in Nm) for the maximum magnitude
        c = 1.5
        d = 9.05
        M0x = 10 ** (c * mmax + d)

        # get moment rate in N-m / yr
        M0r = mu * sqkm2sqm(area) * mm2m(sliprate)

        # equations (16) and (17) solved for N(min_mag) and N(char_mag)
        c1 = exp(-beta * (mmax - mmin - 0.5))
        c2 = exp(-beta * (mmax - mmin - 1.5))
        c3 = beta * c2 / (2 * (1 - c1) + beta * c2)
        c4 = (bval * (10 ** (-c / 2)) / (c - bval)) + \
             (bval * exp(beta) * (1 - (10 ** (-c / 2))) / c)
        n_min_mag = (1 - c1) * M0r / ((1 - c3) * c1 * M0x * c4)
        n_char_mag = c3 * n_min_mag

        A0 = log10((n_min_mag - n_char_mag) \
                / (10 ** (-bval * mmin) - 10 ** (-bval * (mchar - 0.25))))

        return A0, mmax, n_min_mag, n_char_mag

    except:
        print "\nRequires kwargs 'area' or 'crust_thikness\n"

#def truncated_mag_rec(M0r, cmb, mmax, bval):
#    return log10(M0r * cmb) - (cmb * mmax) - log10(bval) - d

def solve_A0(Ai, mmin, mmax, bval, M0r):
    from numpy import arange
#    print Ai, mmin, mmax, bval, M0r
    c = 1.5 # Hanks & Kanamori - logM0 = 1.5 Mw + 9.05 (M0 in N-m)
    d = 9.05
    m0rate = []
    N = []
    bin_width = 0.01
    mrange = arange(mmin, mmax, bin_width)
    for mw in mrange:
        mag_lo = mw - bin_width / 2.0
        mag_hi = mw + bin_width / 2.0
        tmpN = 10**(Ai - bval*mag_lo) - 10**(Ai - bval*mag_hi)
        N.append(tmpN[0])
        m0 = 10**(c*mw + d)
        m0rate.append(tmpN[0] * m0)
#    print N

    return abs(M0r - sum(m0rate))

def truncated_mag_rec_from_slip(sliprate, mu, bval, mmin, mchar, **kwargs):
    from numpy import nan, isnan
    from scipy.optimize import fmin

    '''
#    from Andersom & Luco 1983:
    '''
    # set kwarg variables
    area = nan
    for key in ('crust_thikness', 'area'):
        if key in kwargs:
            # set fault type
            if key == 'crust_thikness':
                crust_thikness = kwargs[key]

            # set fault area
            if key == 'area':
                area = kwargs[key]

    model, ftype = set_kwarg_defaults(kwargs)

    try:
        mmax = mchar

        if isnan(area):
            rup_len = mag2ruplen(mmax, model=model, ftype=ftype)
            area = crust_thikness * rup_len

        # set constants
        c = 1.5 # Hanks & Kanamori - logM0 = 1.5 Mw + 9.05 (M0 in N-m)
        d = 9.05
        M0x = 10 ** (c * mmax + d)

        # get moment rate in N-m / yr
        M0r = mu * sqkm2sqm(area) * mm2m(sliprate)

        # solve for A0
        Ai = 3.0 # initial guess
        A0 = fmin(solve_A0, Ai, args=(mmin, mmax, bval, M0r))

        return A0[0], mchar # return mchar as mmax
    except:
        print "\nRequires kwargs 'area' or 'crust_thikness\n"

def get_cummulative_stats(A0, bval, mc, curvetype, mrange, n_char_mag, name, **kwargs):
    from numpy import round
#    from fault_tools import mag2disp_WC94, mag2disp_L10

    model, ftype = set_kwarg_defaults(kwargs)

    # get cummulative eq num
    mrate = []
    bin_width = round((mrange[1] - mrange[0]), decimals=3)
    mmin = mrange[0]
    if curvetype == 'truncated':
        mx = mc
    else:
        mx = mc + 0.25

    for mw in mrange:
        tmpN = 0.0
        mag_lo = mw - bin_width / 2.0
        mag_hi = mw + bin_width / 2.0

        # assume exponential distribution for bounded & truncated model
        if curvetype == 'bounded' or curvetype == 'truncated':
            if mw <= mx:
                tmpN = 10**(A0 - bval*mag_lo) - 10**(A0 - bval*mag_hi)

        # get rate according to exponential distribution for characteristic model
        elif curvetype == 'characteristic':

            if mw < mc-0.25:
                tmpN = 10**(A0 - bval*mag_lo) - 10**(A0 - bval*mag_hi)

            # else get characteristic rate
            elif mw <= mx:
                tmpN = (n_char_mag / 0.5) * bin_width

        mrate.append(tmpN)

    # sum slip and cummulative rate
    slip_per_mag = []
    cumrate = []
    k = 0
    for mw in mrange:
        cumrate.append(sum(mrate[k:]))
        if  mw >= mmin and mw <= mx:
            if model == 'WC94':
                slip_per_mag.append(mag2disp_WC94(mw, 'ss') * 1000 * mrate[k])
            elif model == 'L10':
                slip_per_mag.append(mag2disp_L10(mw, 'ss') * 1000 * mrate[k])
        k += 1

    #print name, curvetype, str(sum(slip_per_mag)), 'mm/yr'

    return mrate, slip_per_mag, cumrate

def get_betaplt_frisk(beta, N0, mmin, mmax):
    from numpy import arange, exp

    mrange = arange(mmin, mmax-0.01, 0.01)

    betacurve = N0 * exp(-beta  *mrange) * (1 - exp(-beta * (mmax - mrange))) \
                 / (1 - exp(-beta * mmax))

    return betacurve, mrange

def get_oq_incrementalMFD(beta, N0, mmin, mmax, binwid):
    from numpy import arange, exp

    mrange = arange(mmin+0.5*binwid, mmax, binwid)

    betacurve = N0 * exp(-beta  *mrange) * (1 - exp(-beta * (mmax - mrange))) \
                / (1 - exp(-beta * mmax))

    return betacurve, mrange


def sliprate_from_GR_params(aval, bval, mmin, mmax, **kwargs):
    '''
    where:
        aval = log N0
        bval = b-value
    kwargs
        curvetype = distribution (i.e. bounded or characteristic)
        model = fault scaling relation (gets area from mmax)
        ftype = fault type (all, ss, rs, or ns)
        mu = crustal rigidity (in Pa)
        area = fault area in km**2
    '''
    from numpy import nan, isnan, log10, exp

    # first set default kwargs
    model = 'WC94'
    ftype = 'all'
    curvetype = 'bounded'
    mu = 3E10 # Pa
    area = nan

    kwarg_keys = ('model', 'ftype', 'curvetype', 'mu', 'area')

    # set kwarg variables
    for key in kwarg_keys:
        if key in kwargs:
            # set fault type
            if key == 'ftype':
                ftype = kwargs[key]

            # set scaling relation
            if key == 'model':
                model = kwargs[key]

            # set scaling relation
            if key == 'curvetype':
                curvetype = kwargs[key]

            # set scaling relation
            if key == 'mu':
                mu = kwargs[key]

            # set scaling relation
            if key == 'area':
                area = kwargs[key]

    # if not set, get fault area
    mchar = mmax - 0.25
    if isnan(area):
        if model == 'WC94':
            area = mag2area_WC94(mmax, ftype) # in km**2
    area = sqkm2sqm(area)

    # next get moment rate
    c = 1.5
    d = 9.05

    if curvetype == 'bounded':
        cmb = c - bval

        # original code - replaced below
#        ai = aval + log10(bval) + d
#        M0r = ((10**ai) / cmb) * (10 ** (cmb * mmax) - 10 ** (cmb * mmin))

        # alt code for negative beta!
        ai10 = 10**(aval + d) * bval
        M0r = (ai10 / cmb) * (10 ** (cmb * mmax) - 10 ** (cmb * mmin))

    # note - this is not quite right!
    if curvetype == 'characteristic':
        beta = bval2beta(bval)
        M0x = 10 **(c * mmax + d)
        c1 = exp(-beta * (mmax - mmin - 0.5))
        c2 = exp(-beta * (mmax - mmin - 1.5))
        c3 = beta * c2 / (2 * (1 - c1) + beta * c2)
        c4 = (bval * (10 ** (-c / 2)) / (c - bval)) + \
             (bval * exp(beta) * (1 - (10 ** (-c / 2))) / c)

        # get n_min_mag and n_char_mag
        if mmin <= mchar - 0.25:
            n_min_mag = 10**(aval - bval*mmin)

            M0r = n_min_mag * ((1 - c3) * c1 * M0x * c4) / (1 - c1)

    # now get slip rate from moment
    sliprate = m2mm(M0r / (mu * area))

    return M0r, sliprate

# get cummulative slip from negative beta
def cumslip_from_rate(mrange, cumrate, mu, area):
    from mag_tools import mw2m0
    from numpy import nan, nansum, array

    area = sqkm2sqm(area)
    M0array = []
    for i, mw in enumerate(mrange[0:-1]):
        # get eq rate
        rate = cumrate[i] - cumrate[i+1]
        if rate < 0.0 or cumrate[i+1] < 0.0:
            tmpM0r = nan
        else:
            tmpM0r = mw2m0(mw) * rate

        M0array.append(tmpM0r)

    # now sum moment rate
    M0r = nansum(M0array)

    # clac total slip
    sr = m2mm(M0r / (mu * area))

    # get slip per mag bin
    srpb = m2mm(array(M0array) / (mu * area))

    # get cummulative slip
    cumsr = []
    for i in range(0, len(srpb)):
        cumsr.append(nansum(srpb[i:]))

    return M0r, sr, srpb, cumsr

def fault_normal_from_shp(shpfile, fault, dip, ztor):
    from numpy import tan, radians
    import shapefile
    from mapping_tools import get_field_data, get_line_parallels

    # read shapefile with fault length and sliprate
    sf = shapefile.Reader(shpfile)
    fcode = get_field_data(sf, 'CODE', 'str')
    shapes = sf.shapes()

    # get down-dip range
    rng = ztor / tan(radians(dip)) # surface proj to top of rupt

    for i in range(0,len(fcode)):
        if fcode[i] == fault:
            pts = shapes[i].points
            posazpts, negazpts = get_line_parallels(pts, rng)

    return posazpts, negazpts

def get_downdip_depth(lat1, lon1, lat2, lon2, dip):
    from mapping_tools import distance
    from numpy import sin, radians

    rupwid = distance(lat1, lon1, lat2, lon2)[0]

    return rupwid * sin(radians(dip))
























