# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 14:17:35 2013

@author: tallen
"""

def print_functions():
    txt = open('U:/Code/pycode/mag_tools.py').readlines()
    for line in txt:
        if line.find('def') >= 0:
            print line.strip('def ').strip('\n')

# calculate seismic moment (in N-m) from mw
def mw2m0(mw):
    return 10**(3. * mw / 2.  + 9.05)

# calculate mw from seismic moment (in N-m)
def m02mw(m0):   
    from numpy import log10
    return 2 * log10(m0) / 3 - 6.03
    
    

'''
Begin magnitude conversions here
'''

def nsha18_ml2mw(ml):
    # use ODR polynomial of simulated data from Ghasemi & Allen (2017)
    a = 0.04160769
    b = 0.48058286
    c = 1.39485216
    
    # get Mw
    return a*ml**2 + b*ml + c
    
def nsha18_bilin_ml2mw(ml):
    from numpy import zeros_like
    
    # add bi-linear coefs - not used in NSHA18
    a_bl = 0.66053496
    b_bl = 1.20883045
    c_bl = 0.98659071
    hx_bl = 4.25
    hy_bl =  a_bl * hx_bl + b_bl
    
    if ml <= hx_bl:
        mw_bl = a_bl*ml + b_bl
    else:
        mw_bl = c_bl*(ml - hx_bl) + hy_bl
        
    return mw_bl

def nsha18_bilin_mw2ml(mw):
    from numpy import zeros_like
    
    # add bi-linear coefs - not used in NSHA18
    a_bl = 0.66053496
    b_bl = 1.20883045
    c_bl = 0.98659071
    hx_bl = 4.25
    hy_bl =  a_bl * hx_bl + b_bl
    
    if mw <= hy_bl:
        ml_bl = (mw - b_bl) / a_bl
    else:
        ml_bl = hx_bl + (mw - hy_bl) / c_bl
        
    return ml_bl


# from Mueller (pers comm) - from Sipkin SRL 2003
def sipkin_mb2mw(mb):
    if mb > 5.3:
        mw = 1.46 * mb - 2.42
    else:
        mw = mb
    return mw
    
# from Mueller (pers comm) - eyball fits to curves in Utsu IASPEI 2002
def utsu_ms2mw(ms):
    if ms < 5.8:
        mw = 0.75 * (ms + 1.93)
    elif ms > 7.8:
        mw = 1.50 * (ms - 2.60)
    else:
        mw = ms
    return mw
    
def utsu_ml2mw(ml):
    if ml > 6.5:
        mw = 1.67 * (ml - 2.60)
    else:
        mw = ml
    return mw
    
# do Ghasemi 2017 fixed-hinge bilinear GOR ML2MW as used for the NSHA18
def ghasemi_bl_ml2mw(ml):
    a1 = 0.66199378
    a2 = 1.2156352
    a3 = 1.07488336
    mx = 4.5
    my = a1 * mx + a2
    
    if ml <= mx:
        mw = a1 * ml + a2
    else:
        mw = a3 * (ml - mx) + my
        
    return mw

# this is what was used in the DRAFT NSHA!!!    
def ghasemi_bl_ml2mw_stupid_trevor(ml):
    a1 = 0.66199378
    a2 = 1.2156352
    a3 = 1.2156352
    mx = 4.5
    my = a1 * mx + a2
    
    if ml <= mx:
        mw = a1 * ml + a2
    else:
        mw = a3 * (ml - mx) + my
        
    return mw

def ghasemi_bl_mw2ml(mw):
    a1 = 0.66199378
    a2 = 1.2156352
    a3 = 1.07488336
    mx = 4.5
    my = a1 * mx + a2
    
    ml = (mw - a2) / a1
    
    if ml > mx:
        ml = mx + (mw - my) / a3
        
    return ml
    
def ghasemi_bl_mb2mw(mb):
    c1 = 1.1438907424442797
    c2 = 1.0509630268292971
    
    return c1 * mb + c2
    

def ghasemi_bl_ms2mw(ms):    
    c1 = 0.84896727404297323
    c2 = -0.87192285009579173
    
    return c1 * mb + c2

    