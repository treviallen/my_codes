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

    