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

    