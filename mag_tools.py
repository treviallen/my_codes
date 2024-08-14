# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 14:17:35 2013

@author: tallen
"""

def print_functions():
    txt = open('U:/Code/pycode/mag_tools.py').readlines()
    for line in txt:
        if line.find('def') >= 0:
            print(line.strip('def ').strip('\n'))

# calculate seismic moment (in N-m) from mw
def mw2m0(mw):
    return 10**(3. * mw / 2.  + 9.05)

# calculate mw from seismic moment (in N-m)
def m02mw(m0):   
    from numpy import log10
    return 2 * log10(m0) / 3 - 6.03
    
# calculate mw from seismic moment (in dyne-cm)
def m02mw_dyne_cm(m0):   
    from numpy import log10
    return 2 * log10(m0) / 3 - 10.7    
    

'''
Begin magnitude conversions here
'''

# ML 2 MW conversion used in the NSHA18
def nsha18_ml2mw(ml):
    # use ODR polynomial of simulated data from Ghasemi & Allen (2017) and Allen et al (2018)
    a = 0.04160769
    b = 0.48058286
    c = 1.39485216
    sig = 0.195
    
    # get Mw
    return a*ml**2 + b*ml + c
    
def solve_nsha18_mw2ml(mw_list):
    '''
    mw_list = list of mws to convert to ml - numpy array
    '''
    from numpy import arange, interp
    from mag_tools import nsha18_ml2mw
    
    # make finely discrestised ml array
    ml_array = arange(0,7.8, 0.001)
    
    # calculate equivalent mws
    mw_array = nsha18_ml2mw(ml_array)

    # interpolate over mw_list
    return interp(mw_list, mw_array, ml_array, left=0, right=7.8)
    
def nsha18_fixed_bilin_ml2mw(ml):
    from numpy import zeros_like
    
    # add bi-linear coefs - not used in NSHA18
    a_bl = 0.66053496
    b_bl = 1.20883045
    c_bl = 0.98659071
    hx_bl = 4.25
    sig = 0.144
    
    hy_bl =  a_bl * hx_bl + b_bl
    
    if ml <= hx_bl:
        mw_bl = a_bl*ml + b_bl
    else:
        mw_bl = c_bl*(ml - hx_bl) + hy_bl
        
    return mw_bl, sig

def nsha18_fixed_bilin_mw2ml(mw):
    from numpy import zeros_like
    
    # add bi-linear coefs - not used in NSHA18
    a_bl = 0.66053496
    b_bl = 1.20883045
    c_bl = 0.98659071
    hx_bl = 4.25
    sig = 0.144
    
    hy_bl =  a_bl * hx_bl + b_bl
    
    if mw <= hy_bl:
        ml_bl = (mw - b_bl) / a_bl
    else:
        ml_bl = hx_bl + (mw - hy_bl) / c_bl
        
    return ml_bl, sig
    
def nsha18_mb2mw(mb):
    return 1.20 * mb - 1.176
    
def nsha18_ms2mw(ms):
    return 0.075 * ms**2 - 3.357
    
def nsha23_mb2mw(mb):
    
    mhx = 5.642
    if mb <= mhx:
        mw = 1.065 * mb - 0.6433
    else:
        mhy = 1.065 * mhx - 0.6433
        mw = 1.973 * (mb - mhx) + mhy
    
    return mw
    
def nsha23_ms2mw(ms):
    return 0.0723 * ms**2 - 3.482

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
    
    return c1 * ms + c2

# Ghasemi et al (2020) ISC mb2MW conversion for PNG region
def ghasemi20_mb2mw_png(mb):
    return 1.085 * mb - 0.207
    
##########################################################
# ISC-GEM conversions
##########################################################

def digiacomo15_mb2mw_linear(mb):
    return 1.38 * mb - 1.79
    	
def digiacomo15_mb2mw_exp(mb):
    from numpy import exp
    return exp(-4.664 + 0.859 * mb) + 4.555

def get_au_ml_zone(eqlos, eqlas):
    '''
    eqlos & eqlas are a list of lats/lons
    
    returns list of regions
    '''
    
    import shapefile
    from shapely.geometry import Point, Polygon
    from mapping_tools import get_field_data
    from os import getcwd
    from numpy import array
    
    if getcwd().startswith('C:') or getcwd().startswith('Z:') or getcwd().startswith('M:'):
        #shpfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Magnitudes/NEAC/australia_ml_regions.shp'
        try:
            shpfile = 'Z:\\Magnitudes\\NEAC\\australia_ml_regions.shp'
            sf = shapefile.Reader(shpfile)
        except:
            shpfile = 'C:\\Users\\u56903\\OneDrive - Geoscience Australia\\Magnitudes\\shapefiles\\australia_ml_regions.shp'
            sf = shapefile.Reader(shpfile)
        
    else:
        shpfile = '/Users/trev/Documents/Manuscripts/manuscripts/2021/ml_adjustmets/shapefile/australia_ml_regions.shp'
        sf = shapefile.Reader(shpfile)
    
    shapes = sf.shapes()
    polygons = []
    for poly in shapes:
        polygons.append(Polygon(poly.points))
        
    ml_reg = get_field_data(sf, 'ML_REGION', 'str')
    
    # loop thru events and polygons
    eq_region = []
    for lo, la in zip(eqlos, eqlas):
        zone = ''
            
        for poly, reg in zip(polygons, ml_reg):
            pt = Point(lo, la)
            if pt.within(poly):
                zone = reg
                
        eq_region.append(zone)
            
    return array(eq_region)
