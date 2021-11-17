# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 11:28:54 2015

@author: tallen
"""

def cmpsps2g(cmpsps):
    from scipy.constants import g
    return cmpsps / (100. * g)

##########################################################################################
# GMICE
##########################################################################################

# do Wald et al 1999 GMICE
def pgm2mmi_wald99(value, pgm):
    from numpy import array, log10
    
    # first check if is list
    if isinstance(value, list):
        value = array([value])
        
    if pgm == 'pga': # in cm/s
        a = 3.66 
        b = -1.66
        c = 2.20
        d = 1.00
        mmisig = 1.08
        
    elif pgm == 'pgv': # in cm/s
        a = 3.47
        b = 2.35
        c = 2.10
        d = 3.40
        mmisig = 0.98
        
    mmi = a * log10(value) + b
    
    idx = mmi < 5.
    mmi[idx] = c * log10(value[idx]) + d
        
    return mmi, mmisig
    
# do Wald et al 1999 IGMICE - use in reverse
def mmi2pgm_wald99(mmi, pgm):
    from numpy import array, where
    mmi = array(mmi)
    
    # first check if is list
    if isinstance(mmi, list):
        mmi = array([mmi])

    if pgm == 'pga': # in cm/s
        a = 3.66 
        b = -1.66
        c = 2.20
        d = 1.00
        logpgmsig = 1.08 # this is not correct
        
    elif pgm == 'pgv': # in cm/s
        a = 3.47
        b = 2.35
        c = 2.10
        d = 3.40
        logpgmsig = 0.98
        
    value = 10**((mmi - b) / a)
    # now check if mmi < 5
    idx = where(mmi < 5.)[0]
    if len(idx) > 0:
        value[idx] = 10**((mmi[idx] - d) / c)
        
    return value, logpgmsig
        
# do Worden et al 2012 IGMICE
def mmi2pgm_worden12(mmi, pgm, mag, rrup):
    from numpy import array, where, isnan, log10
    
    # first check if is list
    if isinstance(mmi, list):
        mmi = array([mmi])
        
    if pgm == 'pgv': # in cm/s
        t2 = 4.56
        c1 = 3.78
        c2 = 1.47
        c3 = 2.89
        c4 = 3.16
        c5 = 0.90
        c6 = 0.00
        c7 = -0.18
        logpgmsig = 0.40
        
    if pgm == 'pga': # in cm/s**2
        t2 = 4.22
        c1 = 1.78
        c2 = 1.55
        c3 = -1.60
        c4 = 3.70
        c5 = -0.91
        c6 = 1.02
        c7 = -0.17
        logpgmsig = 0.39
        
    if pgm == 'sa03': # in cm/s**2
        t2 = 4.99
        c1 = 1.26
        c2 = 1.69
        c3 = -4.15
        c4 = 4.14
        c5 = -1.05
        c6 = 0.60
        c7 = 0.00
        logpgmsig = 0.46
        
    if pgm == 'sa10': # in cm/s**2
        t2 = 4.98
        c1 = 2.50
        c2 = 1.51
        c3 = 0.20
        c4 = 2.90
        c5 = 2.27
        c6 = -0.49
        c7 = -0.29
        logpgmsig = 0.51
        
    if pgm == 'sa30': # in cm/s**2
        t2 = 4.96
        c1 = 3.81
        c2 = 1.17
        c3 = 1.99
        c4 = 3.01
        c5 = 1.91
        c6 = -0.57
        c7 = -0.21
        logpgmsig = 0.69
    
    # do calc - first calc all
    value = 10**((mmi - c1) / c2)
    # now check if mmi > t2
    idx = where(mmi > t2)[0]
    value[idx] = 10**((mmi[idx] - c3) / c4)
    
    # do mag & dist corrections
    if not isnan(mag):
        value += c5 + c6 * log10(rrup) + c7 * mag
        
        value = 10**((mmi - c1 - c5 - c6*log10(rrup) - c7*mag) / c2)
        # now check if mmi > t2
        idx = where(mmi > t2)[0]
        value[idx] = 10**((mmi[idx] - c3 - c5 - c6*log10(rrup) - c7*mag) / c4)
        
    return value, logpgmsig

# do Worden et al 2012 GMICE
def pgm2mmi_worden12(value, pgm, mag, rrup):
    from numpy import array, log10, isnan, shape
    
    # first check if is list
    if not isinstance(value, list):
        value = array(value)
        
    if not type(value) == array:
        value = array(value)
    
    rrup = array(rrup)
        
    if pgm == 'pgv': # in cm/s
        t1 = 0.53
        c1 = 3.78
        c2 = 1.47
        c3 = 2.89
        c4 = 3.16
        c5 = 0.90
        c6 = 0.00
        c7 = -0.18
        mmisig = 0.73
        
    if pgm == 'pga': # in cm/s**2
        t1 = 1.57
        c1 = 1.78
        c2 = 1.55
        c3 = -1.60
        c4 = 3.70
        c5 = -0.91
        c6 = 1.02
        c7 = -0.17
        mmisig = 0.65
        
    if pgm == 'sa03': # in cm/s**2
        t1 = 2.21
        c1 = 1.26
        c2 = 1.69
        c3 = -4.15
        c4 = 4.14
        c5 = -1.05
        c6 = 0.60
        c7 = 0.00
        mmisig = 0.84
        
    if pgm == 'sa10': # in cm/s**2
        t1 = 1.65
        c1 = 2.50
        c2 = 1.51
        c3 = 0.20
        c4 = 2.90
        c5 = 2.27
        c6 = -0.49
        c7 = -0.29
        mmisig = 0.80
        
    if pgm == 'sa30': # in cm/s**2
        t1 = 0.99
        c1 = 3.81
        c2 = 1.17
        c3 = 1.99
        c4 = 3.01
        c5 = 1.91
        c6 = -0.57
        c7 = -0.21
        mmisig = 0.95
    
    # do calc - first calc all
    mmi = c1 + c2 * log10(value)
    # now check if value > t1
    idx = log10(value) > t1
    mmi[idx] = c3 + c4 * log10(value[idx])
    
    # do mag & dist corrections
    if not isnan(mag):
        mmi += c5 + c6 * log10(rrup) + c7 * mag

    return mmi, mmisig

# do Caprio et al 2015 Global IGMICE
def mmi2pgm_caprio15(mmi, pgm):
    from numpy import array
    
    # first check if is list
    if isinstance(mmi, list):
        mmi = array([mmi])
        
    if pgm == 'pgv': # in cm/s
        tint = 4.92
        a1 = 4.424
        b1 = 1.589
        a2 = 4.018
        b2 = 2.671
        logpgmsig = 0.60
        
    if pgm == 'pga': # in cm/s**2
        tint = 4.87
        a1 = 2.270
        b1 = 1.647
        a2 = -1.361
        b2 = 3.822
        logpgmsig = 0.40
            
    # do calc - first calc all
    pgm = 10**((mmi - a1) / b1)
    # now check if mmi > tint
    idx = mmi > tint
    pgm[idx] = 10**((mmi[idx] - a2) / b2)
    return pgm, logpgmsig
    
# do Caprio et al 2015 Global IGMICE
def pgm2mmi_caprio15(value, pgm):
    from numpy import array, log10, where
    
    # first check if is list
    if isinstance(value, list):
        value = array([value])
        
    if pgm == 'pgv': # in cm/s
        tpgm = 0.3
        a1 = 4.424
        b1 = 1.589
        a2 = 4.018
        b2 = 2.671
        mmisig = 0.90
        
    if pgm == 'pga': # in cm/s**2
        tpgm = 1.6
        a1 = 2.270
        b1 = 1.647
        a2 = -1.361
        b2 = 3.822
        mmisig = 0.70
            
    # do calc - first calc all
    mmi = a1 + b1 * log10(value)
    # now check if mmi > tint
    idx = where(value > tpgm)[0]
    if len(idx) > 0:
        mmi[idx] = a2 + b2 * log10(value[idx])
    return mmi, mmisig    

# do Atkinson & Kaka IGMICE    
def mmi2pgm_atkinson07(mmi, pgm):
    from numpy import array, log10
    
    # first check if is list
    if isinstance(mmi, list):
        mmi = array([mmi])
        
    if pgm == 'pgv': # in cm/s
        yI5 = 0.48
        c1 = 4.37
        c2 = 1.32
        c3 = 3.54
        c4 = 3.03
        mmisig = 0.80
        
    if pgm == 'pga': # in cm/s**2
        yI5 = 1.69
        c1 = 2.65
        c2 = 1.39
        c3 = -1.91
        c4 = 4.09
        mmisig = 1.01
        
    if pgm == 'sa20': # in cm/s**2
        yI5 = 1.00
        c1 = 3.72
        c2 = 1.29
        c3 = 1.99
        c4 = 3.00
        mmisig = 0.86
        
    if pgm == 'sa10': # in cm/s**2
        yI5 = 1.50
        c1 = 3.23
        c2 = 1.18
        c3 = 0.57
        c4 = 2.95
        mmisig = 0.84
       
    if pgm == 'sa03': # in cm/s**2
        yI5 = 1.92
        c1 = 2.40
        c2 = 1.36
        c3 = -1.83
        c4 = 3.56
        mmisig = 0.88
        
    # do calc - first calc all
    pgm = 10**((mmi - c1) / c2)
    # now check if mmi > t1
    idx = log10(pgm) > yI5
    pgm[idx] = 10**((mmi[idx] - c3) / c4)
    
    return pgm, mmisig
    
# do Dangkua & Cramer 2011 ENA IGMICE    
def mmi2pgm_dangkua11_ena(mmi, pgm, mag, rrup):
    from numpy import array, log10, nan
    
    # first check if is list
    if isinstance(mmi, list):
        mmi = array([mmi])
        
    if pgm == 'pgv': # in cm/s
        yt = 9999 # dummy value
        c1 = 5.13
        c2 = 1.59
        c3 = nan
        c4 = nan
        mmisig = 0.63
        
    if pgm == 'pga': # in cm/s**2
        yt = 1.95
        c1 = 2.60
        c2 = 1.58
        c3 = -1.89
        c4 = 3.89
        mmisig = 0.90
        
    if pgm == 'sa03': # in cm/s**2
        yt = 2.52
        c1 = 2.33
        c2 = 1.57
        c3 = 0.11
        c4 = 2.44
        mmisig = 0.77
        
    if pgm == 'sa10': # in cm/s**2
        yt = 9999 # dummy value
        c1 = 4.26
        c2 = 1.62
        c3 = nan
        c4 = nan
        mmisig = 0.74
       
    if pgm == 'sa20': # in cm/s**2
        yt = 9999 # dummy value
        c1 = 5.00
        c2 = 1.27
        c3 = nan
        c4 = nan
        mmisig = 0.51
        
    # do calc - first calc all
    pgm = 10**((mmi - c1) / c2)
    # now check if mmi > t1
    idx = log10(pgm) > yt
    pgm[idx] = 10**((mmi[idx] - c3) / c4)
    
    return pgm, mmisig
    
# do Newmark & Rosenblueth, 1971 - pgv in mm/s
def mmi2pgv_newmark_rosenblueth(mmi):
    from numpy import array, log2, nan
    
    pgv = ((2**mmi) * 5./7.) / 10. # comvert from mm/s to cm/s
    
    return pgv

# pgv in mm/s
def pgv2mmi_newmark_rosenblueth(pgv):
    from numpy import array, log2, nan
    
    mmi = log2((7./5.) * pgv)
    
    return mmi

# do Gaull 1979 in m/s**2
def mmi2pga_gaull(mmi):
    from numpy import array, log2, nan
    from scipy.constants import g
    
    pga = 10**((mmi / 3.1) - 2.3) * 100 # comvert from m/s**2 to cm/s**2
    
    return pga

##########################################################################################
# IPEs
##########################################################################################

def allen_etal_2012_rrup_ipe(mag, rrup, dep):
    from openquake.hazardlib.gsim.allen_2012_ipe import AllenEtAl2012
    from openquake.hazardlib.gsim.base import RuptureContext, SitesContext, DistancesContext
    from openquake.hazardlib.imt import MMI
    from openquake.hazardlib.const import StdDev
    
    from numpy import array, log10, logspace, sqrt
    
    ipe = AllenEtAl2012()
    
    sites = SitesContext()
    
    rup = RuptureContext()
    rup.mag = mag
    rup.hypo_depth = dep
    
    dists = DistancesContext()
    dists.rrup = array(rrup)
    
    mmi, sig = ipe.get_mean_and_stddevs(sites, rup, dists, MMI(), [StdDev.TOTAL])
    
    return mmi, sig[0]
    
def allen_etal_2012_rhypo_ipe(mag, rhypo, dep):
    from openquake.hazardlib.gsim.allen_2012_ipe import AllenEtAl2012Rhypo
    from openquake.hazardlib.gsim.base import RuptureContext, SitesContext, DistancesContext
    from openquake.hazardlib.imt import MMI
    from openquake.hazardlib.const import StdDev
    
    from numpy import array, log10, logspace, sqrt
    
    ipe = AllenEtAl2012Rhypo()
    
    sites = SitesContext()
    
    rup = RuptureContext()
    rup.mag = mag
    rup.hypo_depth = dep
    
    dists = DistancesContext()
    dists.rhypo = array(rhypo)
    
    mmi, sig = ipe.get_mean_and_stddevs(sites, rup, dists, MMI(), [StdDev.TOTAL])
    
    return mmi, sig[0]
    
def atkinson_wald_ceus_ipe(mag, rrup):
    from numpy import array, log10, sqrt, where, zeros_like
    
    rrup = array(rrup)
    
    c1 = 11.72
    c2 = 2.36
    c3 = 0.1155
    c4 = -0.44
    c5 = -0.002044
    c6 = 2.31
    c7 = -0.479
    h  = 17.
    Rt = 80.
    sig = 0.4
    
    R = sqrt(rrup**2 + h**2)
    B = zeros_like(R)
    idx = where(R > Rt)[0]
    B[idx] = log10(R[idx] / Rt)
        
    mmi = c1 + c2*(mag - 6) + c3*(mag - 6)**2 + c4*log10(R) + c5*R + c6*B + c7*mag*log10(R)
    
    return mmi, sig
    
def atkinson_wald_cal_ipe(mag, rrup):
    from numpy import array, log10, sqrt, where, zeros_like
    
    rrup = array(rrup)
    
    c1 = 12.27
    c2 = 2.27
    c3 = 0.1304
    c4 = -1.30
    c5 = -0.0007070
    c6 = 1.95
    c7 = -0.577
    h  = 14.
    Rt = 30.
    sig = 0.4
    
    R = sqrt(rrup**2 + h**2)
    B = zeros_like(R)
    idx = where(R > Rt)[0]
    B[idx] = log10(R[idx] / Rt)
        
    mmi = c1 + c2*(mag - 6) + c3*(mag - 6)**2 + c4*log10(R) + c5*R + c6*B + c7*mag*log10(R)
    
    return mmi, sig

def atkinson_worden_wald14_cal_ipe(mag, rhyp):
    from numpy import array, log10, sqrt, where, zeros_like
    
    rhyp = array(rhyp)
    
    c1 = 0.309
    c2 = 1.864
    c3 = -1.672
    c4 = -0.00219
    c5 = 1.77
    c6 = -0.383
    h  = 14.
    sig = 0.15
        
    R = sqrt(rhyp**2 + h**2)
    Rt = log10(R / 50.)
    B = zeros_like(R)
    idx = where(Rt > 0)[0]
    B[idx] = Rt[idx]
        
    mmi = c1 + c2*mag + c3*log10(R) + c4*R + c5*B + c6*mag*log10(R)
    
    return mmi, sig
    
def atkinson_worden_wald14_ceus_ipe(mag, rhyp, repi):
    from numpy import array, log10, sqrt, where, ones_like, zeros_like, min, max
    
    rhyp = array(rhyp)
    repi = array(repi)
    
    mmi_w = atkinson_worden_wald14_cal_ipe(mag, rhyp)[0]
    
    R1 =  ones_like(repi) * 150.
    idx = where(repi < R1)[0]
    R1[idx] = repi[idx]
    
    R2 = zeros_like(repi)
    Rtmp = 0.8 * log10(R1/50.)
    for i in range(0, len(R2)):
        if Rtmp[i] > R2[i]:
            R2[i] = Rtmp[i]
            
    mmi_e = mmi_w + 0.7 + 0.001 * repi + R2
    sig = 0.15
    
    return mmi_e, sig

def atkinson_worden_wald14_ceus_oq(mag, rhypo, dep):
    from openquake.hazardlib.gsim.atkinson_2014_ipe import AtkinsonEtAl2014CEUS
    from openquake.hazardlib.gsim.base import RuptureContext, SitesContext, DistancesContext
    from openquake.hazardlib.imt import MMI
    from openquake.hazardlib.const import StdDev
    
    from numpy import array, log10, logspace, sqrt
    
    ipe = AtkinsonEtAl2014CEUS()
    
    sites = SitesContext()
    
    rup = RuptureContext()
    rup.mag = mag
    rup.hypo_depth = dep
    
    dists = DistancesContext()
    dists.rhypo = array(rhypo)
    
    mmi, sig = ipe.get_mean_and_stddevs(sites, rup, dists, MMI(), [StdDev.TOTAL])
    
    return mmi, sig[0]

def leonard15_ipe(mw, rrup):
    from numpy import log, exp, sqrt
    
    c0 = 3.5
    c1 = 1.05
    c2 = -1.09
    c3 = 1.1
    
    return c0 + c1 * mw + c2 * log(sqrt(rrup**2 + (1+c3*exp(mw-5))**2))

##########################################################################################
# DYFI
##########################################################################################
    
# files accessed from DYFI Downloads     
def parse_usgs_dyfi_geocoded(dyfifile):
    #dyfifile = 'C:\\Users\\tallen\\Dropbox\\Moe Data\\MMI\\usgs_geocoded_cdi.txt'
    
    lines = open(dyfifile).readlines()[1:]
    
    dyfidict = []
    for line in lines:
        dat = line.split('"')
        subdat = dat[2].strip(',').split(',')
        dyfidict.append({'cdi': float(subdat[0]), 'nresp': int(subdat[1]), \
                        'repi': float(subdat[2]), 'lat': float(subdat[3]), \
                        'lon': float(subdat[4])})
                        
    return dyfidict

# files accessed from DYFI Downloads 
def parse_usgs_dyfi_zip(dyfifile):
    #dyfifile = 'C:\\Users\\tallen\\Dropbox\\Moe Data\\MMI\\usgs_zip_cdi.txt'
    
    lines = open(dyfifile).readlines()[1:]
    
    dyfidict = []
    for line in lines:
        dat = line.split('"')
        subdat = dat[2].strip(',').split(',')
        dyfidict.append({'loc': dat[1].split('::'), 'cdi': float(subdat[0]), \
                         'nresp': int(subdat[1]), 'repi': float(subdat[2]), \
                         'lat': float(subdat[3]), 'lon': float(subdat[4])})
                         	
    return dyfidict
    
# parse geojson
def parse_usgs_dyfi_geojson(dyfifile):
    import json
    from shapely.geometry import Polygon
    
    jsonfile = 'dyfi_geo_1km.geojson'
    with open(jsonfile) as f:
        data = json.load(f)
    
    dyfidict = []
    for feature in data['features']:
        points = feature['geometry']['coordinates'][0]
        
        ref_polygon = Polygon(points)
        # get the x and y coordinate of the centroid
        centtxt = ref_polygon.centroid.wkt.strip('POINT (').strip(')').split()
        
        dyfidict.append({'loc': feature['properties']['name'].split('<br>')[-1], 'cdi': feature['properties']['cdi'], \
                         'nresp': feature['properties']['nresp'], 'repi': feature['properties']['dist'], \
                         'lat': float(centtxt[1]), 'lon': float(centtxt[0])})
    
    return dyfidict

##########################################################################################
# Parse GA MMI data
##########################################################################################

def return_au_mmi_data():
    import shapefile
    #from os import path
    from numpy import array
    import datetime as dt
    from mapping_tools import get_field_data
    shpfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Ground_Motion/MMI/data/iso_p_ASCMM.shp'
    
    print('Reading MMI shapefile...')
    sf = shapefile.Reader(shpfile)
    
    # get data fields
    mmi = array(get_field_data(sf, 'INTERP_MMI', 'float'))
    eqmag = array(get_field_data(sf, 'ML_I', 'float'))
    mmilon = array(get_field_data(sf, 'IP_LONG', 'float'))
    mmilat = array(get_field_data(sf, 'IP_LAT', 'float'))
    mmisc = array(get_field_data(sf, 'WII_NEHRP', 'str'))
    eqname = array(get_field_data(sf, 'EQ_NAME', 'str'))
    eqlon = array(get_field_data(sf, 'EQ_LONG', 'float'))
    eqlat = array(get_field_data(sf, 'EQ_LAT', 'float'))
    eqdep = array(get_field_data(sf, 'DEPTH_KM', 'float'))
    eqdate = array(get_field_data(sf, 'EQ_DATE', 'str'))
    eqtime = array(get_field_data(sf, 'EQ_TIME', 'str'))
    
    eqdt = []
    # get datetime
    for eqd, eqt in zip(eqdate, eqtime):
        d = '-'.join([str(eqd[0]), str('%02d' % eqd[1]), str('%02d' % eqd[2])])
        if eqd[0] > 1800:
            if eqt.strip() == '':
                eqdt.append(dt.datetime.strptime(d+' 00:00:00', '%Y-%m-%d %H:%M:%S'))
            else:
                eqdt.append(dt.datetime.strptime(d+' '+eqt[0:8], '%Y-%m-%d %H:%M:%S'))
        else:
            eqdt.append(dt.datetime(2599, 1, 1, 1, 1))
            
    # make data dictionary
    return {'mmi':mmi, 'mmilon':mmilon, 'mmilat':mmilat, 
            'mmisc':mmisc, 'eqname':eqname, 'eqlon':eqlon, 'eqlat':eqlat,
            'eqdep':eqdep, 'eqmag':eqmag, 'datetime':array(eqdt)}

def export_historic_event_mmi(Y,m,d,H,M, outcsv): 
    import datetime
    from numpy import where, array
    #from mmi_tools import return_au_mmi_data
    
    '''
    # test 1918 Gladstone Data
    Y = 1918
    m = 6
    d = 6
    H = 18
    M = 14
    '''
    dt = datetime.datetime(Y,m,d,H,M)
    
    mmiDat = return_au_mmi_data()
    
    # find data
    idx = where((array(mmiDat['datetime']) > dt-datetime.timedelta(hours=1)) \
                & (array(mmiDat['datetime']) < dt+datetime.timedelta(hours=1)))[0]
    
    # make output csv txt
    outtxt = 'ORIGINTIME,EQLO,EQLA,EQDP,EQMAG,MMILO,MMILA,MMI,MMISC,EQNAME\n'
    for i in idx:
        outtxt += ','.join((mmiDat['datetime'][i].strftime('%Y%m%dT%H%M'), \
                            str('%0.3f' % mmiDat['eqlon'][i]), \
                            str('%0.3f' % mmiDat['eqlat'][i]), \
                            str('%0.1f' % mmiDat['eqdep'][i]), \
                            str('%0.1f' % mmiDat['eqmag'][i]), \
                            str('%0.3f' % mmiDat['mmilon'][i]), \
                            str('%0.3f' % mmiDat['mmilat'][i]), \
                            str('%0.1f' % mmiDat['mmi'][i]), \
                            mmiDat['mmisc'][i], \
                            mmiDat['eqname'][i])) + '\n'
                            
    # write to file
    f = open(outcsv,'wb')
    f.write(outtxt)
    f.close()

    





















