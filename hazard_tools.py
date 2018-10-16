# -*- coding: utf-8 -*-
"""
Created on Mon May 07 11:42:07 2018

Useful functions for interrogating hazard data

@author: u56903
"""
def beta2bval(beta):
    from numpy import log10, exp
    return log10(exp(beta))

def bval2beta(bval):
    from numpy import log
    return log(10**bval)

# eqns taken from: https://earthquake.usgs.gov/hazards/learn/basics.php
def get_probability_from_percent_chance(percent_chance, investigation_time):
     from numpy import log

     p0 = 1 - (percent_chance / 100.)
     n = -log(p0)
     probability = n / investigation_time
     return_period = 1. / probability

     return return_period, probability

def get_percent_chance_from_return_period(return_period, investigation_time):
    from numpy import exp

    n = (1. / return_period) * investigation_time
    p0 = exp(-n)
    percent_chance = 100*(1 - p0)

    return percent_chance

def get_probability_from_obs_and_return_period(return_period, investigation_time):
    from numpy import exp
    return 1 - exp(-investigation_time/return_period)

def poe_invtime_to_annual(imtls, poe, investigation_time):
    from numpy import array, log
    
    # now get annualised curves
    P0 = 1 - array(imtls)
    n = -1*log(P0)
    annual_probs = n / investigation_time
    
    tmpdict[ut+'_probs_annual'] = annual_probs
    tmpdict[ut+'_probs_invtime'] = hazcurve


def get_nsha12_hazard_curve(lon, lat, spectral_period):
    '''
    spectral_periods: list of strings, i.e. '0.0', '0.2'
    
    '''
    from os import path, system
    from numpy import array
    
    # set grid return periods for hazard curve
    return_periods = ['100', '250', '500', '800', '1000', '1500', '2000', '2500', \
                      '3000', '5000', '10000', '20000', '50000', '100000']
    
    # write gmt lon/lat file
    f = open('lola.txt', 'wb')
    f.write('\t'.join((str(lon), str(lat))))
    f.close()
    
    # set base path for grd location on the NAS
    basepath = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/national/Version_13/output/grd_files/combined'
    
    hazArray = []
    probabilities = []
    return_period_nums = []
    exceedances = []
    for return_period in return_periods:
        grdfile = ''.join(('avg_',spectral_period,'s_',return_period,'yr_180_30km_60km.grd'))
        grdpath = path.join(basepath, grdfile)

        # do grdtrack to extract hazard value
        system(''.join(('grdtrack lola.txt -G', grdpath, ' > lolahaz.txt')))
        
        try:
            # parse in hazard value
            hazval = float(open('lolahaz.txt').read().strip().split('\t')[-1])
            
            # append to hazArray
            hazArray.append(hazval)
            
            # get probability
            investigation_time = 50.
            percent_chance = get_percent_chance_from_return_period(float(return_period), investigation_time)
            return_period_num, probability = get_probability_from_percent_chance(percent_chance, investigation_time)
            
            # append to arrays
            return_period_nums.append(return_period_num)
            probabilities.append(probability)
            exceedances.append(percent_chance)
            
        except:
            print 'File not found:', grdpath
            
        
    return array(hazArray), array(return_period_nums), array(exceedances)

def get_nsha12_hazard_spectra(lon, lat, return_period):
    '''
    return_period: list of strings, i.e. 100, 250, 500, etc
    
    
    '''
    from os import path, system
    from numpy import array
    
    # set grid spectral periods for hazard curve
    spectral_periods = ['0.0', '0.05', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', \
                        '0.8', '0.9', '1.0', '1.2', '1.5', '1.7', '2.0', '2.5', \
                        '3.0', '3.5', '4.0', '4.5', '5.0']
    
    return_period = str(int(round(return_period)))
    
    # write gmt lon/lat file
    f = open('lola.txt', 'wb')
    f.write('\t'.join((str(lon), str(lat))))
    f.close()
    
    # set base path for grd location on the NAS
    basepath = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/national/Version_13/output/grd_files/combined'
    
    uhs = []
    periods = []
    
    uhstxt = "Uniform Hazard Spectra based on Geoscience Australia's 2012 National Seismic Hazard Model\n"
    uhstxt += 'Spectral acceleration for a return period of '+return_period+' years in units of g (where PGA=0.0s)\n'
    uhstxt += 'PERIOD,ACCEL\n'
    for spectral_period in spectral_periods:
        grdfile = ''.join(('avg_',spectral_period,'s_',return_period,'yr_180_30km_60km.grd'))
        grdpath = path.join(basepath, grdfile)

        try:
            # do grdtrack to extract hazard value
            system(''.join(('grdtrack lola.txt -G', grdpath, ' > lolahaz.txt')))
        
        
            # parse in hazard value
            hazval = float(open('lolahaz.txt').read().strip().split('\t')[-1])
            
            # append to hazArray
            uhs.append(hazval)
            periods.append(float(spectral_period))
            
            # add text to output
            uhstxt += ','.join((spectral_period, str('%0.4f' % hazval))) + '\n'
            
        except:
            print 'File not found:', grdpath
            
    # write to file
    uhsFile = '_'.join(('NSHM12_UHS',return_period,str(lon),str(lat)))+'.csv'
    f = open(uhsFile, 'wb')
    f.write(uhstxt)
    f.close()    
        
    return array(periods), array(uhs)
    
def return_AS1170_4_shape(periods, siteclass):
    '''
    siteclass = A-B
    '''
    from numpy import array
    
    shp1170 = []
    
    if siteclass.upper() == 'A':
        for t in periods:
            if t <= 0.1:
                shp1170.append(0.8 + 15.5*t)
            elif t > 0.1 and t <= 1.5:
                shp1170.append(min(0.704/t, 2.35))
            else:
                shp1170.append(1.056 / t**2)
                
    if siteclass.upper() == 'B':
        for t in periods:
            if t <= 0.1:
                shp1170.append(1.0 + 19.4*t)
            elif t > 0.1 and t <= 1.5:
                shp1170.append(min(0.88/t, 2.94))
            else:
                shp1170.append(1.32 / t**2)
    
    if siteclass.upper() == 'C':
        for t in periods:
            if t <= 0.1:
                shp1170.append(1.3 + 23.8*t)
            elif t > 0.1 and t <= 1.5:
                shp1170.append(min(1.25/t, 3.68))
            else:
                shp1170.append(1.874 / t**2)
                
    if siteclass.upper() == 'D':
        for t in periods:
            if t <= 0.1:
                shp1170.append(1.1 + 25.8*t)
            elif t > 0.1 and t <= 1.5:
                shp1170.append(min(1.98/t, 3.68))
            else:
                shp1170.append(2.97 / t**2)
                
    if siteclass.upper() == 'E':
        for t in periods:
            if t <= 0.1:
                shp1170.append(1.1 + 25.8*t)
            elif t > 0.1 and t <= 1.5:
                shp1170.append(min(3.08/t, 3.68))
            else:
                shp1170.append(4.62 / t**2)
                
    return array(shp1170)
           

# function to get hazard curves for a list of cities
def get_nsha18_city_curve(citylist, hazcurvefile):
    '''
    # make path to hazcurvefile
    #hazcurvefile = '/Users/tallen/Documents/Geoscience_Australia/NSHA2018/source_models/complete_model/final/results_fractiles/hazard_curve-mean-PGA_1.csv'
    '''
    from tools.oq_tools import return_annualised_haz_curves
    from numpy import around
    
    ##############################################################################
    # parse site file
    ###############################################################################
    sitelistfile = '/Users/tallen/Documents/Geoscience_Australia/NSHA2018/shared/nsha_cities.csv'
    lines = open(sitelistfile).readlines()
    places = []
    place_lat = []
    place_lon = []
    
    for line in lines:
        dat = line.strip().split(',')
        place_lon.append(float(dat[0]))
        place_lat.append(float(dat[1]))
        places.append(dat[2])
    
    ###############################################################################
    # parse first job file to define plotting order
    ###############################################################################
    
    # get data from first job
    siteDict, imls, investigation_time = return_annualised_haz_curves(hazcurvefile)
    
    # loop thru sites in first job file and fill pltDict
    pltDict = []
    for sd in siteDict:
        pltTrue = False
        
        ###############################################################################
        # loops thru places to get title - check if want to plot
        ###############################################################################
        for place, plon, plat in zip(places, place_lon, place_lat):
            if around(plon, decimals=2) == around(sd['lon'], decimals=2) \
               and around(plat, decimals=2) == around(sd['lat'], decimals=2):
                
                # now loop through places we want to plot
                for pp in citylist:
                    if pp == place:
                        sd['place'] = place
                        sd['imls'] = imls
                        #h1 = plt.semilogy(imls, sd1['poe_probs_annual'], color=cs[p*2], lw=2.0, label=label_place+' (F)')
                        
                        pltDict.append(sd)
    
    return pltDict

# partial pythonisation of Nico's code
# does not do iteration - just for plotting purposes
def calc_risk_integral(RTGM, beta, SAs, Probs):
    from scipy.stats import norm, lognorm
    from numpy import array, arange, exp, log, trapz, interp
    from scipy import interpolate
    from misc_tools import extrap1d
    
    FRAGILITY_AT_RTGM = 0.10 
    BETA = 0.6
    AFE4UHGM = - log( 1 - 0.02 )/ 50 # exceedance frequency for 1/2475 yrs
    TARGET_RISK = - log( 1 - 0.01 ) / 50
    
    '''
    SAs = array([ 0.1613, 0.1979, 0.2336, 0.3385, 0.4577, 0.5954, 0.7418, 0.7905, 0.9669, 1.1697])
    Probs = array([0.02, 0.01375, 0.01, 0.00445, 0.0021, 0.001, 0.0005, 0.000404, 0.0002, 0.0001])
    '''
    # get uniform hazard at 1/2475
    UHGM = exp((interp(log(AFE4UHGM), log(Probs[::-1]), log(SAs[::-1]))))
    
    # up sample hazard curve
    UPSAMPLING_FACTOR = 1.05
    SMALLEST_SA = min([min(SAs), UHGM/10])
    LARGEST_SA = max([max(SAs), UHGM*10])
    
    
    upSAs = exp(arange(log(SMALLEST_SA),log(LARGEST_SA),log(UPSAMPLING_FACTOR)))
    f_i = interpolate.interp1d(log(SAs), log(Probs))
    f_x = extrap1d(f_i)
    upProbs = exp(f_x(log(upSAs)))
    '''
    upSAs = SAs
    upProbs = Probs
    '''
    # get fragility curve
    FragilityCurve = {}
    FragilityCurve['Median'] = RTGM / exp( norm.ppf( FRAGILITY_AT_RTGM ) * BETA )  
    FragilityCurve['PDF'] = lognorm.pdf(upSAs,BETA,scale=(FragilityCurve['Median']))
    FragilityCurve['CDF'] = lognorm.cdf(upSAs,BETA,scale=(FragilityCurve['Median']))
    FragilityCurve['SAs'] = upSAs
    FragilityCurve['Beta'] = BETA 
    
    # do risk integral
    Integrand = FragilityCurve['PDF'] * upProbs  
    Risk = trapz(Integrand, upSAs)
    
    # calculate collapse probability
    CollapseProb = 1 - exp(-50 * Risk)
    
    RiskCoefficient = RTGM / UHGM
    
    return upProbs, upSAs, FragilityCurve, Integrand, CollapseProb
