# -*- coding: utf-8 -*-
"""
Created on Mon May 07 11:42:07 2018

@author: u56903
"""
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
    basepath = '/nas/gemd/ehp/georisk_earthquake/modelling/national/Version_13/output/grd_files/combined'
    
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
    spectral_periods: list of strings, i.e. '0.0', '0.2'
    
    '''
    from os import path, system
    from numpy import array
    
    # set grid spectral periods for hazard curve
    spectral_periods = ['0.0', '0.05', '0.1', '0.3', '0.4', '0.5', '0.6', '0.7', \
                        '0.8', '0.9', '1.0', '1.2', '1.5', '1.7', '2.0', '2.5', \
                        '3.0', '3.5', '4.0', '4.5', '5.0']
    
    return_period = str(int(round(return_period)))
    
    # write gmt lon/lat file
    f = open('lola.txt', 'wb')
    f.write('\t'.join((str(lon), str(lat))))
    f.close()
    
    # set base path for grd location on the NAS
    basepath = '/nas/gemd/ehp/georisk_earthquake/modelling/national/Version_13/output/grd_files/combined'
    
    uhs = []
    periods = []
    
    uhstxt = "Uniform Hazard Spectra based on Geoscience Australia's 2012 National Seismic Hazard Model\n"
    uhstxt += 'Spectral acceleration for a return period of '+return_period+' years in units of g (where PGA=0.0s)\n'
    uhstxt += 'PERIOD,ACCEL\n'
    for spectral_period in spectral_periods:
        grdfile = ''.join(('avg_',spectral_period,'s_',return_period,'yr_180_30km_60km.grd'))
        grdpath = path.join(basepath, grdfile)

        # do grdtrack to extract hazard value
        system(''.join(('grdtrack lola.txt -G', grdpath, ' > lolahaz.txt')))
        
        try:
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

