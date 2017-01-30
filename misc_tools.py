# -*- coding: utf-8 -*-
"""
Created on Fri Feb 01 13:44:56 2013

@author: tallen
"""
def print_functions():
    txt = open('U:/Code/pycode/fault_tools.py').readlines()
    for line in txt:
        if line.find('def') >= 0:
            print line.strip('def ').strip('\n')
            
# returns unknown int value from a string
def getintval(line,start,stop):
    from numpy import nan
    tmpstr = line[start:stop].strip()
    try:
        intval = int(tmpstr)
    except:
        intval = nan
    else:
        intval = nan
        
    return intval
    
# returns unknown float value from a string
# numdec = number of decimal places
def getfloatval(line,start1,stop1,start2,stop2):
    from numpy import nan, isnan, ceil, log10
    
    tmpval1 = getintval(line,start1,stop1)
    tmpval2 = getintval(line,start2,stop2)
    # get length of decimal
    if tmpval2 <= 0:
        logdec = 1
    else:
        logdec = ceil(log10(tmpval2))
    
    if isnan(tmpval1) == True:
        floatval = nan
    else:
        if isnan(tmpval2) == True:
            floatval = tmpval1 * 1.
        else:
            floatval = tmpval1 + tmpval2 / 10.**logdec
            
    return floatval

    
def checkfloat(floatstr):
    try:
        return float(floatstr)
    except:
        from numpy import nan
        return nan
        
def checkint(intstr):
    try:
        return int(intstr)
    except:
        from numpy import nan
        return nan   
    
# calculate seismic moment (in N-m) from mw
def mw2m0(mw):
    return 10**(3. * mw / 2.  + 9.05)

# calculate mw from seismic moment (in N-m)
def m02mw(m0):   
    from numpy import log10
    return 2 * log10(m0) / 3 - 6.03

# converts m/s^2 to g    z
def mpss2g(mpss):
    return mpss / 9.81
    
# converts m/s^2 to %g        
def mpss2percentg(mpss):
    return mpss *100 / 9.81
    
# converts cm/s^2 to g    
def cmpss2g(cmpss):
    return (cmpss / 100) / 9.81

# converts cm/s^2 to %g    
def cmpss2percentg(cmpss):
    return (cmpss / 100) * 100 / 9.81
    
# converts cm/s^2 to %g    
def percentg2cmpss(percentg):
    return (percentg / 100) * 9.81 * 100  
    
# create flat dict with same keys
def flatten_dict(indict):
    # creat dict with same keys
    newdict = {}
    
    for key in indict[0].keys():
        newdict[key] = []
        for rec in indict:
            newdict[key].append(rec[key])
            
    return newdict
    
# find indexes in an array of dictionaries (like matlab find)
def find_eq(dictarray, key, value):
    from numpy import array
    
    return array([i for i, dic in enumerate(dictarray) if dic[key] == value])
    
def find_eq_date(dictarray, date):
    # date = integer in format YYYYMMDDHHMN

    from numpy import array
    import datetime as dt
    
    evtime = dt.datetime.strptime(str(date), '%Y%m%d%H%M')
    
    deltat = dt.timedelta(minutes=1.)
    
    return array([i for i, dic in enumerate(dictarray) \
                  if dic['datetime'] > evtime-deltat \
                  and dic['datetime'] < evtime+deltat])

# get data & stats between in bins
'''
eg:
	binsize = 5.
	bins = arange(2.5, 45., binsize)
	xdat = numpy array of x values
	yres = numpy array of y residuals
'''
def get_binned_stats(bins, xdat, yres):
    from numpy import array, diff, nanstd, where, isfinite
    from scipy.stats import nanmedian
    
    binsize = diff(bins)[0]

    medres = []
    stdres = []
    medx = []
    nperbin = []

    yres = array(yres)

    halfbin = binsize / 2.0

    for bin in bins:
        index = array(where((xdat >= bin-halfbin) & (xdat < bin+halfbin))[0])
        medres.append(nanmedian(yres[index]))
        stdres.append(nanstd(yres[index]))
        medx.append(nanmedian(xdat[index]))
        nperbin.append(len(index))
        
    idx = where(isfinite(medres))[0]

    return array(medres)[idx], array(stdres)[idx], array(medx)[idx], bins[idx], array(nperbin)[idx]

def get_binned_stats_mean(bins, xdat, yres):
    from numpy import array, diff, nanstd, where, isfinite
    from scipy.stats import nanmean
    
    binsize = diff(bins)[0]

    medres = []
    stdres = []
    medx = []

    yres = array(yres)

    halfbin = binsize / 2.0

    for bin in bins:
        index = array(where((xdat >= bin-halfbin) & (xdat < bin+halfbin))[0])
        medres.append(nanmean(yres[index]))
        stdres.append(nanstd(yres[index]))
        medx.append(nanmean(xdat[index]))
        
    idx = where(isfinite(medres))[0]

    return array(medres)[idx], array(stdres)[idx], array(medx)[idx], bins[idx]

# get weighted std, from http://stackoverflow.com/questions/2413522/weighted-standard-deviation-in-numpy
def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    from numpy import average, sqrt
    
    wtaverage = average(values, weights=weights)
    variance = average((values-wtaverage)**2, weights=weights)  # Fast and numerically precise
    return (wtaverage, sqrt(variance))

# plot text in figures
def plttext(plt ,xpos, ypos, txt, **kwargs):
    from numpy import array
    
    fsize = 11
    col = 'k'
    
    for key in ['fsize', 'col']:
        if key in kwargs:
            # set font size
            if key == 'fsize':
                fsize = kwargs[key]
            if key == 'col':
                col = kwargs[key]

    #plt.annotate(txt, (0,0), fontsize=fsize, color=col, xytext=(xpos, ypos), xycoords='data')

# to extrapolate linear interp
def extrap1d(interpolator):
    '''
    To use:
    
    from scipy.interpolate import interp1d
    x = arange(0,10)
    y = exp(-x/3.0)
    f_i = interp1d(x, y)
    f_x = extrap1d(f_i)
    
    print f_x([9,10])
    '''
    
    from scipy.interpolate import interp1d
    from scipy import arange, array, exp
    
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]:
            return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else:
            return interpolator(x)

    def ufunclike(xs):
        return array(map(pointwise, array(xs)))

    return ufunclike

# convert YYYYMMDD to Julian day
def ymd2doy(yyyy, mm, dd):
    from datetime import datetime
    
    dt = datetime(yyyy, mm, dd) 
    
    return dt.strftime('%j')
    
# convert Julian day 2 YYYYMMDD 
def doy2ymd(yyyy, doy):
    from datetime import datetime
    
    dt = datetime.strptime(str(yyyy)+str(doy), '%Y%j')
    
    return dt.strftime('%Y%m%d')
    
# time delta to days, hours, mins
def timedelta2days_hours_minutes(td):
            return td.days, td.seconds//3600, (td.seconds//60)%60

# list files with a given extension
def listdir_extension(folder, extension):
    from os import listdir
    
    files = []
    #extension = extension.strip('.')
    files.append([each for each in listdir(folder) if each.endswith('.'+extension)])
    
    return files[0]
    
# for getting log xy plotting locations
def get_log_xy_locs(lims, percent_loc):
   from numpy import log10, diff
   
   loglims = log10(lims)
   logdiff = diff(loglims)[0]
   
   return 10**(min(loglims) + logdiff * percent_loc) 
   
# from http://stackoverflow.com/questions/2413522/weighted-standard-deviation-in-numpy
def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    from numpy import average, sqrt
    
    wtaverage = average(values, weights=weights)
    variance = average((values-wtaverage)**2, weights=weights)  # Fast and numerically precise
    return (wtaverage, sqrt(variance))
   

