# -*- coding: utf-8 -*-
"""
Created on Fri Feb 01 13:44:56 2013

@author: tallen
"""
def print_functions():
    txt = open('U:/Code/pycode/fault_tools.py').readlines()
    for line in txt:
        if line.find('def') >= 0:
            print(line.strip('def ').strip('\n'))
            
# flattens a single key value from a list of dictionaries
def dictlist2array(dictList, key):
    from numpy import array
    
    flatList = []
    for dl in dictList:
        flatList.append(dl[key])
        
    return array(flatList)

# check if key exists in dict
def check_key_exists(dict, key):
    if key in dict.keys():
        return True
    else:
        return False

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
    #from numpy import array, diff, where, isfinite, nanmedian, nanstd, isnan
    from numpy import array, diff, where, isfinite, median, std, isnan, delete
    
    yres = array(yres)
    xdat = array(xdat)
    
    # strip nan values
    idx = where(isnan(yres))[0]
    yres = delete(yres, idx)
    xdat = delete(xdat, idx)
    
    binsize = diff(bins)[0]

    medres = []
    stdres = []
    medx = []
    meanx = []
    nperbin = []

    halfbin = binsize / 2.0

    for bin in bins:
        index = array(where((xdat >= bin-halfbin) & (xdat < bin+halfbin))[0])
        #medres.append(nanmedian(yres[index]))
        #stdres.append(nanstd(yres[index]))
        #medx.append(nanmedian(xdat[index]))
        medres.append(median(yres[index]))
        stdres.append(std(yres[index]))
        medx.append(median(xdat[index]))
        meanx.append(median(xdat[index]))
        nperbin.append(len(index))
        
    idx = where(isfinite(medres))[0]

    return array(medres)[idx], array(stdres)[idx], array(medx)[idx], bins[idx], array(nperbin)[idx]

def get_binned_stats_mean(bins, xdat, yres):
    #from numpy import array, diff, nanstd, where, isfinite, nanmean, isnan
    from numpy import array, diff, std, where, isfinite, mean, median, isnan, delete
    
    yres = array(yres)
    xdat = array(xdat)
    
    # strip nan values
    idx = where(isnan(yres))[0]
    yres = delete(yres, idx)
    xdat = delete(xdat, idx)   
    binsize = diff(bins)[0]

    meanres = []
    stdres = []
    medx = []
    nperbin = []

    halfbin = binsize / 2.0

    for bin in bins:
        index = array(where((xdat >= bin-halfbin) & (xdat < bin+halfbin))[0])
        #medres.append(nanmean(yres[index]))
        #stdres.append(nanstd(yres[index]))
        #medx.append(nanmean(xdat[index]))
        meanres.append(mean(yres[index]))
        stdres.append(std(yres[index]))
        medx.append(mean(xdat[index]))
        nperbin.append(len(index))
        
    idx = where(isfinite(meanres))[0]
    return array(meanres)[idx], array(stdres)[idx], array(medx)[idx], bins[idx], array(nperbin)[idx]
    
def get_binned_stats_meanx(bins, xdat, yres):
    #from numpy import array, diff, where, isfinite, nanmedian, nanstd, isnan
    from numpy import array, diff, where, isfinite, median, std, isnan, delete
    
    yres = array(yres)
    xdat = array(xdat)
    
    # strip nan values
    idx = where(isnan(yres))[0]
    yres = delete(yres, idx)
    xdat = delete(xdat, idx)
    
    binsize = diff(bins)[0]

    medres = []
    stdres = []
    medx = []
    meanx = []
    nperbin = []

    halfbin = binsize / 2.0

    for bin in bins:
        index = array(where((xdat >= bin-halfbin) & (xdat < bin+halfbin))[0])
        #medres.append(nanmedian(yres[index]))
        #stdres.append(nanstd(yres[index]))
        #medx.append(nanmedian(xdat[index]))
        medres.append(median(yres[index]))
        stdres.append(std(yres[index]))
        medx.append(median(xdat[index]))
        meanx.append(median(xdat[index]))
        nperbin.append(len(index))
        
    idx = where(isfinite(medres))[0]

    return array(medres)[idx], array(stdres)[idx], array(meanx)[idx], bins[idx], array(nperbin)[idx]

def trimmed_lin_reg(x, y, nstd = 2):
    from scipy.stats import linregress
    from numpy import nanstd, nan, isnan, where
    
    nidx = where(isnan(y) == False)[0]
    init_fit = linregress(x[nidx], y[nidx])
    init_slope = init_fit[0]
    init_inter = init_fit[1]
    
    # do trimmed lin reg
    mod_pred = init_inter + init_slope * x
    mod_res = y - mod_pred
    mod_std = nanstd(mod_res)
    
    # remove vals outside nstd
    idx = where(abs(mod_res) <= nstd * mod_std)[0]
    
    # re-regress
    if len(idx) > 0:
       lts_fit = linregress(x[idx], y[idx])
       lts_slope = lts_fit[0]
       lts_inter = lts_fit[1]
       lts_rval = lts_fit[2]
       
    else:
       lts_slope = nan
       lts_inter = nan
       lts_rval  = nan
       
    # return fit and data
    return lts_slope, lts_inter, lts_rval, x[idx], y[idx]

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
    
    print(f_x([9,10])
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

##############################################################################
# datetime functions
##############################################################################

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
    
def datetime2epochtime(Y,m,d,H,M,S):
    from datetime import datetime
    
    print(datetime(Y,m,d,H,M,S).strftime('%s'))


# time delta to days, hours, mins
def timedelta2days_hours_minutes(td):
    return td.days, td.seconds//3600, (td.seconds//60)%60

# from https://stackoverflow.com/questions/6451655/python-how-to-convert-datetime-dates-to-decimal-years
def toYearFraction(date):
    from datetime import datetime as dt
    import time

    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction

# list files with a given extension
def listdir_extension(folder, extension):
    from os import listdir
    
    files = []
    #extension = extension.strip('.')
    files.append([each for each in listdir(folder) if each.endswith('.'+extension)])
    
    return files[0]
    
def get_folder_list(infolder):
    from os import path, listdir
    
    # get all objects in folder
    objects = listdir(infolder)
    
    # now, check which of the objects is a directory
    folder_list = []
    for obj in objects:
        obj_path = path.join(infolder, obj)
        
        if path.isdir(obj_path):
            folder_list.append(obj)
            
    return folder_list 

# list files with a common file prefix
def listdir_file_prefix(folder, file_prefix):
    from os import listdir
    
    files = []
    #extension = extension.strip('.')
    files.append([each for each in listdir(folder) if each.startswith(file_prefix)])
    
    for f in files[0]:
        print(f)
    
    return files[0]
    
# for getting log xy plotting locations
def get_log_xy_locs(lims, fraction_loc):
   from numpy import log10, diff
   
   loglims = log10(lims)
   logdiff = diff(loglims)[0]
   
   return 10**(min(loglims) + logdiff * fraction_loc) 
   
# set one decimal point on x & y axes
def fmt_axes_tick_decimals(ax):
    '''
    e.g., ax = plt.subplot(111)
    '''
    
    from matplotlib.ticker import FormatStrFormatter
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

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
    
# code from here: https://gist.github.com/bertspaan/8220892
def dbf2csv(dbffile):
    import csv
    from dbfpy import dbf
    import sys
    
    filename = sys.argv[1]
    if filename.endswith('.dbf'):
        print("Converting %s to csv" % filename)
        csv_fn = filename[:-4]+ ".csv"
        with open(csv_fn,'wb') as csvfile:
            in_db = dbf.Dbf(filename)
            out_csv = csv.writer(csvfile)
            names = []
            for field in in_db.header.fields:
                names.append(field.name)
            out_csv.writerow(names)
            for rec in in_db:
                out_csv.writerow(rec.fieldData)
            in_db.close()
            print("Done...")
    else:
      print("Filename does not end with .dbf")

def remove_last_cmap_colour(cmap):
    from matplotlib.colors import LinearSegmentedColormap
    cmap_list = []
    inc = 1. / (cmap.N-1)
    for i in range(0, cmap.N-1):
        cmap_list.append(cmap(i*inc))
        
    return LinearSegmentedColormap.from_list('cmap2', cmap_list, N=cmap.N-1)
    
def remove_first_cmap_colour(cmap):
    from matplotlib.colors import LinearSegmentedColormap
    cmap_list = []
    inc = 1. / (cmap.N-1)
    for i in range(1, cmap.N):
        cmap_list.append(cmap(i*inc))
        
    return LinearSegmentedColormap.from_list('cmap2', cmap_list, N=cmap.N-1)


def manual_colour_list(ncolours):
    from numpy import array
    colours =  array(['#a020f0', '#8a2be2', '#483d8b', '#000080', '#0000ff', '#4682b4', \
                      '#00ced1', '#00ffff', '#66cdaa', '#2e8b57', '#006400', '#228b22', \
                      '#7cfc00', '#ffff00', '#ffd700', '#ff8c00', '#ff4500', '#ff0000', \
                      '#b22222', '#cd5c5c', '#ff69b4', '#ff1493', '#ffc0cb', '#bebebe', \
                      '#708090', '#2f4f4f', '#000000', '#eedd82'])
              
    # get steps
    step = int(round((len(colours)-1) / (ncolours-1)))
    
    # get indexes
    idx = range(0, len(colours), step)
    
    return colours[::-1]

# from: http://scipy.github.io/old-wiki/pages/Cookbook/SavitzkyGolay
def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    """Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial
    
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except(ValueError, msg):
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

# from https://scipy-cookbook.readthedocs.io/items/SignalSmooth.html
def smooth(x,window_len=11,window='hanning'):
    import numpy
    
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise(ValueError, "smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise(ValueError, "Input vector needs to be bigger than window size.")


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise(ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    s=numpy.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='valid')
    
    return y


# from: https://gist.github.com/jakevdp/91077b0cae40f8f8244a
def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:
    import matplotlib.pyplot as plt
    import numpy as np
    
    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    
    return base.from_list(cmap_name, color_list, N)

# returns N x 3 RGB array
def cmap2rgb(cmap, ncolours):
    '''
    ncolours = N colours in rgb array (integer)
    '''
    import matplotlib.pyplot as plt
    from numpy import linspace, newaxis
    
    x = linspace(0.0, 1.0, ncolours)
    
    return plt.get_cmap(plt.get_cmap(cmap))(x)[newaxis, :, :3]

"""
def shakemap_rgb2cmap():
    from matplotlib.colors import ListedColormap
    
    MMI_COLORS_24BIT = [
    (255, 255, 255),
    (191, 204, 255),
    (160, 230, 255),
    (128, 255, 255),
    (122, 255, 147),
    (255, 255, 0),
    (255, 200, 0),
    (255, 145, 0),
    (255, 0, 0),
    (180, 0, 0)] 
    
    MMI_COLORS = [tuple(x / 255 for x in color) for color in MMI_COLORS_24BIT]
    cmap = ListedColormap(MMI_COLORS)
    
    # usage
    #mmi_data.plot(ax=ax, cmap=cmap, column='intensity', vmin=0.5, vmax=10.5, legend=True)
    
    return cmap
"""
    
def get_mpl2_colourlist():
    return ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
            '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
            '#bcbd22', '#17becf']
            
def get_ga_master_colours():
    return ['#006983', '#a33f1f', '#72c7e7', '#5e6a71']
    
def get_ga_secondary_colours():
    return ['#6e7645', '#ca7700', '#988642', '#a5d867', '#003145']

def get_ga_master_colours_2022():
    return ['#00718b', '#082e41', '#e6edef']
    
def get_ga_secondary_colours_2022():
    return ['#606f74', '#773775', '#637c6b', '#0b5e4a', '#cb6c37', '#b43b3b']
    
def get_line_styles():
    return ['-', '--', '-.', (0, (3, 5, 1, 5, 1, 5)), (0, (3, 1, 1, 1, 1, 1)), (0, (3, 5, 1, 5))]
    
'''    
# function to convert csv to json
# copied from: https://pythonexamples.org/python-csv-to-json/
'''
def csv2json(csvFilePath, jsonFilePath):
    '''
    csvFilePath = input csv file
    jsonFilePath = output json file
    '''
    
    import csv 
    import json 
    jsonArray = []
      
    #read csv file
    with open(csvFilePath, encoding='utf-8') as csvf: 
        #load csv file data using csv library's dictionary reader
        csvReader = csv.DictReader(csvf) 

        #convert each csv row into python dict
        for row in csvReader: 
            #add this python dict to json array
            jsonArray.append(row)
  
    #convert python jsonArray to JSON String and write to file
    with open(jsonFilePath, 'w', encoding='utf-8') as jsonf: 
        jsonString = json.dumps(jsonArray, indent=4)
        jsonf.write(jsonString)

# from https://stackoverflow.com/questions/27327513/create-pdf-from-a-list-of-images
def pnglist2pdf(pnglist, outpdf):
    from fpdf import FPDF
    pdf = FPDF()
    # imagelist is the list with all image filenames
    for image in pnglist:
        pdf.add_page()
        pdf.image(image,x,y,w,h)
    pdf.output(outpdf, "F")             