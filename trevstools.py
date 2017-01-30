# -*- coding: utf-8 -*-
"""
Created on Fri Feb 01 13:44:56 2013

@author: tallen
"""

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

# converts xyz file to point shapefile
def xyz2shp(xyzfile, shpfile):
    import shapefile
    import re

    # now parse ceef
    print 'Reading xyz...'
    data = open(xyzfile).readlines()

    # now output shapefile
    print 'Making shapefile...'
    w = shapefile.Writer(shapefile.POINT)
    w.field('LON','F', 13, 6)
    w.field('LAT','F', 13, 6)
    w.field('Z','F', 13, 6)

    # now loop through points
    for line in data:
        # split space, tab or comma delimitered
        dat = re.split('\s+', line)
        if len(dat) != 3:
            dat = dat.split('\t')
        if len(dat) != 3:
            dat = dat.split(',')

        w.point(float(dat[0]), float(dat[1]))
        w.record(float(dat[0]), float(dat[1]), float(dat[2]))

    print 'Writing shapefile...'
    w.save(shpfile)

# calculate seismic moment (in N-m) from mw
def mw2m0(mw):
    return 10**(3. * mw / 2.  + 9.05)

# calculate mw from seismic moment (in N-m)
def m02mw(m0):
    from numpy import log10
    return 2 * log10(m0) / 3 - 6.03









