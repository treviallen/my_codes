# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 10:29:43 2013

Requires input xyz file - outputs it to point shapefile

Expects input as x, y, z columns

@author: tallen
"""

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
#    w.field('LON','F', 13, 6)
#    w.field('LAT','F', 13, 6)
    w.field('Z','F', 13, 6)
    
    # now loop through points
    for line in data:
        # split space, tab or comma delimitered
        dat = re.split('\s+', line)
#        if len(dat[0]) != 3:
#            dat = dat[0].split('\t')
#        if len(dat[0]) != 3:
#            dat = dat[0].split(',')
            
        w.point(float(dat[0]), float(dat[1]))
#        w.record(float(dat[0]), float(dat[1]), float(dat[2]))
        w.record(float(dat[2]))
            
    print 'Writing shapefile...'
    w.save(shpfile)