# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 12:29:06 2013

1) get_grd_extent(grdfile)
2) parse_gmt_polys(polyfile)

@author: tallen
"""

def get_grd_extent(**kwargs):
    from scipy.io.netcdf import netcdf_file
    from numpy import nan

    grdfile = ''
    nc = nan

    for key in ('grdfile', 'nc'):
        if key in kwargs:
            # min value
            if key == 'grdfile':
                grdfile = kwargs[key]

            # set scaling relation
            if key == 'nc':
                nc = kwargs[key]

    # read file
    if grdfile != '':
        nc = netcdf_file(grdfile, 'r')

    # get keys
    keys = nc.variables.keys()

    # loop through keys assumming x, y, z in degrees
    for key in keys:
        if key == 'x' or key == 'lon':
            xrng = []
            keyvals = nc.variables[key]
            xrng = keyvals.actual_range

        if key == 'y' or key == 'lat':
            yrng = []
            keyvals = nc.variables[key]
            yrng = keyvals.actual_range

        if key == 'z':
            zrng = []
            keyvals = nc.variables[key]
            zrng = keyvals.actual_range

    return xrng, yrng, zrng

# do GMT grdtrack instead
def gmt_grdtrack(lon, lat, grdfile):
    from os import system
    from numpy import nan
    
    # write lon/lat file
    f = open('lonlat.txt', 'wb')
    f.write(','.join((str(lon), str(lat))))
    f.close()
    
    # call grdtrack
    system('gmt grdtrack lonlat.txt -G' + grdfile + ' > lonlatzval.txt')
    
    # get zval from file
    zline = open('lonlatzval.txt').readline()
    if len(zline) > 0:
        zval = float(zline.split()[-1])
    else:
        zval = nan
    
    return zval
    
# do GMT grdtrack instead
def gmt_grdtrack_list(lons, lats, grdfile):
    from os import system
    from numpy import array, nan
    
    # frist make text
    txt = ''
    for lo, la in zip(lons, lats):
        txt += ','.join((str(lo), str(la))) + '\n'
        
    # write lon/lat file
    f = open('lonlats.txt', 'wb')
    f.write(txt)
    f.close()
    
    # call grdtrack
    system('gmt grdtrack lonlats.txt -G' + grdfile + ' -N > lonlatzval.txt')
    
    # get zvals from file
    zlines = open('lonlatzval.txt').readlines()
    zvals = []
    lons = []
    lats = []
    for line in zlines:        
        if len(line) > 0:
            zvals.append(float(line.split()[-1]))
            lons.append(float(line.split()[0]))
            lats.append(float(line.split()[1]))
        else:
            zvals.append(nan)
    
    return array(lons), array(lats), array(zvals)

# python version of GMT grdtrack - very dodgy, but works ok
def pygrdtk(alon, alat, **kwargs):
    '''
    alon & alat must be two arrays of equal length
    kwargs:
    	grdfile: name of netcdf file
    	nc:  pre-read netcdf data

    '''

    '''
    test data
    '''
    #alat = 33.682
    #alon = 136.204

    from netCDF4 import Dataset as NetCDFFile
    from scipy.interpolate import interp2d, interp1d
    from numpy import shape, where, hstack, nanmean, argsort, sort, nan, isnan, array
    from os import system, remove, path

    # get kwargs
    grdfile = ''
    nc = nan

    #grdfile = '/Users/trev/Documents/Data/Slab_1.0/izu_slab1.0_dipclip.grd'

    for key in ('grdfile', 'nc'):
        if key in kwargs:
            # min value
            if key == 'grdfile':
                grdfile = kwargs[key]

            # set scaling relation
            if key == 'nc':
                nc = kwargs[key]

    # force to pixel node if GMT installed
    pixelnode = False
    if grdfile != '':
        try:
            system('grdsample ' + grdfile + ' -Gtmp.grd -F')
            nc = NetCDFFile('tmp.grd')
            pixelnode = True
        except:
            nc = NetCDFFile(grdfile)

        xrng, yrng, zrng = get_grd_extent(grdfile=grdfile)
    else:
       xrng, yrng, zrng = get_grd_extent(nc=nc)

    # if float, convert to numpy array
    try:
       length = len(alat)
    except:
       tmplat = []
       tmplon = []
       tmplat.append(alat)
       tmplon.append(alon)
       alat = array(tmplat)
       alon = array(tmplon)

    zdat = nc.variables['z'][:]
    xdat = nc.variables['x'][:]
    ydat = nc.variables['y'][:]

    # get grid info
    grdshp = shape(zdat)
    if pixelnode == False:
        xres = (xrng[1] - xrng[0]) / (grdshp[1]-1)
        yres = (yrng[1] - yrng[0]) / (grdshp[0]-1)
    else:
        xres = (xrng[1] - xrng[0]) / (grdshp[1])
        yres = (yrng[1] - yrng[0]) / (grdshp[0])

    # 2D interpolation gives poor results, so use mean of 1D in both directions
    # tests suggest average differences of about 0.1% of GMT grdtrack
    '''
    xbuff = xres * 5
    ybuff = yres * 5

    xind = where(logical_and(xdat > lon-xbuff, xdat < lon+xbuff) == True)
    yind = where(logical_and(ydat > lat-ybuff, ydat < lat+ybuff) == True)

    ztmp = zdat[yind[0][0]:yind[0][-1]+1,xind[0][0]:xind[0][-1]+1]

    zval = interp2d(ydat[yind], xdat[xind], ztmp, kind='linear')

    zint2d.append(zval(lon, lat)[0])
    '''

    xbuff = xres * 1
    ybuff = yres * 1
    zval = []

    # now loop thru lat/lon array
    for i in range(0, len(alon)):
        lat = alat[i]
        lon = alon[i]

        if pixelnode == False:
            criterion = (xdat >= lon-xbuff) & (xdat <= lon+xbuff*2)
            xind = where(criterion)[0]

            criterion = (ydat >= lat-ybuff) & (ydat <= lat+ybuff*2)
            yind = where(criterion)[0]
        else:
            criterion = (xdat >= lon-xbuff) & (xdat <= lon+xbuff)
            xind = where(criterion)[0]

            criterion = (ydat >= lat-ybuff) & (ydat <= lat+ybuff)
            yind = where(criterion)[0]

        try:

            # revise indicies if point on edge of map
            if yind[0] < 0:
            	yind[0] = yind[0][1:].append(yind[0][-1]+1)

            if yind[-1] >= grdshp[0] - 1:
            	yind[0] = sort(yind[0][0:-2].append(yind[0][0]-1))

            if xind[0] < 0:
            	xind[0] = xind[0][1:].append(xind[0][-1]+1)

            if xind[-1] >= grdshp[1] - 1:
            	xind[0] = sort(xind[0][0:-1].append(xind[0][0]-1))

            # make data arrays
            y = ydat[yind]
            for j in range(0,len(xind)-1):
                y = hstack((y, ydat[yind]))

            x = xdat[xind]
            for j in range(0,len(yind)-1):
                x = hstack((x, xdat[xind]))

            # if way outside return nan, otherwise continue

            atmp = zdat[yind[0]:yind[-1]+1,xind[0]:xind[-1]+1]
            z = atmp.data.flatten()

            # remove z=nan
            znan = ~isnan(z)
            z = z[znan]

            # get indices to sort
            xa = argsort(x[znan])
            ya = argsort(y[znan])

            # try 1d interp - cubic interpolation tends to crash
            zinterp = interp1d(x[xa], z[xa], kind='linear')
            xz = zinterp(lon)
            zinterp = interp1d(y[ya], z[ya], kind='linear')
            yz = zinterp(lat)

            zval.append(nanmean([xz, yz]))

        except:
            zval.append(nan)

    # delete tmp file if it exists
    if path.isfile('tmp.grd') == True:
        remove('tmp.grd')

    return zval
    
# get correct Slab1.0 grid based on eq location
def get_slab_gridnames(lat, lon):
    from numpy import nan
    # find correct region
    if lon < 0:
        lon = 360 + lon

    if lon >= 167 and lon <= 216 \
       and lat >= 50 and lat <= 65:
       dipfile = 'alu_slab1.0_dipclip.grd'
       stkfile = 'alu_slab1.0_strclip.grd'
       depfile = 'alu_slab1.0_clip.grd'
       avdip = 17.
       maxdep = 53.
       seiswid = 146.
    elif lon >= 136 and lon <= 148 \
       and lat >= 11 and lat <= 34.5:
       dipfile = 'izu_slab1.0_dipclip.grd'
       stkfile = 'izu_slab1.0_strclip.grd'
       depfile = 'izu_slab1.0_clip.grd'
       avdip = 16.
       maxdep = 34.
       seiswid = 85.
    elif lon >= 174 and lon <= 188 \
       and lat >= -39 and lat <= -14:
       dipfile = 'ker_slab1.0_dipclip.grd'
       stkfile = 'ker_slab1.0_strclip.grd'
       depfile = 'ker_slab1.0_clip.grd'
       avdip = 20.
       maxdep = 47.
       seiswid = 109.
    elif lon >= 129 and lon <= 164 \
       and lat >= 32 and lat <= 56.5:
       dipfile = 'kur_slab1.0_dipclip.grd'
       stkfile = 'kur_slab1.0_strclip.grd'
       depfile = 'kur_slab1.0_clip.grd'
       avdip = 17.
       maxdep = 55.
       seiswid = 153.
    elif lon >= 254 and lon <= 279 \
       and lat >= 7 and lat <= 21:
       dipfile = 'mex_slab1.0_dipclip.grd'
       stkfile = 'mex_slab1.0_strclip.grd'
       depfile = 'mex_slab1.0_clip.grd'
       avdip = 22.
       maxdep = 38.
       seiswid = 75.
    elif lon >= 122 and lon <= 128 \
       and lat >= 7 and lat <= 15:
       dipfile = 'phi_slab1.0_dipclip.grd'
       stkfile = 'phi_slab1.0_strclip.grd'
       depfile = 'phi_slab1.0_clip.grd'
       avdip = 26.
       maxdep = 46.
       seiswid = 79.
    elif lon >= 122 and lon <= 138 \
       and lat >= 22 and lat <= 36:
       dipfile = 'ryu_slab1.0_dipclip.grd'
       stkfile = 'ryu_slab1.0_strclip.grd'
       depfile = 'ryu_slab1.0_clip.grd'
       avdip = 18.
       maxdep = 46.
       seiswid = 116.
    elif lon >= 278 and lon <= 300 \
       and lat >= -45 and lat <= 5:
       dipfile = 'sam_slab1.0_dipclip.grd'
       stkfile = 'sam_slab1.0_strclip.grd'
       depfile = 'sam_slab1.0_clip.grd'
       avdip = 16.
       maxdep = 48.
       seiswid = 137.
    elif lon >= 328 and lon <= 337 \
       and lat >= -61 and lat <= 55:
       dipfile = 'sco_slab1.0_dipclip.grd'
       stkfile = 'sco_slab1.0_strclip.grd'
       depfile = 'sco_slab1.0_clip.grd'
       avdip = 17.
       maxdep = 35.
       seiswid = 68.
    elif lon >= 145 and lon <= 165 \
       and lat >= -12 and lat <= -2:
       dipfile = 'sol_slab1.0_dipclip.grd'
       stkfile = 'sol_slab1.0_strclip.grd'
       depfile = 'sol_slab1.0_clip.grd'
       avdip = 31.
       maxdep = 56.
       seiswid = 83.
    elif lon >= 91 and lon <= 124 \
       and lat >= -12 and lat <= 11:
       dipfile = 'sum_slab1.0_dipclip.grd'
       stkfile = 'sum_slab1.0_strclip.grd'
       depfile = 'sum_slab1.0_clip.grd'
       avdip = 15.
       maxdep = 51.
       seiswid = 154.
    elif lon >= 164 and lon <= 173 \
       and lat >= -23.5 and lat <= -9:
       dipfile = 'van_slab1.0_dipclip.grd'
       stkfile = 'van_slab1.0_strclip.grd'
       depfile = 'van_slab1.0_clip.grd'
       avdip = 25.
       maxdep = 42.
       seiswid = 70.
    elif lon >= 231.5 and lon <= 239.5 \
       and lat >= 39 and lat <= 52:
       dipfile = 'cas_slab1.0_dipclip.grd'
       stkfile = 'cas_slab1.0_strclip.grd'
       depfile = 'cas_slab1.0_clip.grd'
       avdip = nan
       maxdep = nan
       seiswid = nan
    else:
       dipfile = ''
       stkfile = ''
       depfile = ''
       avdip = nan
       maxdep = nan
       seiswid = nan

    return dipfile, stkfile, depfile, avdip, maxdep, seiswid

# parses GMT polygon file into python dictionary
def parse_gmt_polys(polyfile):
    from shapely.geometry import Polygon

    # set variables
    polys = []
    points = []
    finish = True

    # read file
    txt = open(polyfile).readlines()

    # strip leading '>'
    if txt[0].find('>') >= 0:
        txt = txt[1:]

    # now loop thru mask file
    for line in txt:
        if line[0].strip() != '>':
            dat = line.rstrip(' \t\n\r').split()
            tmppt = [float(dat[0]), float(dat[1])]
            points.append(tmppt)
            finish = True

        else:
            # close polygon if necessary
            if points[0][0] != points[-1][0] and points[0][1] != points[-1][1]:
                points.append(points[0])

            polys.append(Polygon(points))
            points = []
            finish = False

    # finish last polygon if no trailing line
    if finish == True:
        if points[0][0] != points[-1][0] and points[0][1] != points[-1][1]:
            points.append(points[0])

        polys.append(Polygon(points))

    return polys

# converts xyz file to point shapefile
def xyz2shp(xyzfile, headerlines, shpfile, shptype):
    import shapefile
    import re

    # now parse ceef
    print 'Reading xyz...'
    data = open(xyzfile).readlines()[headerlines:]

    # now output shapefile
    print 'Making shapefile...'
    if shptype == 'polyline':
        w = shapefile.Writer(shapefile.POLYLINE)
    else:
        w = shapefile.Writer(shapefile.POINT)
    w.field('LON','F', 13, 6)
    w.field('LAT','F', 13, 6)
    w.field('Z','F', 13, 6)

    # now loop through points
    for line in data:
        print line
        # split space, tab or comma delimitered
        #dat = re.split('\s+', line.strip())
        dat = line.strip().split(',')
        print dat
        if len(dat) != 3:
            dat = dat.split('\t')
        if len(dat) != 3:
            dat = dat.split(',')

        w.point(float(dat[0]), float(dat[1]))
        w.record(float(dat[0]), float(dat[1]), float(dat[2]))

    print 'Writing shapefile...'
    w.save(shpfile)

# converts grd file to point shapefile
def grd2shp(grdfile, shpfile):
    from os import system, remove
    try:
        # first output to tmp xyz file
        system('grd2xyz '+grdfile+' > tmp.xyz')
        print '\nWriting file: ', shpfile
        xyz2shp('tmp.xyz', shpfile)
        remove('tmp.xyz')
    except:
        print '\nGMT module not found\n'

# converts shp polygons or lines to GMT-friendly file
def shp2gmt(shpfile, outfile, **kwargs):
    '''
    kwargs decides what to to with the following header:
        header="none"
        header="num", field="field name"
        header="str", field="field name"
    '''

    import shapefile

    print 'Reading shapefile...'
    sf = shapefile.Reader(shpfile)
    shapes = sf.shapes()
    records = sf.records()
    nrec = len(records)

    # now make output text file
    # first overwirte
    f = open(outfile,'w')
    f.close()
    all_str = ''

    tmpfile = open(outfile,'a')
    dosimple = True

    # loop through polygons
    for k in range(0,nrec):
    # write polygon to temp file
    # if want to include name

    # if want to include quantity
    #all_str = '> -Z' + str(norm_complete[k,0]) + '\n'

        # check to see if shape has multiple parts
        p = 0
        if len(shapes[k].points) != 0:
            parts = shapes[k].parts
            parts.append(len(shapes[k].points)-1)
            for part in range(0,len(parts)-1):
                pt_str = ''
                all_str = ''
                if dosimple == True:
                    all_str = '>' + '\n'

                while p <= parts[part+1]:
                    pt_str = pt_str + str("%0.5f" % shapes[k].points[::-1][p][0])+"\t" \
                         +str("%0.5f" % shapes[k].points[::-1][p][1])+"\n"    # kluge to reverse points
                    p += 1

                all_str = all_str + pt_str
                tmpfile.write(all_str)

    print 'Writing to file...'

    tmpfile.close()


# converts shp points to psxy friently files with labels
def shp_pt2gmt_pt(shpfile, outfile, labelfield):
    '''
    kwargs decides what to to with the following header:
        header="none"
        header="num", field="field name"
        header="str", field="field name"
    '''

    import shapefile
    from mapping_tools import get_field_data

    print 'Reading shapefile...'
    sf = shapefile.Reader(shpfile)

    # get labels
    labels = get_field_data(sf, labelfield, 'str') # assume always use string

    shapes = sf.shapes()
    records = sf.records()

    all_str = ''

    # loop through points
    for i, rec in enumerate(records):
        all_str = all_str + '\t'.join((str("%0.5f" % shapes[i].points[0][0]), \
                  str("%0.5f" % shapes[i].points[0][1]), '11', '0', '1', 'ML', \
                  labels[i])) + '\n'

    print 'Writing to psxy file...'
    f = open(outfile,'wb')
    f.write(all_str)
    f.close()

'''
code below stolen from:
http://wiki.scipy.org/Cookbook/Matplotlib/Loading_a_colormap_dynamically
'''

def cpt2colormap(fileName, ncolours, **kwargs):

    import colorsys
    from numpy import array, interp, linspace
    from pylab import matplotlib

    # get kwargs
    rev = False
    for key in ['rev']:
        if key in kwargs:
            # set fault type
            if key == 'rev':
                rev = kwargs[key]

    try:
        f = open(fileName)
    except:
        print "file ",fileName, "not found"
        return None

    lines = f.readlines()
    f.close()

    x = []
    r = []
    g = []
    b = []
    colorModel = "RGB"
    for l in lines:
        ls = l.split()
        if l[0] == "#":
           if ls[-1] == "HSV":
               colorModel = "HSV"
               continue
           else:
               continue
        if ls[0] == "B" or ls[0] == "F" or ls[0] == "N":
           pass
        else:
            x.append(float(ls[0]))
            r.append(float(ls[1]))
            g.append(float(ls[2]))
            b.append(float(ls[3]))
            xtemp = float(ls[4])
            rtemp = float(ls[5])
            gtemp = float(ls[6])
            btemp = float(ls[7])

    x.append(xtemp)
    r.append(rtemp)
    g.append(gtemp)
    b.append(btemp)

#    nTable = len(r)
    x = array( x )
    r = array( r )
    g = array( g )
    b = array( b )
    if colorModel == "HSV":
       for i in range(r.shape[0]):
           rr,gg,bb = colorsys.hsv_to_rgb(r[i]/360.,g[i],b[i])
           r[i] = rr ; g[i] = gg ; b[i] = bb
    if colorModel == "HSV":
       for i in range(r.shape[0]):
           rr,gg,bb = colorsys.hsv_to_rgb(r[i]/360.,g[i],b[i])
           r[i] = rr ; g[i] = gg ; b[i] = bb
    if colorModel == "RGB":
        r = r/255.
        g = g/255.
        b = b/255.

    # reverse order
    if rev == True:
        r = r[::-1]
        g = g[::-1]
        b = b[::-1]

    # interpolate to ncolours
    xx = linspace(x[0], x[-1], ncolours)
    r = interp(xx, x, r)
    g = interp(xx, x, g)
    b = interp(xx, x, b)
    x = xx

    xNorm = (x - x[0])/(x[-1] - x[0])

    red = []
    blue = []
    green = []
    for i in range(len(x)):
        red.append([xNorm[i],r[i],r[i]])
        green.append([xNorm[i],g[i],g[i]])
        blue.append([xNorm[i],b[i],b[i]])

    colorDict = {"red":red, "green":green, "blue":blue}

    return matplotlib.colors.LinearSegmentedColormap('my_colormap',colorDict,ncolours), xx
    
def makecpt (zvals, colmat, cptout):
    from numpy import array
    """
    #colmat = 3 x n matrix
    colmat = [[255, 255, 229],
              [255, 247, 188],
              [254, 227, 145],
              [254, 196, 79],
              [254, 153, 41],
              [236, 112, 20],
              [204, 76, 2],
              [140, 45, 4]]

    colmat = array(colmat)/256.
    """
    
    cptstr = ''
    for i, c in enumerate(colmat):
        cs = (str(x) for x in c)
        csj = '\t'.join(cs)   
        cptstr += '\t'.join((str(zvals[i]), csj, str(zvals[i+1]), csj)) + '\n'
        
    # write to file
    f = open(cptout, 'wb')
    f.write(cptstr)
    f.close()

def make_netcdf_map(ax, cnrs, ncfile, cmap, norm, vmin, vmax, mapres, grdres, grdsize, zscale, **kwargs):
    from numpy import arange, percentile, mean
    from netCDF4 import Dataset as NetCDFFile
    from mpl_toolkits.basemap import Basemap
    from mapping_tools import mask_outside_polygons, get_map_polygons
    import matplotlib.pyplot as plt
    
    #get kwargs
    lightsource = False
    filloceans = False

    for key in ('lightsource', 'filloceans'):
        if key in kwargs:
            # min value
            if key == 'lightsource':
                lightsource = kwargs[key]

            # set scaling relation
            if key == 'filloceans':
                filloceans = kwargs[key]

    
    llcrnrlon = cnrs[0]
    urcrnrlon = cnrs[1]
    llcrnrlat = cnrs[2]
    urcrnrlat = cnrs[3]
    
    lon_0 = mean([llcrnrlon, urcrnrlon])
    lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
    lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

    plt.tick_params(labelsize=12)
    
    m = Basemap(projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,\
                llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
                urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
                rsphere=6371200.,resolution=mapres,area_thresh=300)

    # draw coastlines, state and country boundaries, edge of map.
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    #m.drawmapboundary(fill_color='0.8', zorder=100)
    m.drawparallels(arange(-90.,90.,grdsize), labels=[1,0,0,0],fontsize=12, dashes=[2, 2], color='0.5', linewidth=0.75)
    m.drawmeridians(arange(0.,360.,grdsize), labels=[0,0,0,1], fontsize=12, dashes=[2, 2], color='0.5', linewidth=0.75)
    
    
    ##########################################################################################
    # plot grids
    ##########################################################################################
    
    print 'Reading netCDF file...'
    print ncfile
    nc = NetCDFFile(ncfile)
    
    data = nc.variables['z'][:] #/ zscale
    try:
        lons = nc.variables['lon'][:]
        lats = nc.variables['lat'][:]
    except:
        lons = nc.variables['x'][:]
        lats = nc.variables['y'][:]
    
    # transform to metres
    nx = int((m.xmax-m.xmin)/grdres)+1
    ny = int((m.ymax-m.ymin)/grdres)+1
    
    topodat = m.transform_scalar(data,lons,lats,nx,ny)
    
    
    # make shading
    print 'Making map...'
    if lightsource == True:
        from matplotlib.colors import LightSource
        ls = LightSource(azdeg = 180, altdeg = 45)
        rgb = ls.shade(topodat, cmap=cmap, norm=norm)
        
        im = m.imshow(rgb)   
    else:
        im = m.imshow(topodat, cmap=cmap, norm=norm, vmin=vmin, vmax=vmax, interpolation='nearest')
    
    ##########################################################################################
    # get land & lake polygons for masking
    ##########################################################################################
    if filloceans == True:
        polys = get_map_polygons(m)
        
        #mask_outside_polygon(polys[1][::-1], ax=None)
        mask_outside_polygons(polys, 'lightskyblue', plt)
        
        # get lake ploygons
        polygons = []
        for polygon in m.lakepolygons:
            poly = polygon.get_coords()
            plt.fill(poly[:,0], poly[:,1], 'lightskyblue')
            polygons.append(poly)
        
    return m
