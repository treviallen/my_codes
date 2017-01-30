
# gets distance and az along depth contour from previous point of origin
def get_dist_az_vect(lastlat, lastlon, latvect, lonvect):
   from numpy import array
   from mapping_tools import distance

   rng = []
   az  = []
   baz = []

   # assume opposite dirn ~180 degrees
   for l in range(0,len(latvect)):
       rngtmp, aztmp, baztmp = distance(lastlat, lastlon, latvect[l], lonvect[l])
       rng.append(rngtmp)
       az.append(aztmp)
       baz.append(baztmp)

   return array(rng), array(az), array(baz)


def get_next_point(lat, lon, latvect, lonvect, res):
    from numpy import argwhere
    from mapping_tools import reckon

    rng, az, baz = get_dist_az_vect(lat, lon, latvect, lonvect)

    # if points < res, remove
    ltind = argwhere(rng < res)
    if len(ltind) > 0:
        latvect = latvect[ltind[-1]+1:]
        lonvect = lonvect[ltind[-1]+1:]
        rng     = rng[ltind[-1]+1:]

    if rng[0] < rng[-1]:
        return latvect, lonvect, lonvect[0], latvect[0], rng[0]
    else:
        return latvect, lonvect, lonvect[-1], latvect[-1], rng[-1]


def get_pts_along_wid(flon, flat, fdep, lastlon, lastlat, wid, ncdip, ncstk, ncdep, res, maxdep):
	from gmt_tools import pygrdtk
	from numpy import cos, radians, isnan
	from mapping_tools import reckon

	wlastlon1 = lastlon
	wlastlat1 = lastlat
	wlastlon2 = lastlon
	wlastlat2 = lastlat
	totwid = 0.0
	doW1 = True
	doW2 = True

	# get dep of origin
	#print ncdip, ncstk, ncdep
	#fdep.append(pygrdtk(lastlon, lastlat, nc=ncdep)[0])
	stopwid = False

	while totwid < wid and stopwid == False:
		# stop getting width if cannot do any more
		if doW1 == False and doW2 == False:
			stopwid = True

		if doW1 == True:
			try:
				# get dip & stk at current point
				dip = pygrdtk(wlastlon1, wlastlat1, nc=ncdip)[0]
				stk = pygrdtk(wlastlon1, wlastlat1, nc=ncstk)[0]

				# get horiz distance
				hwid = abs(res * cos(radians(dip)))

				# get down dip points
				lonlat = reckon(wlastlat1, wlastlon1, hwid, stk+90)
				tmpdep = pygrdtk(lonlat[0], lonlat[1], nc=ncdep)[0]


				# depth cannot exceed 75th percentile depth from BOR dataset or from Hayes etal 2012
				if isnan(maxdep):
					maxdep = -52
				else:
					maxdep = -1*maxdep

				if tmpdep > maxdep and isnan(tmpdep) == False:
					flon.append(lonlat[0])
					flat.append(lonlat[1])
					fdep.append(tmpdep)
					wlastlon1 = lonlat[0]
					wlastlat1 = lonlat[1]
					totwid += res
				else:
					doW1 = False
			except:
				doW1 = False

		if doW2 == True and totwid < wid:
			try:
				dip = pygrdtk(wlastlon2, wlastlat2, nc=ncdip)[0]
				stk = pygrdtk(wlastlon2, wlastlat2, nc=ncstk)[0]

				# get horiz distance
				hwid = abs(res * cos(radians(dip)))

				# get up dip points
				wlastlat2, wlastlon2, hwid, stk-90
				lonlat = reckon(wlastlat2, wlastlon2, hwid, stk-90)

				tmpdep = pygrdtk(lonlat[0], lonlat[1], nc=ncdep)[0]

				if isnan(tmpdep) == False:
					#print 'W2', tmpdep
					flon.append(lonlat[0])
					flat.append(lonlat[1])
					fdep.append(tmpdep)
					wlastlon2 = lonlat[0]
					wlastlat2 = lonlat[1]
					totwid += res
				else:
					doW2 = False

			except:
				doW2 = False
				#print tmpdep,'False2'

	return flon, flat, fdep

# get correct grid based on eq location
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

'''
Start main
'''

def plot_faults(evlat, evlon, evdep, flat, flon, fdep):
	#import matplotlib.pyplot as plt
	#plt.plot(flon, flat, 'g-')
	from mpl_toolkits.mplot3d import axes3d, Axes3D
	from matplotlib.mlab import griddata
	import matplotlib.pyplot as plt
	fig = plt.figure()
	plt.plot(evlon, evlat, 'ko')
	ax = Axes3D(fig)
	ax.view_init(30, 30)
	ax.scatter(evlon, evlat,evdep, 'bo')
	ax.scatter(flon, flat, fdep, marker='+', color = 'r')
	# define grid
	#step = res/111.
	#xi = arange(min(flon), max(flon)+step, step)
	#yi = arange(min(flat), max(flat)+step, step)
	#ZI = griddata(flon,flat,fdep,xi,yi,interp='linear')
	#XI, YI = meshgrid(xi, yi)
	#
	#ax.plot_surface(XI, YI, ZI)
	plt.show()


'''
Second try
'''
def make_slab_matrix(lat, lon, mw, res, grdpath):
    from gmt_tools import pygrdtk
    from mapping_tools import get_line_parallels, distance, reckon
    from netCDF4 import Dataset as NetCDFFile
    from numpy import meshgrid, argsort, argwhere, argmin, logical_and, logical_or, abs, radians, cos, arange, nan, hstack
    import scipy.ndimage as ndimage
    from os import path
    import matplotlib.pyplot as plt
    from make_slab_fault import get_slab_gridnames

    import warnings
    warnings.filterwarnings("ignore")

    # testing data
    '''
    lat = -5.
    lon = -78
    mw = 8.5
    res = 10
    lat = 16.012
    lon = -96
    mw = 7.8
    res = 10
    grdpath = '/Users/trev/Documents/Data/Slab_1.0'
    '''

    #try:
    # set grids
    dipfile, stkfile, depfile, avdip, maxdep, seiswid = get_slab_gridnames(lat, lon)
    print 'maxdep', maxdep
    '''
    if lon < 0 and dipfile[0:3] != 'cas':
        lon = 360 + lon
    '''
    if lon < 0 and dipfile[0:3] != 'cas':
        lon = 360 + lon

    dipfile = path.join(grdpath, dipfile)
    stkfile = path.join(grdpath, stkfile)
    depfile = path.join(grdpath, depfile)

    #read grids
    ncdep = NetCDFFile(depfile)
    zdat = ncdep.variables['z'][:]
    xdat = ncdep.variables['x'][:]
    ydat = ncdep.variables['y'][:]

    # read dipgrd
    ncdip = NetCDFFile(dipfile)

    # read stkgrd
    ncstk = NetCDFFile(stkfile)

    dep = pygrdtk(lon, lat, nc=ncdep)
    dip = pygrdtk(lon, lat, nc=ncdip)
    stk = pygrdtk(lon, lat, nc=ncstk)

    # set hypo
    if lon > 180:
        hlon = lon - 360
    else:
        hlon = lon
    hypo = [hlon, lat, dep[0]]

    # resample grids
    '''
    zdat = ndimage.zoom(zdat, 0.2, order=1)
    xdat = ndimage.zoom(xdat, 0.2, order=1)
    ydat = ndimage.zoom(ydat, 0.2, order=1)
    '''

    X, Y = meshgrid(xdat, ydat)

    # get allen len for given mag - updated 2015-04-15
    a = -2.623
    b =  0.6103
    sig = 0.277

    flen  = 10**(a + b * mw)

    # get wid from allen - updated 2015-04-15
    a   = -2.247
    b   = 0.539
    c   = 2.318
    fwid = 10**(min([a + b*mw, c]))

    print 'dims',flen,fwid

    # smooth data before contouring
    #zdatsmooth = ndimage.gaussian_filter(zdat, sigma=.5, order=0)

    # contour at epicentre - slab interface
    cs = plt.contour(X, Y, zdat, dep, color='g')
    #plt.plot(lon, lat, 'ro')
    #plt.show()
    plt.clf()

    lonvect = []
    latvect = []
    for i in range(0, len(cs.collections[0].get_paths())):
        p = cs.collections[0].get_paths()[i]
        v = p.vertices
        lonvect = hstack((lonvect, v[:,0]))
        latvect = hstack((latvect, v[:,1]))

    # discretize lat and lon vectors
    '''
    if len (lonvect) > 25:
        lonvect = lonvect[range(0, len(lonvect), 25)]
        latvect = latvect[range(0, len(latvect), 25)]
    '''
    # get last pos
    try:
        lastlat1 = lat[0]
    except:
        lat = [lat]
        lon = [lon]
        lastlat1 = lat[0]

    lastlon1 = lon[0]
    lastlat2 = lat[0]
    lastlon2 = lon[0]

    # now step out for length of fault
    # assume azimths in each direction don't differ by more then +/- 90deg
    flat = []
    flon = []
    fdep = []
    flat.append(lat[0])
    flon.append(lon[0])
    fdep.append(dep[0])

    # now get points along the width for origin
    flon, flat, fdep = get_pts_along_wid(flon, flat, fdep, lon[0], lat[0], fwid, ncdip, ncstk, ncdep, res, maxdep)
    #print flon[1], flat[1], lastlat1, lastlon1, '1st'

    # for contour, split into two arrays
    rng, az, baz = get_dist_az_vect(lat[0], lon[0], latvect, lonvect)

    # get min dist index
    if len(rng) > 0:
        mi = argmin(rng)

        # split into two lat/lon arrays
        latd1 = latvect[0:mi]
        lond1 = lonvect[0:mi]
        rngd1 = rng[0:mi]

        latd2 = latvect[mi+1:]
        lond2 = lonvect[mi+1:]
        rngd2 = rng[mi+1:]

        # sort by distance
        rngind = argsort(rngd1)
        rngd1  = rngd1[rngind]
        latd1  = latd1[rngind]
        lond1  = lond1[rngind]

        rngind = argsort(rngd2)
        rngd2  = rngd2[rngind]
        latd2  = latd2[rngind]
        lond2  = lond2[rngind]

        # now loop through contour points until equal fault len
        i = 0
        totlen = 0
        nxtlon1 = lon[0]
        nxtlat1 = lat[0]
        nxtlon2 = lon[0]
        nxtlat2 = lat[0]
        doD1 = True
        doD2 = True
        stoploop = False

        while totlen < flen and stoploop == False:
            if doD1 == True:
                try:
                    # get next point in dirn1
                    #print totlen, nxtlat1, nxtlon1, len(latd1), latd1, lond1, 'L1'
                    latd1, lond1, nxtlon1, nxtlat1, dres = get_next_point(nxtlat1, nxtlon1, latd1, lond1, res)
                    flon.append(nxtlon1)
                    flat.append(nxtlat1)
                    fdep.append(pygrdtk(nxtlon1, nxtlat1, nc=ncdep)[0])
                    totlen += dres
                    print totlen, nxtlat1, len(latd1), 'L1'

                    flon, flat, fdep = get_pts_along_wid(flon, flat, fdep, nxtlon1, nxtlat1, fwid, ncdip, ncstk, ncdep, res, maxdep)

                except:
                    doD1 = False


            if doD2 == True and totlen < flen:
                try:
                    # get next point in dirn 2
                    latd2, lond2, nxtlon2, nxtlat2, dres = get_next_point(nxtlat2, nxtlon2, latd2, lond2, res)
                    flon.append(nxtlon2)
                    flat.append(nxtlat2)
                    fdep.append(pygrdtk(nxtlon2, nxtlat2, nc=ncdep)[0])
                    totlen += dres
                    print totlen, nxtlat2, len(latd2), 'L2'

                    flon, flat, fdep = get_pts_along_wid(flon, flat, fdep, nxtlon2, nxtlat2, fwid, ncdip, ncstk, ncdep, res, maxdep)

                except:
                    doD2 = False

            if doD1 == False and doD2 == False:
                stoploop = True

        return flat, flon, fdep, hypo

