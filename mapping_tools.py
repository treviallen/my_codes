# -*- coding: utf-8 -*-
"""
Created on Tue Jun 04 14:55:21 2013

@author: tallen
"""
def print_functions():
    txt = open('U:/Code/pycode/mapping_tools.py').readlines()
    for line in txt:
        if line.find('def') >= 0:
            print(line.strip('def ').strip('\n'))
            

def load_shape(shpfile):
    import shapefile
    return shapefile.Reader(shpfile)
    
def get_shp_polygons(shpfile):
    import shapefile
    from shapely.geometry import Polygon
    
    sf = shapefile.Reader(shpfile)
    shapes = sf.shapes()
    polygons = []
    for poly in shapes:
        polygons.append(Polygon(poly.points))

    return sf, polygons

def get_field_index(sf, field):
    fields = sf.fields[1:]
    findex = -1
    for i, f in enumerate(fields):
        if f[0] == field:
            findex = i
            
    return findex
    
def get_field_data(sf, field, datatype):
    '''
    datatype = float, str
    '''
    from misc_tools import checkfloat, checkint    
    from numpy import array
    
    # get index
    index = get_field_index(sf, field)
    
    # get records
    recs = sf.records()
    
    # now loop thru recs and get data
    data = []
    for rec in recs:
        if datatype == 'str':
            data.append(rec[index])
        elif datatype == 'float':
            data.append(checkfloat(rec[index]))
        elif datatype == 'int':
            data.append(checkint(rec[index]))
            
    return array(data)
        
def shapepoly2shapley(shpfile):
    import shapefile
    from shapely.geometry import Polygon    
    
    sf = shapefile.Reader(shpfile)
    shapes = sf.shapes()
    polygons = []
    for poly in shapes:
        polygons.append(Polygon(poly.points))
        
    return polygons

# returns list of shapely polygons of different shape parts    
def split_shape_parts(sf):
    from shapely.geometry import Polygon
    from numpy import array
    
    shapes = sf.shapes()
    polygons = []

    for i, shape in enumerate(shapes):
        parts = shape.parts
        print(parts)
        parts.append(len(shape.points)-1) # adding last point
        
        for i in range (0, len(parts)-1):
            polygons.append(Polygon(shape.points[parts[i]:parts[i+1]]))
            print(array(polygons[-1].bounds))
    return polygons

'''
# get colour map and vector of colour indexes
    sf     = shapefile object
    field  = relevant field in shapefile
    colmap = colourmap to be used
    ncolours = number of colour ranges
kwargs:
    zmin
    zmax    
'''
def getshapecolour(sf, field, colmap, ncolours, **kwargs):
    import matplotlib.cm as cm
    from numpy import arange, isnan, nan
    
    # if assigned, set kwargs
    zmin = nan
    zmax = nan
    for key in ('zmin', 'zmax'):
        if key in kwargs:
            # min value
            if key == 'zmin':
                zmin = kwargs[key]

            # set scaling relation
            if key == 'zmax':
                zmax = kwargs[key]
    
    # get field index
    findex = get_field_index(sf, field)
    
    # get max/min value in field
    fvect = []
    recs = sf.records()
    for rec in recs:
        fvect.append(float(rec[findex]))
    
    if isnan(zmin):    
        zmin = min(fvect)
    if isnan(zmax):
        zmax = max(fvect)
    zrng = zmax - zmin
    
    # make colourmap
    cmap = cm.get_cmap(colmap, ncolours)
    cs = (cmap(arange(ncolours)))
    
    # get colour index
    ci = []
    for f in fvect:
        print(f)
        ci.append(round((ncolours-1.) / zrng * (f-zmin)))
        
    return cs, ci, cmap, zmin, zmax
    
def drawshapepoly(m, plt, sf, label='null', fillcolor='none', edgecolor='k', alpha=1, cmap=-99, zorder=20000, **kwargs):
    from numpy import arange, isnan, nan
    
    # get kwargs
    ncolours = 256
    lw = 0.5
    ls = '-'
    fillshape = False
    polyline = False # do not close polygon

    for key in ('col', 'cindex', 'cmap', 'ncolours', 'lw', 'ls', 'fillshape', 'polyline'):
        if key in kwargs:
            if key == 'col':
                col = kwargs[key]
                
            if key == 'cindex':
                cindex = kwargs[key]

            if key == 'cmap':
                cmap = kwargs[key]

            if key == 'ncolours':
                ncolours = kwargs[key]

            if key == 'ls':
                ls = kwargs[key]

            if key == 'lw':
                lw = kwargs[key]

            if key == 'fillshape':
                fillshape = kwargs[key]
                 
            if key == 'polyline':
                polyline = kwargs[key]
    
    shapes = sf.shapes()

    for i, shape in enumerate(shapes):
        newfill = True
        # get colour
        try:
            cs = (cmap(arange(ncolours)))
            if isnan(cindex[i]) == False and cindex[i] != -1: 
                col = [cs[int(cindex[i])][0],cs[int(cindex[i])][1],cs[int(cindex[i])][2]]
            else:
                newfill = False
                col = '0.9'
                #col = 'none'
            
        except:
            try:
                col = col
            except:
                col = 'w'
                
        '''
        if fillshape == True:
            fillcol = col
            linecol = 'k'
        else:
            linecol = col
        '''
        # get xy coords
        x = []
        y = []

        p = 0
        parts = shape.parts
        parts.append(len(shape.points))
        
        for prt in parts[1:]:
            #print('part',prt
            while p < prt:
                x.append(shape.points[p][0])
                y.append(shape.points[p][1])
                p += 1
            
            # close polygon if not polyline
            if polyline == False:
                if x[0] != x[-1] or y[0] != y[-1]:
                    x.append(x[0])
                    y.append(y[0])

            try:
              # plot each polygon
                
                xx, yy = m(x,y)
                #print(edgecolor)
                if fillshape == True and cmap != -99:
                    plt.fill(xx,yy, facecolor=col, edgecolor=edgecolor, linewidth=lw, label=label, alpha=alpha)  
                    #m.plot(xx, yy, linewidth=lw, color=edgecolor, ls=ls, zorder=1)
                elif fillshape == True and newfill == True:
                    #print(fillcolor)
                    if label == 'null':
                        plt.fill(xx,yy, facecolor=fillcolor, edgecolor=edgecolor, linewidth=lw, alpha=alpha)
                        m.plot(xx, yy, linewidth=lw, color=edgecolor, ls=ls, zorder=1)
                    else:
                        plt.fill(xx,yy, facecolor=fillcolor, edgecolor=edgecolor, linewidth=lw, label=label, alpha=alpha)  
                        m.plot(xx, yy, linewidth=lw, color=edgecolor, ls=ls, zorder=1)
                        label = 'null' # set to null after first plot
                elif label == 'null':
                    m.plot(xx, yy, linewidth=lw, color=edgecolor, linestyle=ls, zorder=zorder)
                else:
                    m.plot(xx, yy, linewidth=lw, color=edgecolor, linestyle=ls, label=label, zorder=zorder)
                    label = 'null' # set to null after first plot
            except:
                print('Skipping polygon...')

            x = []
            y = []
            
    '''
    end of drawshapepoly
    '''
    
def drawoneshapepoly(m, plt, sf, field, value, **kwargs):
    from numpy import arange
    
    # get kwargs
    ncolours = 256
    col = 'k'
    lw = 1.0
    ls = '-'
    fillshape = False
    polyline = False # do not close polygon

    for key in ('col', 'cindex', 'cmap', 'ncolours', 'lw', 'ls', 'fillshape', 'polyline'):
        if key in kwargs:
            if key == 'col':
                col = kwargs[key]
                
            if key == 'cindex':
                cindex = kwargs[key]

            if key == 'cmap':
                cmap = kwargs[key]

            if key == 'ncolours':
                ncolours = kwargs[key]

            if key == 'ls':
                ls = kwargs[key]

            if key == 'lw':
                lw = kwargs[key]

            if key == 'fillshape':
                fillshape = kwargs[key]
                 
            if key == 'polyline':
                polyline = kwargs[key]


    shapes = sf.shapes()
    recs = sf.records()
    
    
    # get field index
    findex = get_field_index(sf, field)
    #print('findex', findex)

    for i, shape in enumerate(shapes):
        #print('i', i)
        if recs[i][findex] == value:
            # get colour
            try:
                cs = (cmap(arange(ncolours)))
                col = [cs[int(cindex[i])][0],cs[int(cindex[i])][1],cs[int(cindex[i])][2]]
            except:
                try:
                    col = col
                except:
                    col = 'k'
    
            if fillshape == True:
                fillcol = col
                linecol = 'k'
            else:
                linecol = col
    
            # get xy coords
            x = []
            y = []
    
            p = 0
            parts = shape.parts
            parts.append(len(shape.points)-1)
            
            for prt in range(0,len(parts)-1):
                while p <= parts[prt+1]:
                    x.append(shape.points[p][0])
                    y.append(shape.points[p][1])
                    p += 1
                
                # close polygon if not polyline
                if polyline == False:
                    if x[0] != x[-1] or y[0] != y[-1]:
                        x.append(x[0])
                        y.append(y[0])
    
                # plot each polygon
                xx, yy = m(x,y)
                if fillshape == True:
                    plt.fill(xx,yy,color=fillcol)
                m.plot(xx, yy, linewidth=lw, color=linecol, linestyle=ls, zorder=1)
    
                x = []
                y = []
            
    '''
    end of drawshapepoly
    '''
    
def drawmanualshapepoly(m, plt, sf, label='null', fillcolor='none', edgecolor='k', alpha=1, cmap=-99, zorder=20000, **kwargs):
    from numpy import arange, isnan, nan
    
    # get kwargs
    ncolours = 256
    lw = 0.5
    ls = '-'
    fillshape = False
    polyline = False # do not close polygon

    for key in ('col', 'cindex', 'cmap', 'ncolours', 'lw', 'ls', 'fillshape', 'polyline'):
        if key in kwargs:
            if key == 'col':
                col = kwargs[key]
                
            if key == 'cindex':
                cindex = kwargs[key]

            if key == 'cmap':
                cmap = kwargs[key]

            if key == 'ncolours':
                ncolours = kwargs[key]

            if key == 'ls':
                ls = kwargs[key]

            if key == 'lw':
                lw = kwargs[key]

            if key == 'fillshape':
                fillshape = kwargs[key]
                 
            if key == 'polyline':
                polyline = kwargs[key]
    
    shapes = sf.shapes()

    for i, shape in enumerate(shapes):
        newfill = True
        # get colour
        try:
            cs = (cmap(arange(ncolours)))
            if isnan(cindex[i]) == False and cindex[i] != -1: 
                col = [cs[int(cindex[i])][0],cs[int(cindex[i])][1],cs[int(cindex[i])][2]]
            else:
                newfill = False
                col = '0.9'
                #col = 'none'
            
        except:
            try:
                col = col
            except:
                col = 'w'
                
        '''
        if fillshape == True:
            fillcol = col
            linecol = 'k'
        else:
            linecol = col
        '''
        # get xy coords
        x = []
        y = []

        p = 0
        parts = shape.parts
        parts.append(len(shape.points))
        
        for prt in parts[1:]:
            #print('part',prt
            while p < prt:
                x.append(shape.points[p][0])
                y.append(shape.points[p][1])
                p += 1
            
            # close polygon if not polyline
            if polyline == False:
                if x[0] != x[-1] or y[0] != y[-1]:
                    x.append(x[0])
                    y.append(y[0])

            try:
              # plot each polygon
                
                xx, yy = m(x,y)
                #print(edgecolor)
                if fillshape == True :
                    plt.fill(xx,yy, facecolor=col, edgecolor=edgecolor, linewidth=lw, label=label, alpha=alpha)  
                    #m.plot(xx, yy, linewidth=lw, color=edgecolor, ls=ls, zorder=1)
                elif fillshape == True and newfill == True:
                    #print(fillcolor)
                    if label == 'null':
                        plt.fill(xx,yy, facecolor=fillcolor, edgecolor=edgecolor, linewidth=lw, alpha=alpha)
                        m.plot(xx, yy, linewidth=lw, color=edgecolor, ls=ls, zorder=1)
                    else:
                        plt.fill(xx,yy, facecolor=fillcolor, edgecolor=edgecolor, linewidth=lw, label=label, alpha=alpha)  
                        m.plot(xx, yy, linewidth=lw, color=edgecolor, ls=ls, zorder=1)
                        label = 'null' # set to null after first plot
                
            except:
                print('Skipping polygon...')

            x = []
            y = []
            
    '''
    end of drawshapepoly
    '''
        
# label polygon with a field in shapefile
def labelpolygon(m, plt, sf, field, addOutline=False, **kwargs):
    '''
    fstyle  = ['normal', 'italic', 'oblique']
    fweight = ['light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black']
    '''
    
    from shapely.geometry import Polygon
    from mapping_tools import get_field_index
    import matplotlib as mpl
    import matplotlib.patheffects as path_effects
    
    mpl.rcParams['pdf.fonttype'] = 42
    
    xoff = 0.
    yoff = 0.
    fsize = 10
    col = 'k'
    fstyle = 'normal'
    fweight = 'normal'
    value = None
    
    for key in ('xoff', 'yoff', 'fsize', 'fweight', 'col', 'fstyle', 'value'):
        if key in kwargs:
            if key == 'xoff':
                xoff = kwargs[key]
            if key == 'yoff':
                yoff = kwargs[key]
            if key == 'fsize':
                fsize = kwargs[key]
            if key == 'fweight':
                fweight = kwargs[key]
            if key == 'fstyle':
                fstyle = kwargs[key]
            if key == 'col':
                col = kwargs[key]
            if key == 'value':
                value = kwargs[key]
            
    
    shapes = sf.shapes()
    recs = sf.records()
    
    # get field index
    findex = get_field_index(sf, field) 
    
    for i, shape in enumerate(shapes):
        centroid = Polygon(shape.points).centroid.wkt
        centroid = centroid.strip('PIONT').replace('(',' ').replace(')',' ').split()
        tx, ty = m(float(centroid[0]),float(centroid[1]))
        if tx > m.xmin and tx < m.xmax and ty > m.ymin and ty < m.ymax:
            if value == None or value == recs[i][findex]:
                
                if addOutline == False:
                    #path_effect=[path_effects.Normal()]
                    path_effect=[]
                else:
                    lineWidth= 3.
                    backColour = 'w'
                    path_effect=[path_effects.withStroke(linewidth=lineWidth, foreground=backColour)]    
                    
                plt.text(tx + xoff, ty + yoff, recs[i][findex], size=fsize, \
                         weight=fweight, color=col, style=fstyle, va='center', \
                         ha='center', path_effects=path_effect)
                     
        '''
        centroidx = []
            centroidy = []
            for j in range(0,len(shape.points)):
                centroidx.append(shape.points[j][0])
                centroidy.append(shape.points[j][1])
            tx, ty = m(mean(centroidx),mean(centroidy))
        '''

def return_map_shape_points(m, shpfile):
    import shapefile
    sf = shapefile.Reader(shpfile)
    shapes = sf.shapes()
    shplon = []
    shplat = []
    for shape in shapes:
        shplon.append(shape.points[0][0])
        shplat.append(shape.points[0][1])
    
    x, y = m(shplon, shplat)
    
    return x, y

# Add background to text to highlight
def addTextOutline(textHandle, lineWidth, backColour):
    '''
    e.g.:
        textHandle = plt.text(2,2,'This is a test', size=11, color='black')
    '''

    import matplotlib.patheffects as PathEffects
    textHandle.set_path_effects([PathEffects.withStroke(linewidth=lineWidth, foreground=backColour)])


def labelCentroid(plt, m, txt, fsize, xarray, yarray, xoff):
    from shapely.geometry import Polygon
    
    col = 'k'
    
    # make pts array
    pts = []
    for x, y in zip(xarray, yarray):
        pts.append([x+xoff, y])
    
    centroid = Polygon(pts).centroid.wkt
    centroid = centroid.strip('PIONT').replace('(',' ').replace(')',' ').split()
    tx, ty = m(float(centroid[0]),float(centroid[1]))
    plt.text(tx, ty, txt, size=fsize, color=col, va='center', ha='center')

    
# join two epicentres with a line
def joinpoints(m,plt,joinfile):
    data = open(joinfile).readlines()
    
    for line in data:
        pts = map(float, line.split(','))
        
        # draw connecting line
        xx, yy = m([pts[0],pts[2]],[pts[1],pts[3]])
        m.plot(xx,yy,'k-',linewidth=0.5)
        
        # plot points
        x1, y1 = m(pts[0],pts[1])
        x2, y2 = m(pts[2],pts[3])
        p1 = m.plot(x1,y1,'o',color='g',markersize=5) # SHEEF
        p2 = m.plot(x2,y2,'o',color='orange',ms=5) # WHITE
    
    return p1, p2

# map single points        
def mappoints(m,plt,ptfile):
    data = open(ptfile).readlines()
    
    for line in data:
        pts = map(float, line.split(','))
        # plot points
        x, y = m(pts[0],pts[1])
        p = m.plot(x,y,'o',color='purple',markersize=5)
    
    return p

# calculate lat lon from range (km) and bearing (degrees)
# code adapted from: http://stackoverflow.com/questions/7222382/get-lat-long-given-current-point-distance-and-bearing
def reckon(lat1d, lon1d, rngkm, brngd):
    from math import radians, asin, sin, cos, atan2, degrees   
    
    R = 6378.1 #Radius of the Earth
    brng = radians(brngd)
    
    lat1r = radians(lat1d) #Current lat point converted to radians
    lon1r = radians(lon1d) #Current long point converted to radians
    
    lat2r = asin(sin(lat1r)*cos(rngkm/R) + cos(lat1r)*sin(rngkm/R)*cos(brng))

    lon2r = lon1r + atan2(sin(brng)*sin(rngkm/R)*cos(lat1r), \
            cos(rngkm/R)-sin(lat1r)*sin(lat2r))

    lat2d = degrees(lat2r)
    lon2d = degrees(lon2r)
    
    return [lon2d, lat2d]

# check to see if a point is inside a polygon
def checkPointInPoly(shpfile, lon, lat):
    import shapefile
    from shapely.geometry import Point, Polygon

        
    #parese shapefile
    sf = shapefile.Reader(shpfile)
    shapes = sf.shapes()
    poly = Polygon(shapes[0].points)
    
    pt = Point(lon, lat)
    
    if pt.within(poly):
        value = True
    else:
        value = False
    return value
    
def get_line_parallels(pts, rngkm):
    from obspy.geodetics.base import gps2dist_azimuth
    from mapping_tools import reckon
    '''
    pts are an N x [lon, lat] matrix, i.e.:
                   [[lon1, lat1],
                    [lon2, lat2]]
    '''    
    # set outputs
    posazpts = []
    negazpts = []
    
    for j, pt in enumerate(pts):
        # if 1st point
        if j == 0:
            rngm, az, baz = gps2dist_azimuth(pts[j][1], pts[j][0], \
                                            pts[j+1][1], pts[j+1][0])
            
        # if last point
        elif j == len(pts)-1:
            rngm, az, baz = gps2dist_azimuth(pts[j-1][1], pts[j-1][0], \
                                            pts[j][1], pts[j][0])
                                           
        # use points either side (assumes evenly spaced)
        else:
            rngm, az, baz = gps2dist_azimuth(pts[j-1][1], pts[j-1][0], \
                                            pts[j+1][1], pts[j+1][0])
           
        # get azimuth for new points
        azpos = az + 90.
        azneg = az - 90.
        # get points
        posazpts.append(reckon(pts[j][1], pts[j][0], rngkm, azpos))
        negazpts.append(reckon(pts[j][1], pts[j][0], rngkm, azneg))
    
    '''    
    # for testing only
    x=[]
    y=[]
    xx=[]
    yy=[]
    xxx=[]
    yyy=[]
    for j, pt in enumerate(pts):
        x.append(pt[0])
        y.append(pt[1])
        xx.append(posazpts[j][0])
        yy.append(posazpts[j][1])
        xxx.append(negazpts[j][0])
        yyy.append(negazpts[j][1])
        
    plt.plot(x,y,'b-')
    plt.plot(xx,yy,'r-')
    plt.plot(xxx,yyy,'g-')
    plt.show()
    '''
    return posazpts, negazpts

# renames obspy tool to something I remember
# returns rngkm (km), az, baz (degrees)
def distance(lat1, lon1, lat2, lon2):
    #from obspy.core.util.geodetics import gps2DistAzimuth
    from obspy.geodetics.base import gps2dist_azimuth
    
    #print(lat1, lon1, lat2, lon2)
    
    #rngm, az, baz = gps2DistAzimuth(lat1, lon1, lat2, lon2)
    rngm, az, baz = gps2dist_azimuth(lat1, lon1, lat2, lon2)
    
    rngkm = rngm / 1000.
    
    return rngkm, az, baz

def km2deg(rngkm):
    from obspy.geodetics.base import kilometer2degrees
    
    return kilometer2degrees(rngkm, radius=6371.)

# function to get linear profile of lat/lons    
def get_distance_locs(lat1, lon1, lat2, lon2, npts):
    from mapping_tools import reckon, distance
    from numpy import array
    
    '''
    npts = 1 more than number of steps
    '''
    '''
    lon1 = 131.0
    lat1 = -11.5
    lon2 = 134.5
    lat2 = -17.0
    '''
    rngkm, az, baz = distance(lat1, lon1, lat2, lon2)
    
    # get inc distance
    inckm = rngkm / (npts-1)
    
    plats = []
    plons = []
    csvtxt = ''
    
    for i in range(0, npts):
    
        lon2d, lat2d = reckon(lat1, lon1, i*inckm, az)
        plats.append(lat2d)
        plons.append(lon2d)
        
        csvtxt += ','.join((str('%0.4f' % lon2d), str('%0.4f' % lat2d))) + '\n'
        
    f = open('profile.csv', 'wb')
    f.write(csvtxt)
    f.close()
    
    return array(plons), array(plats)


# generates a vector of distance, azimuth & back azimuth using "distance"
# returns rngkm (km), az, baz (degrees)
def dist_vect(lastlat, lastlon, latvect, lonvect):
   from numpy import array
   from mapping_tools import distance

   rng = []
   az  = []
   baz = []

   for l in range(0,len(latvect)):
       rngtmp, aztmp, baztmp = distance(lastlat, lastlon, latvect[l], lonvect[l])
       rng.append(rngtmp)
       az.append(aztmp)
       baz.append(baztmp)

   return array(rng), array(az), array(baz)

# read shapfile and get plotting coords in metres
def get_mapping_extent_from_shp(shpfile, lonoff):
    import shapefile
    from numpy import mean, ceil, floor
    
    sf = shapefile.Reader(shpfile)
    shapes = sf.shapes()

    bbox = [180, 90, -180, -90]
    padboxlat = 1
    padboxlon = 1
    for shape in shapes:
        sbbox = shape.bbox
        if sbbox[0] < bbox[0]:
            bbox[0] = sbbox[0]
        if sbbox[1] < bbox[1]:
            bbox[1] = sbbox[1]
        if sbbox[2] > bbox[2]:
            bbox[2] = sbbox[2]
        if sbbox[3] > bbox[3]:
            bbox[3] = sbbox[3]
    
    # pad bounding box and round to nearest degree
    bbox[0] = floor(bbox[0] - padboxlon)
    bbox[1] = floor(bbox[1] - padboxlat)
    bbox[2] = ceil(bbox[2] + padboxlon)
    bbox[3] = ceil(bbox[3] + padboxlat)
    
    # get central lon/lat
    lon_0 = mean([bbox[0], bbox[2]]) + float(lonoff)
    lat_0 = mean([bbox[1], bbox[3]])
    
    # get xy extents
    xkm, az, baz = distance(lat_0, bbox[0], lat_0, bbox[2])
    ykm, az, baz = distance(bbox[1], lon_0, bbox[3], lon_0)
    
    xm = 1000 * xkm
    ym = 1000 * ykm
    
    return lon_0, lat_0, xm, ym
    
# read shapfile and get plotting coords
def get_mapping_extent_in_metres(bbox, lonoff):
    '''
    bbox is an array of [lonmin, lonmax, latmin, latmax]
    '''
    from numpy import mean, ceil
        
    # get central lon/lat
    lon_0 = mean([bbox[0], bbox[1]]) + float(lonoff)
    lat_0 = mean([bbox[2], bbox[3]])
    
    # get xy extents
    xkm, az, baz = distance(lat_0, bbox[0], lat_0, bbox[1])
    ykm, az, baz = distance(bbox[2], lon_0, bbox[3], lon_0)
    
    xm = 1000 * xkm
    ym = 1000 * ykm
    
    return lon_0, lat_0, xm, ym
    
# get 
def get_map_extent(gmtbounds):
    
    '''
    assumes GMT format lon/lat params
        e.g. gmtbounds = [llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat]
    '''
    
    from numpy import mean
    llcrnrlon = gmtbounds[0]
    urcrnrlon = gmtbounds[1]
    llcrnrlat = gmtbounds[2]
    urcrnrlat = gmtbounds[3]
    lon_0 = mean([gmtbounds[0], gmtbounds[1]])
    lat_0 = mean([gmtbounds[2], gmtbounds[3]])
    
    return llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat, lat_0, lon_0

# sets map extent for Lambert Conic Conformal projection    
def set_lcc_basemap(llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat, resolution):
    from mpl_toolkits.basemap import Basemap
    from numpy import mean, percentile
    
    # get map bounds
    lon_0 = mean([llcrnrlon, urcrnrlon])
    lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
    lat_2 = percentile([llcrnrlat, urcrnrlat], 75)
    
    # set map
    return Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
                urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,
                projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,
                resolution=resolution,area_thresh=100.)
                
# gets land polygons from matplotlib basemap
def get_map_polygons(m):
    polygons = []
    for polygon in m.landpolygons:
        coords = polygon.get_coords()
        tuplelist = []
        for c in coords:
            tuplelist.append((c[0], c[1]))
        polygons.append(tuplelist)
        
    return polygons
    
def fill_lake_polygons(m, plt, facecolor):
    polygons = []
    for polygon in m.lakepolygons:
        poly = polygon.get_coords()
        plt.fill(poly[:,0], poly[:,1], facecolor) # could do outside function?
        polygons.append(poly)
        
    return polygons

# for masking one poly    
def mask_outside_polygon(poly_verts, facecolor, plt):
    """
    Plots a mask on the specified axis ("ax", defaults to plt.gca()) such that
    all areas outside of the polygon specified by "poly_verts" are masked.  

    "poly_verts" must be a list of tuples of the verticies in the polygon in
    counter-clockwise order.

    Returns the matplotlib.patches.PathPatch instance plotted on the figure.
    """
    import matplotlib.patches as mpatches
    import matplotlib.path as mpath
    from numpy import array, vstack

    #if ax is None:
    ax = plt.gca()

    # Get current plot limits
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    # Verticies of the plot boundaries in clockwise order
    bound_verts = [(xlim[0], ylim[0]), (xlim[0], ylim[1]), 
                   (xlim[1], ylim[1]), (xlim[1], ylim[0]), 
                   (xlim[0], ylim[0])]

    # A series of codes (1 and 2) to tell matplotlib whether to draw a line or 
    # move the "pen" (So that there's no connecting line)
    bound_codes = [mpath.Path.MOVETO] + (len(bound_verts) - 1) * [mpath.Path.LINETO]
    print(bound_codes)
    poly_codes = [mpath.Path.MOVETO] + (len(poly_verts) - 1) * [mpath.Path.LINETO]
    print(poly_codes)
    
    # Plot the masking patch
    path = mpath.Path(bound_verts + poly_verts, bound_codes + poly_codes)
    #path = mpath.Path(vstack((bound_verts, poly_verts)), vstack((bound_codes, poly_codes)))
    patch = mpatches.PathPatch(path, facecolor='white', edgecolor='none')
    patch = ax.add_patch(patch)

    # Reset the plot limits to their original extents
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    return patch

# for masking multiple non-conected polys
def mask_outside_polygons(polys, facecolor, plt):
    """
    Plots a mask on the specified axis ("ax", defaults to plt.gca()) such that
    all areas outside of the polygon specified by "poly_verts" are masked.  

    "poly_verts" must be a list of tuples of the verticies in the polygon in
    counter-clockwise order.

    Returns the matplotlib.patches.PathPatch instance plotted on the figure.
    """
    import matplotlib.patches as mpatches
    import matplotlib.path as mpath

    #if ax is None:
    ax = plt.gca()

    # Get current plot limits
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    # Verticies of the plot boundaries in clockwise order
    bound_verts = [(xlim[0], ylim[0]), (xlim[0], ylim[1]), 
                   (xlim[1], ylim[1]), (xlim[1], ylim[0]), 
                   (xlim[0], ylim[0])]

    # A series of codes (1 and 2) to tell matplotlib whether to draw a line or 
    # move the "pen" (So that there's no connecting line)
    bound_codes = [mpath.Path.MOVETO] + (len(bound_verts) - 1) * [mpath.Path.LINETO]
    verts = bound_verts
        
    # loop thru multiple polys
    codes = bound_codes
    for poly_verts in polys:
        verts += poly_verts[::-1]
        codes += [mpath.Path.MOVETO] + (len(poly_verts) - 1) * [mpath.Path.LINETO]

    # Plot the masking patch
    path = mpath.Path(verts, codes)
    #path = mpath.Path(vstack((bound_verts, poly_verts)), vstack((bound_codes, poly_codes)))
    patch = mpatches.PathPatch(path, facecolor=facecolor, edgecolor='none')
    patch = ax.add_patch(patch)

    # Reset the plot limits to their original extents
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    return patch

def get_WGS84_area(geom):
    '''
    geom = shapely polygon
    '''
    import pyproj    
    import shapely.ops as ops
    #from shapely.geometry.polygon import Polygon
    from functools import partial
    
    try:
        geom_area = ops.transform(
            partial(
                pyproj.transform,
                pyproj.Proj(init='EPSG:4326'),
                pyproj.Proj(
                    proj='aea',
                    lat_1=geom.bounds[1],
                    lat_2=geom.bounds[3])),
            geom)
    # different keys for old lat/lon
    except:
        geom_area = ops.transform(
            partial(
                pyproj.transform,
                pyproj.Proj(init='EPSG:4326'),
                pyproj.Proj(
                    proj='aea',
                    lat1=geom.bounds[1],
                    lat2=geom.bounds[3])),
            geom)
    
    # print(the area in km^2)
    #print(geom_area.area / 1000000.
    
    return geom_area.area / 1000000.
    
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
        print("file ",fileName, "not found")
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


# !!!!!!!!!!code to convert projections!!!!!!!!!!
# http://all-geo.org/volcan01010/2012/11/change-coordinates-with-pyproj/
# http://stackoverflow.com/questions/26452972/coordinates-conversion-with-pyproj
"""
nc = NetCDFFile('//home//tallen//DATA//Population//2011_census_pop_canada_lcc83.grd')

data = nc.variables['z'][:] # persons per square km
eastings = nc.variables['x'][:]
northings = nc.variables['y'][:]

# Define a projection with Proj4 notation
#nad1983=pyproj.Proj("+proj=lcc +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=1700000 +y_0=300000 +no_defs +a=6378137 +rf=298.257222101 +to_meter=1")
nad1983=pyproj.Proj("+proj=lcc +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +no_defs +a=6378137 +rf=298.257222101 +to_meter=1") # as per Will's geodatabase
 
# Define some common projections using EPSG codes
wgs84=pyproj.Proj("+init=EPSG:4326") # LatLon with WGS84 datum used by GPS units and Google Earth
#NAD83=pyproj.Proj("+init=EPSG:3347") # NAD83 / Statistics Canada Lambert

# get lat/lon for each point:
print('Converting to lon/lat...')
xx, yy = meshgrid(eastings, northings)
lons, lats = pyproj.transform(nad1983, wgs84, xx.flatten(), yy.flatten())

print('Looping thru pop data...')
# write population file for grdtrack
pop = []
plo = []
pla = []
for lo, la, d in zip(lons, lats, data.flatten()):
    if d > 0.0:
        plo.append(lo)
        pla.append(la)
        pop.append(d)

plo = array(plo).reshape((len(plo),1))
pla = array(pla).reshape((len(pla),1))
pop = array(pop).reshape((len(pop),1))

darray = hstack((plo, pla, pop))
savetxt('can_pop_dat.txt', darray, delimiter='\t', fmt='%0.3f')
"""

def annotate_cities(numCities, plt, m, markerfacecolor='k', markeredgecolor='k', \
                    marker='o', markersize=6, markeredgewidth=0.5, fs=14):
    from numpy import argsort
    from os import getcwd
    import matplotlib.patheffects as PathEffects
    path_effects=[PathEffects.withStroke(linewidth=3, foreground="w")]
    
    # set grid size (in degrees)
    lonrng = m.urcrnrlon - m.llcrnrlon
    print('lonrng', lonrng)
    if lonrng <= 0.4:
        pltbuffer = 0.04
        txtoff = 0.0035 
    elif lonrng <= 0.61:
        pltbuffer = 0.04
        txtoff = 0.005 
    elif lonrng <= 1:
        pltbuffer = 0.04
        txtoff = 0.006
    elif lonrng <= 2:
        pltbuffer = 0.06
        txtoff = 0.012 
    elif lonrng <= 2.5:
        pltbuffer = 0.075
        txtoff = 0.015  
    elif lonrng <= 3:
        pltbuffer = 0.1
        txtoff = 0.017
    elif lonrng > 3 and lonrng < 6:
        pltbuffer = 0.175
        txtoff = 0.03
    elif lonrng >= 6:
        pltbuffer = 0.7
        txtoff = 0.1
    #print(pltbuffer
    
    # parse AU cities
    if getcwd().startswith('/nas'):
        cityFile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/REQIT/GA-shakemap/mapping/cities1000_au_ascii.txt'
    else:
        cityFile = '//Users//trev//Documents//DATA//Places//cities1000_au_ascii.txt'
    
    lines = open(cityFile).readlines()
    
    # label cities
    clatList = []
    clonList = []
    
    # for ordering in plot priority
    pltlo = []
    pltla = []
    pltpop = []
    pltloc = []
    
    for line in lines:
        dat = line.strip().split('\t')
        clat = float(dat[4])
        clon = float(dat[5])
        #loc = dat[1]
        
        if clat > m.llcrnrlat and clat < m.urcrnrlat-txtoff \
           and clon > m.llcrnrlon and clon < m.urcrnrlon:
            
            # add city to list
            pltlo.append(clon)
            pltla.append(clat)
            pltpop.append(float(dat[14]))
            pltloc.append(dat[1])
    
    # order locs in trems of poplation
    sortidx = argsort(pltpop)[::-1] # and reverse
    
    # now loop through ordered cities - use max 15
    i = 0
    for si in sortidx:
        # plt max num locs
        if i < numCities:
            pltCity = True
            print(pltloc[si])
            
            for clol, clal in zip(clonList, clatList):
                if abs(pltla[si] - clal) < pltbuffer and abs(pltlo[si] - clol) < (3.*pltbuffer):
                    #print(pltloc[si]
                    pltCity = False
                    
            # build list of locs
            clatList.append(pltla[si])
            clonList.append(pltlo[si])
            
            if pltCity == True:
                x, y = m(pltlo[si], pltla[si])
                plt.plot(x, y, marker, markerfacecolor=markerfacecolor, \
                         markeredgecolor=markeredgecolor, markeredgewidth=markeredgewidth, \
                         markersize=markersize, zorder=11000)
        
                x, y = m(pltlo[si]-txtoff, pltla[si]+txtoff) # changed for camden fig
                plt.text(x, y, pltloc[si], size=fs, ha='right', weight='normal', path_effects=path_effects, zorder=11000)
                
                i += 1
                
    return clatList, clonList
                                      
# code based on: http://qingkaikong.blogspot.com/2016/06/nice-python-basemap-background.html
def make_street_map(clat, clon, service='ESRI_Imagery_World_2D', ll_buffer = 0.1, \
             xpixels = 1500, plt_inset = True, inset_state = 'nsw', inset_loc = 3, \
             inset_multiplier=0.03, plt_marker = True, marker='*', ms = 14, mew = 0.5, mfc = 'r', mec='k'):
    
    from mpl_toolkits.basemap import Basemap
    import matplotlib.pyplot as plt
    from numpy import array
    
    '''
    Map Services:
        ESRI_Imagery_World_2D
        ESRI_StreetMap_World_2D
        NatGeo_World_Map
        NGS_Topo_US_2D
        Ocean_Basemap
        USA_Topo_Maps
        World_Imagery
        World_Physical_Map
        World_Shaded_Relief
        World_Street_Map
        World_Terrain_Base
        World_Topo_Map
        
    Example:
        from mapping_tools import make_street_map
        plt, m = make_street_map(-34.070702, 150.480499, service='NatGeo_World_Map', marker='t')
    '''
    # set map bounds
    urcrnrlat = clat + ll_buffer
    urcrnrlon = clon + ll_buffer
    llcrnrlat = clat - ll_buffer
    llcrnrlon = clon - ll_buffer

    # EPSG Projection 3577 - GDA94 / Australian Albers
    # EPSG Projection 3112 - GDA94 / Geoscience Australia Lambert
    # EPSG:3308 - GDA94 / NSW Lambert 
    epsg = 3308
    
    # note, you need change the epsg for different region, 
    #US is 4269, and you can google the region you want
    plt.figure(1, figsize = (10, 10))
    ax = plt.subplot(111)
    
    m = Basemap(projection='mill',llcrnrlon=llcrnrlon ,llcrnrlat=llcrnrlat,
        urcrnrlon=urcrnrlon ,urcrnrlat=urcrnrlat, resolution = 'h', epsg = epsg)
    
    # xpixels controls the pixels in x direction, and if you leave ypixels
    # None, it will choose ypixels based on the aspect ratio
    m.arcgisimage(service=service, xpixels = xpixels, verbose= True, dpi=600)
    
    # plot marker
    if plt_marker == True:
        x, y = m(clon, clat)
        plt.plot(x, y, marker=marker, ms=ms, mfc=mfc, mew=mew, mec=mec)
    
    # make inset
    if plt_inset == True:
        if inset_state == 'nsw':
            iurcrnrlat = -27.7
            iurcrnrlon = 154.2
            illcrnrlat = -38.
            illcrnrlon = 140.4
        
        from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
        print(inset_loc)
        axins = zoomed_inset_axes(ax, ll_buffer*inset_multiplier, loc=inset_loc)
        
        m2 = Basemap(projection='merc',llcrnrlon=illcrnrlon ,llcrnrlat=illcrnrlat,
                     urcrnrlon=iurcrnrlon ,urcrnrlat=iurcrnrlat, resolution = 'l')
                    
        m2.drawmapboundary(fill_color='0.8')
        m2.fillcontinents(color='w', lake_color='0.8')
        m2.drawcoastlines()
        m2.drawcountries()
        m2.drawstates()
        
        # fill main area
        #xv = array([llcrnrlon, llcrnrlon, urcrnrlon, urcrnrlon, llcrnrlon])
        #yv = array([llcrnrlat, urcrnrlat, urcrnrlat, llcrnrlat, llcrnrlat])
        #x, y = m2(xv, yv)
        #plt.fill(x, y, 'r-', lw=1.5)
        x, y = m2(clon, clat)
        plt.plot(x, y, 'rs')

    return plt, m, ax