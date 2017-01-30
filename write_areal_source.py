# -*- coding: utf-8 -*-
"""
Created on Wed Apr 03 13:34:10 2013

@author: tallen

List of functions:
    beta2bval: convert beta to b-value
    bval2beta: convert b-value to beta
    read_szonefile*: function to read *.zon files
    read_betafile*: function to read output from betaplMX (*.beta files)
    make_parts: makes lon/lat pairs for outputting shapefiles
    write_src_shape: sets up area and fault source data for writing to shapefile
    fill_shape_values: fills shapefile values extracted from gscfrisk inputs
    write_src_shapewrite_src_shape: writes shapefile from gscfrisk inputs
    write_cmpfile: function to write *.cmp file for input to szonegmt
    write_szonegmt: function to write szonegmt.inp for input to szonegmt
    

    * will be removed
    
"""

'''
# convert beta to b-value
'''
def beta2bval(beta):
    from numpy import log10, exp
    return log10(exp(beta))
    
'''    
# convert b-value to beta
'''
def bval2beta(bval):
    from numpy import log
    return log(10**bval)   

'''
# function to read *.zon files
'''
def read_szonefile(zonefile):
    data = open(zonefile).readlines()
    zoneheader = data[0].strip('\n')
    
    # now get coordinates
    lat = []
    lon = []
    for i in range(1,len(data)-1): # by default, close polygon
        dat = data[i].split()
        lat.append(float(dat[0]))
        lon.append(float(dat[1]))
        
    return zoneheader, lat, lon

'''        
# function to read output from betaplMX (*.beta files)
'''
def read_betafile(betafile):

    data = open(betafile).readlines()
    getnextline = 0
    
    for line in data:   
        # get lower curve parameters
        if line.find('LOWER CURVE:') >= 0:
            dat = line.split()
            lcparam = dat[2:]
            
        # get best curve parameters
        if line.find('BEST ESTIMATE:') >= 0:
            dat = line.split()
            bcparam = dat[2:]
            
        # get upper curve parameters
        if line.find('UPPER CURVE:') >= 0:
            dat = line.split()
            ucparam = dat[2:]
            
        # get gscfrisk input line
        if getnextline == 1:
            gscfrisk = line.strip('\n')
            getnextline = 0
            
        if line.find('GSCFRISK MODEL INPUT') >= 0:
            getnextline = 1
            
    return lcparam, bcparam, ucparam, gscfrisk
        
def make_parts(mod):
    mod1 = []
    for i in range(0, len(mod)):
        mod1.append([mod[i,0], mod[i,1]])
    
    mod2 = []
    mod2.append(mod1)
    return mod2

# this function sets up area and fault source data for writing to shapefile
def write_src_shape(mod, modelfile, altzones):
    import shapefile
    
    # set file name
    tmpname = modelfile.strip().split('.')
    s = ''
    tmpfile = s.join(tmpname[0:-1])
    outfile = tmpfile
    outshp = outfile+'_area'+str(altzones)
    
    # Create a polygon shapefile for areal source zones
    w = shapefile.Writer(shapefile.POLYGON)
    
    # check to see if areas exist
    check_model = True
    i = 0
    while check_model == True:
        if mod[i]['src_type'] == 'area':
            # send data for writing
            fill_shape_values(w, mod, outshp)
            check_model = False
        i += 1
        
        if i == len(mod):
            check_model = False

    
    # Create a line shapefile for faults
    # set file name
    tmpfile = s.join(tmpname[0:-1])
    outfile = tmpfile
    outshp = outfile+'_fault'+str(altzones)
    
    w = shapefile.Writer(shapefile.POLYLINE)
    # check to see if faults exist
    check_model = True
    i = 0
    while check_model == True:
        if mod[i]['src_type'] == 'fault':
            # send data for witing
            fill_shape_values(w, mod, outshp)
            check_model = False
        i += 1
        
        if i == len(mod):
            check_model = False

# this function populates shapefile values and writes to file
def fill_shape_values(w, mod, outshp):
    w.field('CODE','C','100')
    w.field('NAME','C','100')
    w.field('DEP_BEST','F', 13, 6)
    w.field('DEP_LOWER','F', 13, 6)
    w.field('DEP_UPPER','F', 13, 6)
    w.field('MIN_MAG','F', 13, 6)
    w.field('MMAX_BEST','F', 13, 6)
    w.field('MMAX_LOWER','F', 13, 6)
    w.field('MMAX_UPPER','F', 13, 6)
    w.field('A0_BEST','F', 13, 6)
    w.field('A0_LOWER','F', 13, 6)
    w.field('A0_UPPER','F', 13, 6)
    w.field('BETA_BEST','F', 13, 6)
    w.field('BETA_LOWER','F', 13, 6)
    w.field('BETA_UPPER','F', 13, 6)    
    
    for i in range(0, len(mod)):    
        # now loop through polygons
        if mod[i]['src_type'] == 'area':
            poly_parts = make_parts(mod[i]['src_shape'])
            w.poly(parts=poly_parts)
            
            w.record(mod[i]['src_name'], mod[i]['src_reg'], \
                mod[i]['src_dep'][0], mod[i]['src_dep'][2], mod[i]['src_dep'][1], \
                mod[i]['min_mag'], \
                mod[i]['max_mag'][0], mod[i]['max_mag'][1], mod[i]['max_mag'][2], \
                mod[i]['src_A0'][0], mod[i]['src_A0'][2], mod[i]['src_A0'][1], \
                mod[i]['src_bval'][0], mod[i]['src_bval'][2], mod[i]['src_bval'][1])
            
    # now save area shapefile
    print outshp
    w.save(outshp)

# function to write *.cmp file for input to szonegmt    
def write_cmpfile(moddir, params):
    from os import path, mkdir
    
    outfile = path.join(moddir,params['zcode'],params['zcode']+'.cmp')
    
    # now check to make sure output path exists
    try:
        f = open(outfile,'wb')
    except:
        if path.exists(path.join(moddir,params['zcode'])) == False:
            mkdir(path.join(moddir,params['zcode']))
        f = open(outfile,'wb')
    
    # now start preparing text
    cmptxt = 'T\n' + params['zcode'] + ' - ' + params['zname'] + '\n' \
             + '  '.join(('',str("%0.1f" % params['mmin']),str("%0.1f" % params['mxbest']), \
                          str("%0.1f" % params['msep']),str(len(params['mcmp'])))) + '\n'
    for i in range(0,len(params['mcmp'])):
        cmptxt = cmptxt + '  '.join(('',str("%0.1f" % params['mcmp'][i]), \
                               str(params['ycmp'][i]))) + '\n'
    
    cmptxt = cmptxt + '   ' + str(params['ymax']) + "\n'(A2,I4,10X,F3.1)'" + '\n'
             
        
    # write completeness file
    f.write(cmptxt)
    f.close()

# function to write betapl.inp file for input to betaplMRx    
def write_betainp(moddir, params):
    from os import path, mkdir
    from numpy import isnan
    
    outfile = path.join(moddir,params['zcode'],'betapl.inp')
    
    # now check to make sure output path exists
    try:
        f = open(outfile,'wb')
    except:
        if path.exists(path.join(moddir,params['zcode'])) == False:
            mkdir(path.join(moddir,params['zcode']))
        f = open(outfile,'wb')
        
    # now start preparing text
    betatxt = params['zcode'] + '.incmp' + '\n' \
              + params['zcode'] + '.beta' + '\n' \
              + params['zcode'] + '.beta.ps' + '\n'
              
    # check if fixed beta
    if isnan(params['fixedbetaval']):
        fixtxt = 'N\n'
    else:
        fixtxt = 'Y\n' + str(params['fixedbetaval']) + ', ' \
                       + str(params['fixedbetasig']) + '\n'
    
    betatxt = betatxt + fixtxt
    betatxt = betatxt + 'Y\nY\nY\n'
    betatxt = betatxt + str("%0.1f" % params['mxlower']) + ', ' \
                      + str("%0.1f" % params['mxupper']) + '\n'
    betatxt = betatxt + '1\n'
              
    # write betapl.inp file
    f.write(betatxt)
    f.close()
    
# function to write szonegmt.inp for input to szonegmt
def write_szonegmt(moddir, params):
    from os import path, mkdir
    
    outfile = path.join(moddir,params['zcode'],'szonegmt.inp')
    
    # now check to make sure output path exists
    try:
        f = open(outfile,'wb')
    except:
        if path.exists(path.join(moddir,params['zcode'])) == False:
            mkdir(path.join(moddir,params['zcode']))
        f = open(outfile,'wb')
    
    sztxt = params['zcode'] + '\n' \
          + params['zcode'] + '.zon\n' \
          + params['zcode'] + '.cmp\n' \
          + params['catalogue'] + '\n'
    
    # write betapl.inp file
    f.write(sztxt)
    f.close()
    
def szone2shp(zonfile, outshp, srctype):
    from parse_areal_source import read_szonefile
    import shapefile
    '''
    srctype = fault or area
    '''
    zoneheader, lat, lon = read_szonefile(zonfile, srctype)
    
    # make mod
    mod1 = []
    mod2 = []
    for i in range(0, len(lat)):
        mod1.append([lon[i], lat[i]])
    
    # now prep shapefile
    if srctype == 'area':
        w = shapefile.Writer(shapefile.POLYGON)
        # check to see of polygon closes
        if lon[0] != lon[-1] or lat[0] != lat[-1]:
            mod1.append([lon[0], lat[0]])
    else:
        w = shapefile.Writer(shapefile.POLYLINE)
   
    mod2.append(mod1)
    w.field('NAME','C','50')
    w.record(zoneheader)
    w.line(parts=mod2)
    w.save(outshp)
    
    
    































