# -*- coding: utf-8 -*-
"""
List of functions:
    read_gscfrisk_areal: reads the GSC FRISK input data files
    beta2bval: convert beta to b-value
    bval2beta: convert b-value to beta
    read_szonefile: function to read *.zon files
    read_betafile: function to read output from betaplMRX (*.beta files)
    read_betafile: function to read output from betaplMRX (*.beta files)
    read_cmpfile: function to read completeness files input to szonegmt
    read_betaplfile: function to read input to betaplMRX
    read_zonefile: function to read input *.zon
    read_hazparam_lookup: reads csv lookup of recurrence model parameters for all zones
    
"""

def make_parts(mod):
    mod1 = []
    for i in range(0, len(mod)):
        mod1.append([mod[i,0], mod[i,1]])
    
    mod2 = []
    mod2.append(mod1)
    return mod2

# this function populates area shapefile values and writes to file
def fill_shape_values(w, mod, outshp):
    w.field('SRC_NAME','C','100')
    w.field('CODE','C','10')
    w.field('SRC_REGION','C','100')
    w.field('SRC_TYPE','C','10')
    w.field('SRC_WEIGHT','F', 13, 6)
    w.field('DEP_BEST','F', 13, 6)
    w.field('DEP_LOWER','F', 13, 6)
    w.field('DEP_UPPER','F', 13, 6)
    w.field('MIN_MAG','F', 13, 6)
    w.field('MMAX_BEST','F', 13, 6)
    w.field('MMAX_LOWER','F', 13, 6)
    w.field('MMAX_UPPER','F', 13, 6)
    w.field('N0_BEST','F', 13, 6)
    w.field('N0_LOWER','F', 13, 6)
    w.field('N0_UPPER','F', 13, 6)
    w.field('BETA_BEST','F', 13, 6)
    w.field('BETA_LOWER','F', 13, 6)
    w.field('BETA_UPPER','F', 13, 6)
    w.field('SRC_GMPE','C','100')
    
    for i in range(0, len(mod)):    
        # now loop through polygons
        if mod[i]['src_type'] == 'area':
            poly_parts = make_parts(mod[i]['src_shape'])
            w.poly(parts=poly_parts)
            code = mod[i]['src_name'].split()[0]
            
            w.record(mod[i]['src_name'], code, mod[i]['src_reg'], mod[i]['src_type'], \
                mod[i]['src_weight'], mod[i]['src_dep'][0], max(mod[i]['src_dep']), min(mod[i]['src_dep']), \
                mod[i]['min_mag'], \
                mod[i]['max_mag'][0], mod[i]['max_mag'][1], mod[i]['max_mag'][2], \
                mod[i]['src_N0'][0], mod[i]['src_N0'][1], mod[i]['src_N0'][2], \
                mod[i]['src_beta'][0], mod[i]['src_beta'][1], mod[i]['src_beta'][2], \
                mod[i]['gmpe'])
            
    # now save area shapefile
    print outshp
    w.save(outshp)

# this function populates area shapefile values and writes to file
def fill_shape_values_fault(w, mod, outshp):
    import shapefile 
    
    w.field('SRC_NAME','C','100')
    w.field('CODE','C','10')
    w.field('SRC_REGION','C','100')
    w.field('SRC_TYPE','C','10')
    w.field('SRC_WEIGHT','F', 13, 6)
    w.field('DIP_LOWER','F', 13, 6)
    w.field('DIP_UPPER','F', 13, 6)    
    w.field('DEP_UPPER','F', 13, 6)
    w.field('DEP_MIDDLE','F', 13, 6)
    w.field('DEP_LOWER','F', 13, 6)    
    w.field('MIN_MAG','F', 13, 6)
    w.field('MMAX_BEST','F', 13, 6)
    w.field('MMAX_LOWER','F', 13, 6)
    w.field('MMAX_UPPER','F', 13, 6)
    w.field('N0_BEST','F', 13, 6)
    w.field('N0_LOWER','F', 13, 6)
    w.field('N0_UPPER','F', 13, 6)
    w.field('BETA_BEST','F', 13, 6)
    w.field('BETA_LOWER','F', 13, 6)
    w.field('BETA_UPPER','F', 13, 6)
    w.field('SRC_GMPE','C','100')
    
    for i in range(0, len(mod)):    
        # now loop through polygons
        if mod[i]['src_type'] == 'fault':
            poly_parts = make_parts(mod[i]['src_shape'])
            w.line(parts=poly_parts, shapeType=shapefile.POLYLINE)
            code = mod[i]['src_name'].split()[0]
            
            w.record(mod[i]['src_name'], code, mod[i]['src_reg'],mod[i]['src_type'], \
                mod[i]['src_weight'], mod[i]['fault_dip'][0], mod[i]['fault_dip'][1], \
                mod[i]['src_dep'][0], mod[i]['src_dep'][1], mod[i]['src_dep'][2], \
                mod[i]['min_mag'], \
                mod[i]['max_mag'][0], mod[i]['max_mag'][1], mod[i]['max_mag'][2], \
                mod[i]['src_N0'][0], mod[i]['src_N0'][1], mod[i]['src_N0'][2], \
                mod[i]['src_beta'][0], mod[i]['src_beta'][1], mod[i]['src_beta'][2], \
                mod[i]['gmpe'])
            
    # now save area shapefile
    print outshp
    w.save(outshp)
    
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
            fill_shape_values_fault(w, mod, outshp)
            check_model = False
        i += 1
        
        if i == len(mod):
            check_model = False


# this function reads the GSC FRISK input data files and writes out shapefiles
def read_gscfrisk_areal(modelfile):
    from numpy import nan, array, zeros, vstack
    import pickle
    
    # for testing 
    #modelfile = 'U:\\2015_Hazard_Model\\Source_Zones\\2013_TE\\SECan_T3E.model'
    #modelfile = 'U:\\2015_Hazard_Model\\Source_Zones\\W2012_Trial3_area1_v3\\SWCan_T3E.model'
    # set variables
    skipshape = True    
    
    data = open(modelfile).readlines()
    
    # start reading source file
    i = 0
    altzones = 1
    
    while i < len(data):
        # find where new model is generated
        ind = data[i].find('1  1  1  1')
        if ind >= 0:
            # start writing new model
            model = []
            
            skipshape = False
            # get number of zones in model
            num_zones = data[i-4].strip().split()
            num_zones = int(num_zones[0])
    
            # now, loop thru zones in model
            for z in range(0, num_zones):
                # get number of sub-zones (i.e. for weightings)
                subzones = int(data[i+2].strip())
                
                for sz in range(0, subzones):
                    if subzones > 1:
                        src_name = data[i+3].strip()
                        if sz == 0:
                            src_reg = data[i+1].strip()
                    else:
                        src_name = data[i+1].strip()
                        src_reg = data[i+3].strip()
                                        
                    # set new dictionary
                    mod = {}
                    mod['src_name'] = src_name
                    mod['src_code'] = mod['src_name'].split()[0]
                    mod['src_reg'] = src_reg
                    mod['src_weight'] = float(data[i+4].strip())
                    mod['src_type'] = data[i+5].strip()
                    
                    # test if completeness info included              
                    
                    # check if using old model 4
                    ind = data[i+6].find('********')
                    if ind >= 0:
#                        mod['src_comp_mag'] = nan
#                        mod['src_comp_yr'] = nan
                        # iterate i
                        i += 9
                    elif len(data[i+6].strip().split()) == 1 and mod['src_type'] == 'area':
                        # get completeness
#                        mod['src_comp_mag'] = array(data[i+7].strip().split(), dtype='f')
#                        mod['src_comp_yr'] = array(data[i+8].strip().split(), dtype=int)
                        # iterate i
                        i += 9
                    elif len(data[i+6].strip().split()) > 5 and mod['src_type'] == 'area':
                        # get completeness
#                        mod['src_comp_mag'] = array(data[i+7].strip().split(), dtype='f')
#                        mod['src_comp_yr'] = array(data[i+8].strip().split(), dtype=int)
                        # iterate i
                        i += 9
                    elif mod['src_type'] == 'area':
#                        mod['src_comp_mag'] = nan
#                        mod['src_comp_yr'] = nan
                        # iterate i
                        i += 6
                    elif mod['src_type'] == 'fault':                        
#                        mod['src_comp_mag'] = nan
#                        mod['src_comp_yr'] = nan
                        #mod['src_name'] = mod['src_reg'].split()[0]
    
                        # check if Gen 4 or 5
                        if len(data[i+7].strip().split()) == 3:
                            i += 7
                        else:
                            i += 9
                    
                    # get source depth
                    if mod['src_type'] == 'area':
                        depstr = data[i].strip().split()
                        if len(depstr) == 3:
                            mod['src_dep'] = array(depstr[0:3], dtype='f') # not correct for fault src
                        elif len(depstr) == 2:
                            tmpdep = zeros((3,1))
                            tmpdep[0] = float(depstr[0])
                            mod['src_dep'] = tmpdep
                        mod['fault_dip'] = nan
                        
                    elif mod['src_type'] == 'fault':
                        depstr = data[i-3].strip().split()
                        mod['src_dep'] = array(depstr[2:], dtype='f')
                        mod['fault_dip'] = array(depstr[0:2], dtype='f')
                        
                    
                    # get source polygon
                    numlinesstr = data[i+1].strip().split()
                    numlines = int(numlinesstr[0])
                    src_shape = zeros((numlines, 2))
                    
                    i += 2
                    for j in range(0, numlines):
                        xystr = data[i].strip().split()
                        src_shape[j,0] = float(xystr[0])
                        src_shape[j,1] = float(xystr[1])
                        i += 1
                        
                    # close polygon
                    if mod['src_type'] == 'area':
                        if src_shape[0,0] != src_shape[-1,0] \
                        or src_shape[0,1] != src_shape[-1,1]:
                            src_shape = vstack((src_shape, src_shape[0,:]))
                    
                    # add shape to dict
                    mod['src_shape'] = src_shape
                    
                    # get Mmin
                    magstr = data[i].strip().split()
                    mod['min_mag'] = round(float(magstr[0]), 2)
                    #mod['min_mag'] = '%.2f' % float(magstr[0])
                    
                    # get Mmax
                    magstr = data[i].strip().split()
                    #mod['max_mag'] = round(array(magstr[1:4], dtype='f'), decimals=2)
                    mod['max_mag'] = [round(elem,2) for elem in array(magstr[1:4], dtype='f')]
                    
                    # get recurrence values
                    recstr = data[i+2].strip().replace('!',' ').split()
                    #mod['src_N0'] = array(recstr[0:6:2], dtype='f')
                    mod['src_N0'] = [round(elem,4) for elem in array(recstr[0:6:2], dtype='f')]
                    #mod['src_beta'] = array(recstr[1:7:2], dtype='f')
                    mod['src_beta'] = [round(elem,4) for elem in array(recstr[1:7:2], dtype='f')]
                    
                    # get GMPE
                    mod['gmpe'] = data[i+3].strip().split(',')[0]
                    
                    # check to see if generation 4 or 5
                    if i < len(data)-3:
                        ind = data[i+3].find('.txt')
                        
                        if ind >= 0:
                            gen5 = True
                        else:
                            gen5 = False
                        
                        if subzones > 1 and sz < subzones-1:
                            i += 1
                        else:
                            i += 3
                            
                        if gen5 == False:
                            i -= 1
                
                    # append zone to model
                    model.append(mod)
                    
            # write model to pickle
            pfile = ''.join((modelfile.split('.')[0], str(altzones), '.pkl'))
            pklfile = open(pfile, 'wb')
            pickle.dump(model, pklfile, -1)
            pklfile.close()
        i += 1
    
        # now write model to shapefile    
        if skipshape == False:
            write_src_shape(model, modelfile, altzones)
            altzones += 1
            skipshape = True
        
    return model
    
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
def read_szonefile(zonefile, srctype):
    data = open(zonefile).readlines()
    zoneheader = data[0].strip('\n')
    
    # now get coordinates
    lat = []
    lon = []
    
    if srctype == 'fault':
        tindex = 0
    elif srctype == 'area':
        tindex = 1
        
    for i in range(1,len(data)-tindex): # by default, close polygon
        dat = data[i].split()
        lat.append(float(dat[0]))
        lon.append(float(dat[1]))
        
    return zoneheader, lat, lon

'''        
# function to read output from betaplMX (*.beta files)
'''
def read_betafile(betafile):

    data = open(betafile,'rb').readlines()
    getnextline = 0
    
    # set values to null
    lcparam = 'null'
    bcparam = 'null'
    ucparam = 'null'
    gscfrisk = 'null'
    
    for line in data:
        # get mag range used
        if line.find('LOW AND HIGH MAGNITUDES USED:') >= 0:
            dat = line.split()
            mrng = dat[-2:]
            
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
            
    return lcparam, bcparam, ucparam, gscfrisk, mrng
    
'''        
# function to read completeness files (*.cmp) input to szonegmt
'''
def read_cmpfile(cmpfile):    
    cmpdict = {}
    data = open(cmpfile,'rb').readlines()
    
    # get zone name
    sep = ' - '
    dat = data[1].split('-')
    cmpdict['zcode'] = dat[0].strip()
    # rejoin zname if necessary
    cmpdict['zname'] = sep.join(dat[1:]).strip()
    
    # get mag info
    dat = data[2].split()
    cmpdict['mmin'] = float(dat[0])
    cmpdict['mxbest'] = float(dat[1])
    cmpdict['msep'] = float(dat[2])
    ncmp = int(dat[3])
    
    # now get completeness periods
    mcmp = []
    ycmp = []
    for i in range(3,3+ncmp):
        dat = data[i].split()
        mcmp.append(float(dat[0]))
        ycmp.append(int(dat[1]))
    
    cmpdict['mcmp'] = mcmp
    cmpdict['ycmp'] = ycmp
    
    # get max year
    cmpdict['ymax'] = int(data[-2].strip())
    
    return cmpdict

'''        
# function to read input parameters used for calculating beta
'''   
def read_betaplfile(betaplfile):
    from numpy import nan
    
    betadict = {}
    data = open(betaplfile,'rb').readlines()
    
    # check if using fixed beta value
    betadict['fixedbeta'] = data[3].strip()
    if betadict['fixedbeta'] == 'Y':
        nlines = 10
        dat = data[4].strip().split(',')
        betadict['fixedbetaval'] = float(dat[0])
        betadict['fixedbetasig'] = float(dat[1])
    else:
        nlines = 9
        betadict['fixedbetaval'] = nan
        betadict['fixedbetasig'] = nan
        
    # get mx lower and upper
    dat = data[nlines-2].strip().split(',')
    betadict['mxlower'] = float(dat[0])
    betadict['mxupper'] = float(dat[1])

    return betadict

'''
# script to parse single gscfrisk *.zon file
def read_zonefile(zonefile):
    # read file
    data = open(zonefile,'rb').readlines()
    
    # get header
    header = data[0].strip('99').split()
    
    zcode = header[0]
    
    zname = ' '.join(header[1:])
    
    
    return zcode, zname
'''    

# reads csv lookup of recurrence model parameters for all zones
# Flag 1: original files
#      2: review files produced using "make_gscfrisk_param_lookup.py": *_review.csv
#      3: review files produced using "make_gscfrisk_review.py": *_modelreview.csv
def read_hazparam_lookup(lufile, vflag):
    hazparams = []
    data = open(lufile).readlines()
    
    inc = 0
    if vflag == 2:
        inc = 1
    
    if vflag <= 2:    
        for i in range(1,len(data)):
            tmpdict = {}
            dat = data[i].strip().split(',')
            tmpdict['zcode'] = dat[0]
            tmpdict['zname'] = dat[1]
            tmpdict['mmin'] = float(dat[2])
            tmpdict['mxbest'] = float(dat[3+inc])
            tmpdict['mxlower'] = float(dat[4+inc])
            tmpdict['mxupper'] = float(dat[5+inc])
            tmpdict['msep'] = float(dat[6+inc])
            tmpdict['fixedbetaval'] = float(dat[7+inc])
            tmpdict['fixedbetasig'] = float(dat[8+inc])
            tmpdict['mcmp'] = map(float,dat[9+inc].split(';'))
            tmpdict['ycmp'] = map(int,dat[10+inc].split(';'))
            tmpdict['ymax'] = int(dat[11+inc])
            try:
                tmpdict['catalogue'] = dat[12+inc]
            except:
                tmpdict['catalogue'] = dat[-1]
            
            hazparams.append(tmpdict)
            
    elif vflag == 3:
        for i in range(1,len(data)):
            tmpdict = {}
            dat = data[i].strip().split(',')
            tmpdict['zcode'] = dat[0]
            tmpdict['mmin'] = float(dat[1])
            tmpdict['mxbest'] = float(dat[2])
            tmpdict['mxlower'] = float(dat[3])
            tmpdict['mxupper'] = float(dat[4])
            tmpdict['gmpe'] = dat[-1]
            
            hazparams.append(tmpdict)
                    
    return hazparams
