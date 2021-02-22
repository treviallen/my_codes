# -*- coding: utf-8 -*-
"""
Created on Thu May 16 12:02:57 2013

functions to parse earthquake catalogues of various formats:

    - parse_sheef(sheeffile)
    
@author: tallen
"""

def getintval(line,start,stop):
    from numpy import nan
    tmpstr = line[start:stop].strip()
    if len(tmpstr) > 0:
        try:
            intval = int(tmpstr)
        except:
            intval = nan
    else:
        intval = nan
        
    return intval
    
def getfloatval(line,start1,stop1,start2,stop2):
    from numpy import isnan, nan, ceil, log10
    
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
    from numpy import nan
    try:
        return float(floatstr)
    except:
        return nan
        
def checkint(intstr):
    from numpy import nan
    try:
        return int(intstr)
    except:
        return nan

def forcezero(num):
    from numpy import isnan
    if isnan(num):
        return 0
    else:
        return num
        
def forceone(num):
    from numpy import isnan
    if isnan(num):
        return 1
    else:
        return num
        
def fix_string_len(string, length):
    # if string too long, truncate
    if len(string) > length:
        return string[0:length]
        
    # else if too short, pad with spaces
    elif len(string) < length:
        for i in range(len(string), length):
            string = string + ' '
        return string
    
    # else correct length    
    else:
        return string

def check_time_fmt(datestr):
    '''
    assumes fmt = '%Y%m%d%H%M'
    '''
    # first check month is non-zero
    if int(datestr[4:6]) == 0:
        datestr = ''.join((datestr[0:4],'01',datestr[-6:]))
    # second check day is non-zero
    if int(datestr[6:8]) == 0:
        datestr = ''.join((datestr[0:6],'01',datestr[-4:]))
        
    return datestr
    

# parse SHHEF-fmt catalogue
def parse_sheef(sheeffile):
    from numpy import nan
    from datetime import datetime as dt
    
    # set dictionary
    sheef = []
    
    # now parse sheef
    print('Reading SHEEF...')
    data = open(sheeffile).readlines()

    for line in data:
        # remove leading space for *.incmp files
        if line[0:2] == '  ':
            line = line[1:]
        
        # "if" is to exclude header info from *.incmp files
        if len(line.strip()) > 60: 
            tmpdict = {}
            # get datetime - check date fmt
            tmptime = check_time_fmt(line[0:13].strip())
            tmpdict['datetime'] = dt.strptime(tmptime, '%Y%m%d%H%M')
            tmpdict['datestr'] = line[0:13].strip() # preserve original time string
            
            tmpdict['locsrc'] = line[40:46].strip()
            tmpdict['year'] = checkint(line[1:5])
            tmpdict['month'] = checkint(line[5:7])
            tmpdict['day'] = checkint(line[7:9])
            tmpdict['hour'] = checkint(line[9:11])
            tmpdict['min'] = checkint(line[11:13])
            tmpdict['lat'] = checkfloat(line[31:37])
            tmpdict['lon'] = checkfloat(line[20:28].strip())
            
            # check depth
            if len(line) < 75:
                dep = checkfloat(line[44:51].strip())
                fixdep = line[52:53].strip()
            else:
                dep = checkfloat(line[47:53].strip())
                fixdep = line[54:55].strip()
            '''    
            if dep == 0.0:
                tmpdict['dep'] = nan
            else:
                tmpdict['dep'] = dep
            '''
            tmpdict['dep'] = dep
                
            if fixdep == 'F' or fixdep == 'G':
                tmpdict['fixdep'] = 1
                tmpdict['depflag'] = fixdep
            else:
                tmpdict['fixdep'] = 0
                tmpdict['depflag'] = ' '
            
            # get mags
            tmpdict['prefmag'] = checkfloat(line[15:18])
            tmpdict['mwconv'] = checkfloat(line[71:])
            
            if len(line) > 56 and len(line) < 58:
                tmpdict['prefmagtype'] = line[54:].strip()
            elif len(line) > 58 and len(line) < 75:
                tmpdict['prefmagtype'] = line[54:57].strip()
            elif len(line) >= 75:
                tmpdict['prefmagtype'] = line[56:59].strip()
            else:
                tmpdict['prefmagtype'] = '  '
                
            if tmpdict['prefmagtype'] == 'Mw' or tmpdict['prefmagtype'] == 'Mwp':
                tmpdict['mwtype'] = tmpdict['prefmagtype']
                tmpdict['mw'] = tmpdict['prefmag']
            else:
                tmpdict['mwtype'] = '  '
                tmpdict['mw'] = nan
            
            if len(line) > 59 and len(line) < 75:
                tmpdict['omag'] = checkfloat(line[58:61])
            elif len(line) >= 75:
                tmpdict['omag'] = checkfloat(line[61:64])
            else:
                tmpdict['omag'] = nan
                
            if len(line) > 63 and len(line) < 67:
                tmpdict['omagtype'] = line[62:].strip()
            elif len(line) > 67 and len(line) < 75:
                tmpdict['omagtype'] = line[62:67].strip()
            elif len(line) >= 75:
                tmpdict['omagtype'] = line[65:70].strip()
            else:
                tmpdict['omagtype'] = ''
                
            sheef.append(tmpdict)
    
    return sheef
    
def get_sheef_field(sheef, key):
    '''
    sheef: dictionary generated by parse_sheef()
    key:   field from above
    dtype: string input ('str', 'float')
    '''
    keydat = []        
    for ev in sheef:
        keydat.append(ev[key])
    
    return keydat
    

def sheef2hmtk_csv(sheeffile):
    
    '''
    returns OQ compliant catalogue in csv fmt
    '''
    
    from io_catalogues import parse_sheef
    from numpy import nan
    
    sheef = parse_sheef(sheeffile)
        
    # make oq cat dict
    header = ','.join(('eventID','year', 'month', 'day', 'hour', 'minute', 'second', 'longitude', 'latitude','depth','magnitude','magnitudeType','Agency'))
    oq_dat = header + '\n'
    # loop thru eqs
    for s in sheef:
        line = ','.join((s['datestr'],str(s['year']), str(s['month']),str(s['day']),str(s['hour']),str(s['min']),str(nan),str(s['lon']),str(s['lat']), \
                        str(s['dep']),str(s['prefmag']),s['prefmagtype'],s['locsrc']))
        oq_dat += line + '\n'
        
    #write to OQ out
    print('Writing HMTK csv...\n')
    hmtkfile = sheeffile.split('.')[0] + '_hmtk.csv'
    f = open(hmtkfile, 'wb')
    f.write(oq_dat)
    f.close()
    
    return hmtkfile
    
def hmtk2sheef(hmtkcat, sheeffile):
    '''
    takes HMTK catalogue class and writes to SHEEF-formatted file
    '''
    print('Writing HMTK to SHEEF...')
    newsheef = ''
    cat = hmtkcat.data #just shorten variable name
    for i in range(0, hmtkcat.get_number_events()):
        mag = str('%0.1f' % cat['magnitude'][i])
        
        lat = '  ' + str('%0.3f' % cat['latitude'][i])
        
        lon = str('%0.3f' % cat['longitude'][i])
        if cat['longitude'][i] > -100.:
            lon = ' ' + lon
        
        dep = str('%0.2f' % cat['depth'][i])
        if isnan(cat['depth'][i]):
            dep = '         '
        elif cat['depth'][i] < 10.:
            dep = '     ' + dep
        elif cat['depth'][i] < 100.:
            dep = '    ' + dep
        else:
            dep = '   ' + dep
        
        src = '  ' + cat['Agency'][i]
        if len(src) == 3:
            src += '  '
        elif len(src) == 3:
            src += ' '
            
        mtype = cat['magnitudeType'][i]
        if len(mtype) == 2:
            mtype += ' '
            
        mag2 = str('%0.2f' % cat['magnitude'][i])
        
        newsheef += ' '.join(('', str(cat['eventID'][i]), '', mag, '', lon, lat, src, dep, ' ', mtype, '', mag, 'UK    ', mag2)) + '\n'
        
    f = open(sheeffile, 'wb')
    f.write(newsheef)
    f.close()
    

'''
write_sheef:  writes sheef formatted catalogue from input dictionary

Usage: write_sheef(sheefdict,outfile)

dict must have following items:
    datestr: fmt yyyymmddHHMM ('%Y%m%d%H%M')
    mw
    lon
    lat
    dep
    depflag
    mwtype (how Mw calculated)
    omag
    omagtype
    locsrc
    mwconv (converted mw)
'''
def write_sheef(sheefdict,outfile):
    from numpy import isnan    
    # make file for writing
    sheeftxt = ''
    
    for rec in sheefdict:
        # set string constant width
        if rec['lon'] > -100.0:
            lon = ' ' + str("%0.3f" % rec['lon'])
        else:
            lon = str("%0.3f" % rec['lon'])
            
        lat = str("%0.3f" % rec['lat'])
          
        locsrc = fix_string_len(rec['locsrc'], 6)
        
        omagtype = fix_string_len(rec['omagtype'], 5)
            
        if rec['dep'] >= 100.0:
            dep = str("%0.2f" % rec['dep'])
        elif rec['dep'] >= 10.0:
            dep = ' ' + str("%0.2f" % rec['dep'])
        elif rec['dep'] >= 0.0:
            dep = '  ' + str("%0.2f" % rec['dep'])
        elif isnan(rec['dep']):
            dep = '      '
            
        if isnan(rec['omag']):
            omag = '   '
        else:
            omag = str("%0.1f" % rec['omag'])
            
        if isnan(rec['mw']):
            prefmag = str("%0.1f" % rec['omag'])
            prefmagtype = fix_string_len(omagtype.strip(), 3)
        else:
            prefmag = str("%0.1f" % rec['mw'])
            prefmagtype = fix_string_len(rec['mwtype'], 3)
            
        '''
        print(' ', rec['datestr'], '  ', prefmag, \
                       '  ', lon, '   ', lat, '   ', locsrc, ' ', dep, ' ', \
                       rec['depflag'], ' ', prefmagtype, ' ', \
                       ' ', omag, ' ', omagtype, '   ', str("%0.2f" % rec['mwconv'])
        '''               
        line = ''.join((' ', rec['datestr'], '  ', prefmag, \
                       '  ', lon, '   ', lat, '   ', locsrc, ' ', dep, ' ', \
                       rec['depflag'], ' ', prefmagtype, ' ', \
                       ' ', omag, ' ', omagtype, '  ', str("%0.2f" % rec['mwconv'])))
                       
        sheeftxt = sheeftxt + line + '\n'
                
    f = open(outfile, 'wb')
    f.write(sheeftxt)
    f.close()

# write sheef-format file to shapefile    
def sheef2shp(sheeffile, shpfile):
    '''
    sheeffile: name of ASCII sheef file
    shpfile:   output shapefile
    '''
    
    import shapefile
    from numpy import nan, isnan
    
    # now parse sheef
    print('Reading SHEEF...')
    data = open(sheeffile).readlines()
    
    # set dictionary
    sheef = []
    
    # now loop through events
    for line in data:
        tmpdict = {}
        tmpdict['locsrc'] = line[40:43].strip()
        tmpdict['year'] = checkint(line[1:5])
        tmpdict['month'] = checkint(line[5:7])
        tmpdict['day'] = checkint(line[7:9])
        tmpdict['hour'] = checkint(line[9:11])
        tmpdict['min'] = checkint(line[11:13])
        tmpdict['lat'] = checkfloat(line[31:37])
        tmpdict['lon'] = checkfloat(line[20:28].strip())
        
        # check depth
        if len(line) < 75:
            dep = checkfloat(line[44:51].strip())
            fixdep = line[52:53].strip()
        else:
            dep = checkfloat(line[47:53].strip())
            fixdep = line[54:55].strip()
            
        if dep == 0.0:
            tmpdict['dep'] = nan
        else:
            tmpdict['dep'] =dep
            
        if fixdep == 'F' or fixdep == 'G':
            tmpdict['fixdep'] = 1
        else:
            tmpdict['fixdep'] = 0
        
        # get mags
        tmpdict['prefmag'] = checkfloat(line[15:18])
        
        if len(line) > 56 and len(line) < 58:
            tmpdict['prefmagtype'] = line[54:].strip()
        elif len(line) > 58 and len(line) < 75:
            tmpdict['prefmagtype'] = line[54:57].strip()
        elif len(line) > 75:
            tmpdict['prefmagtype'] = line[56:59].strip()
        
        if len(line) > 59 and len(line) < 75:
            tmpdict['omag'] = checkfloat(line[58:61])
        elif len(line) > 75:
            tmpdict['omag'] = checkfloat(line[61:64])
        else:
            tmpdict['omag'] = nan
            
        if len(line) > 63 and len(line) < 67:
            tmpdict['omagtype'] = line[62:].strip()
        elif len(line) > 67 and len(line) < 75:
            tmpdict['omagtype'] = line[62:66].strip()
        elif len(line) > 75:
            tmpdict['omagtype'] = line[65:69].strip()
        else:
            tmpdict['omagtype'] = ''
            
        sheef.append(tmpdict)
        
    '''    
    # now make point shapefile
    '''
    
    print('Making shapefile...')
    w = shapefile.Writer(shapefile.POINT)
    w.field('LOC_SRC','C','100')
    w.field('YEAR','F', 9, 0)
    w.field('MONTH','F', 9, 0)
    w.field('DAY','F', 9, 0)
    w.field('HOUR','F', 9, 0)
    w.field('MIN','F', 9, 0)
    w.field('LAT','F', 13, 6)
    w.field('LON','F', 13, 6)
    w.field('DEP','F', 13, 6)
    w.field('FIX_DEP','C','100')
    w.field('PREF_MAG','F', 13, 6)
    w.field('PREF_MAG_TYPE','C','100')
    w.field('ORIG_MAG','F', 13, 6)
    w.field('ORIG_MAG_TYPE','C','100')
    
    # now loop thru records
    for i in range(0, len(sheef)):
        if isnan(sheef[i]['lon']) | isnan(sheef[i]['lat']):
            lon = 0.0
            lat = 0.0
        else:
            lon = sheef[i]['lon']
            lat = sheef[i]['lat']
    
        w.point(lon, lat)
        w.record(sheef[i]['locsrc'],sheef[i]['year'],sheef[i]['month'], \
                 sheef[i]['day'],sheef[i]['hour'],sheef[i]['min'], \
                 sheef[i]['lat'],sheef[i]['lon'], \
                 sheef[i]['dep'],sheef[i]['fixdep'], \
                 sheef[i]['prefmag'],sheef[i]['prefmagtype'], \
                 sheef[i]['omag'],sheef[i]['omagtype'])
    
    print('Writing shapefile...')
    w.save(shpfile)    

# code to parse CEEF and dump to pickle - don't think this works, but preserving just in case 2015-02-26
def parse_ceef(ceeffile, pklfile):
    '''
    ceeffile is the eqlist.out.XXXX.XXX file
    '''        
    from pickle import dump
    from datetime import datetime
    from numpy import floor
    
    # now parse ceef
    print('Reading CEEF...')
    data = open(ceeffile).readlines()
    
    # set dictionary
    ceefdat = []
    
    # now loop through events
    readevent = False
    for line in data:
        if readevent == True:
            tmpdict = {}
            tmpdict['locsrc'] = 'GSC'
            tmpdict['year'] = checkint(line[0:4])
            tmpdict['month'] = checkint(line[5:7])
            tmpdict['day'] = checkint(line[8:10])
            tmpdict['hour'] = checkint(line[11:13])
            tmpdict['min'] = checkint(line[13:15])
            tmpdict['sec'] = checkfloat(line[16:20])
            tmpdict['lat'] = checkfloat(line[25:31])
            tmpdict['lon'] = checkfloat(line[32:40])
            tmpdict['dep'] = checkfloat(line[42:47])
            tmpdict['depflag'] = line[47]
            tmpdict['omag'] = checkfloat(line[53:57])
            tmpdict['mw'] = checkfloat(line[59:])
            
            # make datetime object
            tmpmonth = forceone(tmpdict['month'])
            tmpday   = forceone(tmpdict['day'])
            tmphour  = forcezero(tmpdict['hour'])
            tmpmin   = forcezero(tmpdict['min'])
            tmpsec   = forcezero(floor(tmpdict['sec']))
            tmpmsec  = 1000000 * (forcezero(tmpdict['sec']) - tmpsec)
            
            tmpdict['datetime'] = datetime(tmpdict['year'], tmpmonth, tmpday, \
                                           tmphour, tmpmin, int(tmpsec), int(tmpmsec))
            
            ceefdat.append(tmpdict)
        
        # begin reading events    
        elif line.find('*****') == 0:
            readevent = True
    
    # pickle ceef
    #print(ceefdat[-1]
    output = open(pklfile, 'wb')
    dump(ceefdat,output)
    output.close()
    
    return ceefdat

# code to parse CEEF and dump to pickle  
def parse_ceef_out(ceeffile):
    '''
    ceeffile is the eqlist.out.XXXX.XXX file
    '''
    from datetime import datetime
    from numpy import floor, isnan
    
    # now parse ceef
    print('Reading CEEF...')
    data = open(ceeffile).readlines()
    
    # set dictionary
    ceefdat = []
    
    # now loop through events
    readevent = False
    for i, line in enumerate(data):
        if readevent == True:
            stopnonpref = False
            if line[0] == '+': # preferred solution          
                tmpdict = {}
                tmpdict['locsrc'] = line[111:].rstrip('\n\r')
                tmpdict['year'] = checkint(line[2:6])
                tmpdict['month'] = checkint(line[7:9])
                tmpdict['day'] = checkint(line[10:12])
                tmpdict['hour'] = checkint(line[13:15])
                tmpdict['min'] = checkint(line[15:17])
                tmpdict['sec'] = checkfloat(line[18:22])
                tmpdict['lat'] = checkfloat(line[27:33])
                tmpdict['lon'] = checkfloat(line[34:42])
                tmpdict['dep'] = checkfloat(line[44:49])
                tmpdict['depflag'] = line[49]
                tmpdict['mb'] = checkfloat(line[55:59])
                tmpdict['mn'] = checkfloat(line[61:65])
                tmpdict['ml'] = checkfloat(line[67:71])
                tmpdict['ms'] = checkfloat(line[73:77])
                tmpdict['mc'] = checkfloat(line[79:83])
                tmpdict['mw'] = checkfloat(line[85:89])
                tmpdict['otm'] = checkfloat(line[91:95])
                tmpdict['mpref'] = checkfloat(line[97:101])
                tmpdict['mwp'] = checkfloat(line[103:107])
                
                # populate alternative magnitudes from non-pref locations
                j = 1                
                while stopnonpref == False:
                    if i+j < len(data):
                        line2 = data[i+j]

                        if line2[0] == '+':
                            stopnonpref = True
                        else: # fill mags
                            if isnan(tmpdict['mb']):
                                tmpdict['mb'] = checkfloat(line2[55:59])
                            if isnan(tmpdict['mn']):
                                tmpdict['mn'] = checkfloat(line2[61:65])
                            if isnan(tmpdict['ml']):
                                tmpdict['ml'] = checkfloat(line2[67:71])
                            if isnan(tmpdict['ms']):
                                tmpdict['ms'] = checkfloat(line2[73:77])
                            if isnan(tmpdict['mc']):
                                tmpdict['mc'] = checkfloat(line2[79:83])
                            if isnan(tmpdict['mw']):
                                tmpdict['mw'] = checkfloat(line2[85:89])
                            if isnan(tmpdict['otm']):
                                tmpdict['otm'] = checkfloat(line2[91:95])
                        j+=1
                    else:
                        stopnonpref = True

                    
                # get preferred mag type
                if tmpdict['mw'] == tmpdict['mpref']:
                    tmpdict['mpreftype'] = 'MW'
                elif tmpdict['ml'] == tmpdict['mpref']:
                    tmpdict['mpreftype'] = 'ML'
                elif tmpdict['mb'] == tmpdict['mpref']:
                    tmpdict['mpreftype'] = 'mb'
                elif tmpdict['ms'] == tmpdict['mpref']:
                    tmpdict['mpreftype'] = 'MS'
                elif tmpdict['mn'] == tmpdict['mpref']:
                    tmpdict['mpreftype'] = 'Mn'
                elif tmpdict['mc'] == tmpdict['mpref']:
                    tmpdict['mpreftype'] = 'MC'
                elif tmpdict['otm'] == tmpdict['mpref']:
                    tmpdict['mpreftype'] = 'OT'
                elif tmpdict['mwp'] == tmpdict['mpref']:
                    tmpdict['mpreftype'] = 'MWp'
                else:
                    tmpdict['mpreftype'] = ''
                                
                # make datetime object
                tmpmonth = forceone(tmpdict['month'])
                tmpday   = forceone(tmpdict['day'])
                tmphour  = forcezero(tmpdict['hour'])
                tmpmin   = forcezero(tmpdict['min'])
                tmpsec   = forcezero(floor(tmpdict['sec']))
                tmpmsec  = 1000000 * (forcezero(tmpdict['sec']) - tmpsec)
                
                tmpdict['datetime'] = datetime(tmpdict['year'], tmpmonth, tmpday, \
                                               tmphour, tmpmin, int(tmpsec), int(tmpmsec))
                #print(tmpdict
                ceefdat.append(tmpdict)
        
        # begin reading events    
        elif line.find('*****') == 0:
            readevent = True
            
    return ceefdat

# scrape Earthquakes of the last 30 days webpage
def latest_events2shp(outshp):  
    """
    Script to parse latest events webpage and export to shapefile    
    Usage:
        python latest_events2shp.py outshp
        
        where:
            outshp is output file name (e.g. latest_events.shp)
    """    
    from lxml import html
    import requests
    import shapefile
    
    try:
        
        # parse webpage
        page = requests.get('http://www.earthquakescanada.nrcan.gc.ca/recent/maps-cartes/index-eng.php')
        tree = html.fromstring(page.text)
        
        # get event data
        date = tree.xpath('//td[@headers="date"]/text()')
        time = tree.xpath('//td[@headers="time"]/text()')
        lat = tree.xpath('//td[@headers="lat"]/text()')
        lon = tree.xpath('//td[@headers="lon"]/text()')
        dep = tree.xpath('//td[@headers="depth"]/text()')
        mag = tree.xpath('//td[@headers="mag"]/text()')
        #felt = tree.xpath('//span[@style="color"]/text()')
        reg = tree.xpath('//td[@headers="region"]/text()')
        
        # convert relevant fields to float
        lat = [float(x) for x in lat]
        lon = [float(x) for x in lon]
        dep = [float(x) for x in dep]
        mag = [float(x) for x in mag]
        
        # set shapefile headers
        w = shapefile.Writer(shapefile.POINT)
        w.field('DATE','C','12')
        w.field('TIME','C', '10')
        w.field('LAT','F', 13, 6)
        w.field('LON','F', 13, 6)
        w.field('DEP','F', 13, 6)
        w.field('MAG','F', 13, 6)
        #w.field('FELT','C','10')
        w.field('REGION','C','100')
        
        # now loop thru records and add to shapefile
        for i in range(0, len(date)):
            w.point(lon[i], lat[i])
            w.record(date[i], time[i], lat[i], lon[i], dep[i], mag[i], reg[i])
        
        print('Writing shapefile...')
        w.save(outshp)
        
        # write projection file
        prjfile = outshp.strip().split('.')[0]+'.prj'
        f = open(prjfile, 'wb')
        f.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]')
        f.close()
        
    except:
        print('\nUsage: python latest_events2shp.py outshp\n')

# parses GGCat csvfile
def parse_ggcat(ggcatcsv):
    import csv
    from numpy import nan, isnan, floor
    from misc_tools import checkint
    import datetime as dt
    
    # open file
    raw = open(ggcatcsv).readlines()[1:] # exclude header    
    
    # parse csv    
    lines = csv.reader(raw)
    
    # set array to append event dictionaries
    ggcat = []
    
    # loop through events
    for line in lines:
        # set null vals to nans
        for i in range(len(line)):
           if len(str(line[i]))<1:
              line[i] = nan
        
        # fill temp dict
        tmpdict = {'auth':line[0], 'place':line[1],'year':checkint(line[5]), 'month':checkint(line[6]), 'day':checkint(line[7]), \
                   'hour':checkint(line[8]), 'min':checkint(line[9]), 'sec':float(line[10]), 'lon':float(line[11]), 'lat':float(line[12]), 'dep':float(line[13]), \
                   'zcode':line[14], 'prefmagtype':line[15], 'prefmag':float(line[16]), 'ml':float(line[17]), 'mb':float(line[18]), 'ms':float(line[19]), \
                   'mw':float(line[20]), 'md':float(line[21]), 'mp':float(line[22]), 'fixdep':0}
                   	
        # add datetime
        if ~isnan(tmpdict['sec']):
            if int(floor(tmpdict['sec'])) >= 60:
                tmpdict['sec'] = 59
            tmpdict['datetime'] = dt.datetime(tmpdict['year'], tmpdict['month'], tmpdict['day'], tmpdict['hour'], tmpdict['min'], int(floor(tmpdict['sec'])))
        elif ~isnan(tmpdict['hour']):
            tmpdict['datetime'] = dt.datetime(tmpdict['year'], tmpdict['month'], tmpdict['day'], tmpdict['hour'], tmpdict['min'], 0)
        else:
            tmpdict['datetime'] = dt.datetime(tmpdict['year'], tmpdict['month'], tmpdict['day'])

        
        ggcat.append(tmpdict)
        
    return ggcat

# parses USGS format web query csv
def parse_usgs_event_query(usgscsv):
    from obspy.core.utcdatetime import UTCDateTime
    lines = open(usgscsv).readlines()[1:]
    
    # build dict
    evdict = []
    for line in lines:
        dat = line.strip().split(',')
        
        #evdt = dt.datetime.strptime(dat[0], '%Y-%m-%dT%H:%M:%S.%fZ')
        starttime = dat[0][:15] + '0:00.000Z'
        tdict = {'time': UTCDateTime(dat[0]), 'lat': float(dat[1]), 'lon': float(dat[2]), \
                 'dep': float(dat[3]), 'mag': float(dat[4]), 'magType': dat[5], \
                 'timestr': dat[0], 'starttime': UTCDateTime(starttime)}
                 
        evdict.append(tdict)
        
    return evdict

# parses csv export from GA web search - needs updating
def parse_ga_event_query(gacsv):
    '''
    gacsv: path to csv file downloaded from https://earthquakes.ga.gov.au/
    '''
    import datetime as dt
    
    # open and read file
    lines = open(gacsv).readlines()[1:]
    
    # initiate event dictonary
    evdict = []
    
    for line in lines:
        dat = line.strip().split(',')
        
        if len(dat) >= 27 and not line.startswith('"'):
            try:
                try:
                    dateTime = dt.datetime.strptime(dat[27], '%Y-%m-%dT%H:%M:%S.%f')
                except:
                    try:
                        dateTime = dt.datetime.strptime(dat[27], '%Y-%m-%dT%H:%M:%S')
                    except:
                        continue
                oidx = 27
                
            except:
                oidx = 24
            #print oidx
            try:
                dateTime = dt.datetime.strptime(dat[oidx], '%Y-%m-%dT%H:%M:%S.%f')
            except:
                dateTime = dt.datetime.strptime(dat[oidx], '%Y-%m-%dT%H:%M:%S')
            #print(dateTime)
            
            # note - changed pref mag idx from 27 to 30
            tdict = {'datetime': dateTime, \
                     'year': dateTime.year, 'month': dateTime.month, \
                     'day': dateTime.day, 'hour': dateTime.hour, \
                     'minute': dateTime.minute, 'second': dateTime.second, \
                     'lat': float(dat[13]), 'lon': float(dat[14]), \
                     'dep': float(dat[4]), 'mag_ml': checkfloat(dat[17]), \
                     'mag_mb': checkfloat(dat[15]), 'mag_ms': checkfloat(dat[18]),
                     'mag_mw': checkfloat(dat[19]), 'mag_mwp': checkfloat(dat[21]),
                     'mag': checkfloat(dat[30]), 'magType': dat[28], \
                     'timestr': dat[10], 'description': dat[6], 'event_id':dat[11]} 
                     	
            evdict.append(tdict)
            
    return evdict

