# parse ShakeMap info xml
def parse_infoxml(infoxml):
    import xml.etree.ElementTree as e

    # read file
    tree = e.parse(infoxml)
    root = tree.getroot()

    # get data
    infodat = {}
    for child in root.iter('tag'):
        key = child.attrib['name']
        tmptxt = child.attrib['value'].strip().split('::')[-1]
        infodat[key] = tmptxt

    return infodat

# parse ShakeMap event XML file
def parse_eventxml(evxml):
    import xml.etree.ElementTree as e

    # read file
    tree = e.parse(evxml)
    root = tree.getroot()

    # get event info
    event = {}
    event['id'] = root.attrib['id']
    event['lat'] = root.attrib['lat']
    event['lon'] = root.attrib['lon']
    event['dep'] = root.attrib['depth']
    event['mag'] = root.attrib['mag']
    event['year'] = root.attrib['year']
    event['month'] = root.attrib['month']
    event['day'] = root.attrib['day']
    event['hour'] = root.attrib['hour']
    event['min'] = root.attrib['minute']
    event['sec'] = root.attrib['second']
    event['locstring'] = root.attrib['locstring']

    return event

# parse grid.xml and return event details, grid specs and fields and data grid
def parse_dataxml(datxml):
    import xml.etree.ElementTree as e
    from numpy import fromstring
    
    # read file
    tree = e.parse(datxml)
    root = tree.getroot()
    
    # get event details
    for child in root.iter('{http://earthquake.usgs.gov/eqcenter/shakemap}event'):
        event = {'event_id': child.attrib['event_id'], 
                 'event_network': child.attrib['event_network'], 
                 'event_timestamp': child.attrib['event_timestamp'], 
                 'event_description':  child.attrib['event_description'], 
                 'lon': float(child.attrib['lon']), 
                 'lat': float(child.attrib['lat']), 
                 'depth': float(child.attrib['depth']), 
                 'magnitude': float(child.attrib['magnitude'])}
        
    # get grid specs
    for child in root.iter('{http://earthquake.usgs.gov/eqcenter/shakemap}grid_specification'):
        gridspec = {'lon_min': float(child.attrib['lon_min']), 
                 'lat_min': float(child.attrib['lat_min']), 
                 'lon_max': float(child.attrib['lon_max']), 
                 'lat_max': float(child.attrib['lat_max']), 
                 'nominal_lon_spacing': float(child.attrib['nominal_lon_spacing']), 
                 'nominal_lat_spacing': float(child.attrib['nominal_lat_spacing']), 
                 'nlon': int(child.attrib['nlon']),
                 'nlat': int(child.attrib['nlat'])}
    
    # get fields
    fields = {}
    for i, child in enumerate(root.iter('{http://earthquake.usgs.gov/eqcenter/shakemap}grid_field')):
        column = 'field'+child.attrib['index']
        units  = 'units'+child.attrib['index']
        fields[column] = child.attrib['name'].lower()
        fields[units]  = child.attrib['units'].lower()
    fields['cols'] = i+1
    
    # now get data    
    for child in tree.iter('{http://earthquake.usgs.gov/eqcenter/shakemap}grid_data'):
        dataStr = child.text
        
    # reformat text to data array
    griddata = fromstring(dataStr.strip(), sep=' ')
    rows = len(griddata) / fields['cols']
    newshape = (rows, fields['cols'])
    griddata = griddata.reshape(newshape)

    return event, gridspec, fields, griddata

# parse data dat xml
def parse_gridxml(datxml):
    import xml.etree.ElementTree as e

    # read file
    tree = e.parse(datxml)
    root = tree.getroot()
    stndict = []

    for child in root.iter('station'):
        tmpdict = {}
        tmpdict['insttype'] = child.attrib['insttype']
        tmpdict['commtype'] = child.attrib['commtype']
        tmpdict['code'] = child.attrib['code']
        #tmpdict['dist'] = child.attrib['dist']
        tmpdict['name'] = child.attrib['name']
        tmpdict['netid'] = child.attrib['netid']
        tmpdict['lon'] = child.attrib['lon']
        tmpdict['lat'] = child.attrib['lat']
        tmpdict['source'] = child.attrib['source']
        #tmpdict['loc'] = child.attrib['loc']
        # check if intensity exists
        for key in child.keys():
            if key == 'intensity':
                tmpdict['intensity'] = child.attrib['intensity']

                # only add if intensity exists
                stndict.append(tmpdict)

        '''
        # note, this works - but not needed
        for subchild in child.iter('code'):
            print subchild.attrib
        '''
    return stndict

def parse_infojson(jsonFilePath):
    import json
    from numpy import array
    
    jsonFilePath = 'info.json'
    with open(jsonFilePath) as f:
        data = json.load(f)
    
    # return imts
    imts = data['output']['ground_motions'].viewkeys()
    
    # loop thru imts
    bias = []
    max_land = []
    max_grid = []
    units = []
    for imt in imts:
        # get values
        bias.append(data['output']['ground_motions'][imt]['bias'])
        max_land.append(data['output']['ground_motions'][imt]['max'])
        max_grid.append(data['output']['ground_motions'][imt]['max_grid'])
        units.append(data['output']['ground_motions'][imt]['units'])
        
    bias = array(bias)
    max_land = array(max_land)
    max_grid = array(max_grid)
    units = array(units)
    
    return data, bias, max_land, max_grid, units

def station_xml2csv(stationxml, stationcsv):
    from shakemap_tools import parse_dataxml
    
    stadict = parse_dataxml(stationxml)
    
    csvtxt = 'CODE,LON,LAT,MMI,SOURCE\n'
    for sta in stadict:
        csvtxt += ','.join((sta['code'], str('%0.4f' % float(sta['lon'])), str('%0.4f' % float(sta['lat'])), \
                           sta['intensity'], sta['source'])) + '\n'
    
    f = open(stationcsv, 'wb')
    f.write(csvtxt)
    f.close()
        
# parse fault files
def parse_faultdat(faultfile):
    from numpy import nan

    txt = open(faultfile).readlines()
    lat = []
    lon = []
    dep = []
    faultdat = []
    for line in txt:
        if line[0].strip() != '>' and len(line[0].strip()) > 0:
            
            lat.append(float(line.strip().split()[1])) 
            lon.append(float(line.strip().split()[0]))
            if len(line.strip().split()) > 2:
                dep.append(float(line.strip().split()[2]))
            else:
                dep.append(nan)

        else:
            if len(lat) > 0:
                tmpfault = {'lat': lat, 'lon': lon, 'dep': dep}
                faultdat .append(tmpfault)
            lat = []
            lon = []
            dep = []

    if len(lat) > 0:
        tmpfault = {'lat': lat, 'lon': lon, 'dep': dep}
        faultdat .append(tmpfault)

    return faultdat

# make fault mesh
'''
makes fault mesh from shakemap *_fault.txt files
	res is grid resolution in km
'''
def make_fault_mesh(faultfile, res):
    from mapping_tools import get_line_parallels, distance, reckon
    from numpy import array, arcsin, argsort, argwhere, abs, radians, arange, ceil, sqrt, unique, tan
    import matplotlib.pyplot as plt
    #from mpl_toolkits.mplot3d import Axes3D
    #fig = plt.figure()
    #ax = Axes3D(fig)
    #ax.view_init(60, 210)
    #
    faultdat = parse_faultdat(faultfile)

    # loop thru sub-faults
    x = []
    y = []
    z = []
    for sf in faultdat:
        lat = sf['lat']
        lon = sf['lon']
        dep = sf['dep']

        #ax.scatter(lon, lat, -1*array(dep), marker='o', color = 'b')

    #    unidep = unique(dep) # not used, but will be useful


        # check fault has vertical dimensions
        # if not, do vertically dipping fault
        if max(dep) == min(dep):
            for i in range(0, len(lat)-1):
                rngkm, az, baz = distance(lat[i], lon[i], lat[i+1], lon[i+1])
                # get npts
                npts = int(ceil(rngkm / res))

                for j in range(0, npts+1):
                    hinc = j * rngkm / npts
                    xy = reckon(lat[i], lon[i], hinc, az)
                    x.append(xy[0])
                    y.append(xy[1])
                    z.append(dep[0])

        # do dipping fault
        else:
            lat = lat[0:-1]
            lon = lon[0:-1]
            dep = dep[0:-1]

            # first find shallow and deep coords
            dsort = argsort(dep)

            # get side 1 dist & az - assumes top deps the same
            sidekm1, az1, baz1 = distance(lat[dsort[0]], lon[dsort[0]], lat[dsort[0]+3], lon[dsort[0]+3])
            ddep1 = float(abs(dep[dsort[0]] - dep[dsort[0]+3]))
            abskm1 = sqrt(sidekm1**2 + ddep1**2)
            dip1 = arcsin(ddep1 / abskm1)
            vpts = int(ceil(abskm1 / res))
            hinc1 = sidekm1 / vpts
            vinc1 = sqrt(abskm1**2 - sidekm1**2) / vpts

            # get side 2 dist & az - assumes top deps the same
            sidekm2, az2, baz2 = distance(lat[dsort[1]], lon[dsort[1]], lat[dsort[1]+1], lon[dsort[1]+1])
            ddep2 = float(abs(dep[dsort[1]] - dep[dsort[1]+1]))
            abskm2 = sqrt(sidekm2**2 + ddep2**2)
            dip2 = arcsin(ddep2 / abskm2)
            hinc2 = sidekm2 / vpts
            vinc2 = sqrt(abskm2**2 - sidekm2**2) / vpts


            for i in range(0, vpts+1):
                # get down-dip xyz locs
                xy1 = reckon(lat[dsort[0]], lon[dsort[0]], i*hinc1, az1)
                xy2 = reckon(lat[dsort[1]], lon[dsort[1]], i*hinc2, az2)
                z1  = dep[dsort[0]] + i*vinc1
                z2  = dep[dsort[1]] + i*vinc2

                # get horiz range
                rngkm, az, baz = distance(xy1[1], xy1[0], xy2[1], xy2[0])
                abskm = sqrt(rngkm**2 + abs(z1 - z2)**2)
                # get angle between dep0 and dep1
                ddep = float(abs(z1 - z2))
                # get horiz pts
                hpts = int(ceil(abskm / res))

                # now get pts joining
                for j in range(0, hpts+1):
                    hinc = j * rngkm / hpts
                    ainc = j * abskm / hpts
                    xy = reckon(xy1[1], xy1[0], hinc, az)
                    x.append(xy[0])
                    y.append(xy[1])
                    z.append(-1*(sqrt(ainc**2 - hinc**2) + z1))


    #ax.scatter(x, y, z, marker='+', color = 'r')
    #plt.show()
    return x, y, z

def change_grind_param(smpath, configin, event, param, value, numflags):
    '''
    changes all GMPEs in output config path to input
    event file must exist

        smpath = '/usr/local/shake'
        event = <evid>
        param = 'gmpe', 'ipe'
        value = e.g. 'BJF97'
        numflags = e.g. 4, if negative, no flags are used
        configin: 0 = use default from shake/config; 1 = use shake/data/<evid>/config
    '''
    from os import path, mkdir

    # get grind config file
    if configin == 0:
        grindfile = path.join(smpath, 'config', 'grind.conf')
    elif configin == 1:
        grindfile = path.join(smpath, 'data', event, 'config', 'grind.conf')

    # get output file path
    grindout = path.join(smpath, 'data', event, 'config', 'grind.conf')
    configdir = path.join(smpath, 'data', event, 'config')
    try:
        f = open(grindout,'r')
    except:
        try:
            mkdir(configdir)
        except:
            'Folder exists'

    # now read grind.config and change gmpe
    outtxt = ''
    txt = open(grindfile).readlines()
    for line in txt:
        # initalise newline
        newline = line
        dat = line.split()
        if len(dat) > 0:
            if len(dat[0]) >= len(param):
                if dat[0][0:len(param)] == param:
                    if numflags > 0:
                        newline = ' '.join((param.strip('#'), ':', value, ' '.join(dat[-numflags:]), '\n'))
                    else:
                        newline = ' '.join((param.strip('#'), ':', value, '\n'))
        outtxt += newline

    #write to file
    f = open(grindout, 'wb')
    f.write(outtxt)
    f.close()

def write_mmi_obs_raw(source_str, mmidict):
    obstxt = 'Source: '+source_str+'\n'
    
    for mmi in mmidict:
        if mmi['nresp'] > 1: 
            try:
                obstxt += '\t'.join((str(mmi['cdi']), str(mmi['lat']), str(mmi['lon']), mmi['loc'])) + '\n'
            except:
                obstxt += '\t'.join((str(mmi['cdi']), str(mmi['lat']), str(mmi['lon']), 'Unknown')) + '\n'
            	                                                                            
    # write SM raw            
    outfile = 'mmi_obs.dat'
    f = open(outfile, 'wb')   
    f.write(obstxt)           
    f.close()     

# this needs work!!!!
def make_shakemap_stationlist_xml():
    import pickle
    from sys import argv
    from numpy import interp
    
    locstring = argv[1]
    
    # get station data
    cnsntxt = open('canstn20140116.txt','rb').readlines()
    cnsn = []
    for line in cnsntxt:
        tmp = {}
        tmp['stn'] = line.strip().split('|')[1].strip()
        tmp['lon'] = line.strip().split('|')[4].strip()
        tmp['lat'] = line.strip().split('|')[3].strip()
        tmp['netid'] = line.strip().split('|')[6].strip()
        tmp['place'] = line.strip().split('|')[8].strip()
        cnsn.append(tmp)
    
    # open pkl file
    pklfile = open('20140424_2.pkl', 'rb')
    
    # load data
    gmdat = pickle.load(pklfile)
    gmdat = gmdat[0:54] # somehting a bit wonky here!!!!
    
    # get eq details
    eqid  = gmdat[0]['eqtime'].strftime('%Y%m%d%H%M')
    eqlat = str(gmdat[0]['eqlat'])
    eqlon = str(gmdat[0]['eqlon'])
    mag   = str(gmdat[0]['Mw'])
    year  = gmdat[0]['eqtime'].strftime('%Y')
    month  = gmdat[0]['eqtime'].strftime('%m')
    day  = gmdat[0]['eqtime'].strftime('%d')
    hour  = gmdat[0]['eqtime'].strftime('%H')
    minute  = gmdat[0]['eqtime'].strftime('%M')
    sec  = gmdat[0]['eqtime'].strftime('%S')
    dep = str(gmdat[0]['dep'])
    
    # make earthquake header
    smtxt = '<shakemap-data code_version="3.5.1 GSM" map_version="1">\n'
    earthquake = ''.join(('<earthquake id="', eqid, '" lat="', eqlat, '" lon="', \
                           eqlon, '" mag="', mag, '" year="', year, '" month="', month, \
                           '" day="', day, '" hour="', hour, '" minute="', minute, \
                           '" second="', sec, '" timezone="GMT" depth="', dep, \
                           '" locstring="', locstring, '" created="1313677296"/>\n'))
    
    endtxt = '</stationlist>\n</shakemap-data>'
    
    smtxt += earthquake + '<stationlist created="1313677296">\n'
    
    # now loop thru station data and get GM values
    cnt = 0
    for stn in gmdat:
    
        # get stn header
        code = stn['stn']
    
        insttype = stn['zchan'][0:2]
    
        # get stn loc
        for s in cnsn:
            if s['stn'].strip() == stn['stn']:
                stlat = s['lat']
                stlon = s['lon']
    
                source = 'CNSN'
                netid  = s['netid']
                name = s['place']
    
                station = '"'.join(('<station code=', code, ' name=', name, ' insttype=', insttype, \
                                    ' lat=', stlat, ' lon=', stlon, ' source=', source, \
                                    ' netid=', netid, ' commtype="UNK">\n'))
    
                smtxt += station
    
                # make Z comp
                name = stn['zchan']
                acc = str('%0.4f' % (stn['zpga'] * 100.))
                vel = str('%0.4f' % (stn['zpgv']))
                psa03 = str('%0.4f' % (interp(0.3, stn['periods'], stn['zpsa'].flatten()) * 100.))
                psa10 = str('%0.4f' % (interp(1.0, stn['periods'], stn['zpsa'].flatten()) * 100.))
                psa30 = str('%0.4f' % (interp(3.0, stn['periods'], stn['zpsa'].flatten()) * 100.))
    
                zcomp = '"'.join(('<comp name=', name, '>\n<acc value=', acc, '/>\n<vel value=', \
                                  vel, '/>\n<psa03 value=', psa03, '/>\n<psa10 value=', psa10, \
                                   '/>\n<psa30 value=', psa30, '/>\n</comp>\n'))
    
                try:
                   # make E comp
                   name = stn['echan']
                   acc = str('%0.4f' % (stn['epga'] * 100.))
                   vel = str('%0.4f' % (stn['epgv']))
                   psa03 = str('%0.4f' % (interp(0.3, stn['periods'], stn['epsa'].flatten()) * 100.))
                   psa10 = str('%0.4f' % (interp(1.0, stn['periods'], stn['epsa'].flatten()) * 100.))
                   psa30 = str('%0.4f' % (interp(3.0, stn['periods'], stn['epsa'].flatten()) * 100.))
    
                   ecomp = '"'.join(('<comp name=', name, '>\n<acc value=', acc, '/>\n<vel value=', \
                                     vel, '/>\n<psa03 value=', psa03, '/>\n<psa10 value=', psa10, \
                                      '/>\n<psa30 value=', psa30, '/>\n</comp>\n'))
    
                   # make E comp
                   name = stn['nchan']
                   acc = str('%0.4f' % (stn['npga'] * 100.))
                   vel = str('%0.4f' % (stn['npgv']))
                   psa03 = str('%0.4f' % (interp(0.3, stn['periods'], stn['npsa'].flatten()) * 100.))
                   psa10 = str('%0.4f' % (interp(1.0, stn['periods'], stn['npsa'].flatten()) * 100.))
                   psa30 = str('%0.4f' % (interp(3.0, stn['periods'], stn['npsa'].flatten()) * 100.))
    
                   ncomp = '"'.join(('<comp name=', name, '>\n<acc value=', acc, '/>\n<vel value=', \
                                     vel, '/>\n<psa03 value=', psa03, '/>\n<psa10 value=', psa10, \
                                      '/>\n<psa30 value=', psa30, '/>\n</comp>\n'))
    
                   smtxt += ecomp + ncomp + zcomp + '</station>\n'
    
                except:
                   smtxt += zcomp + '</station>\n'
    
    # end text
    smtxt += endtxt
    
    # write to text
    f = open('cnsn_dat.xml', 'wb')
    f.write(smtxt)
    f.close()
    
            
# this needs work!!!!
def csv2stationlist_xml(csvfile, eqla, eqlo, dep, mag, yyyymmddHHMMSS, locstring):
    from numpy import interp
    from roman import toRoman
    from mapping_tools import distance
    from mmi_tools import mmi2pgm_worden12, cmpsps2g
    import time
    
    print yyyymmddHHMMSS
    '''
    csvfile format with one header line:
        LAT, LON, MMI
    '''
    
    # get station data
    lines = open(csvfile).readlines()[1:]
    lat = []
    lon = []
    mmi = []
    for line in lines:
        dat = line.strip().split(',')
        lat.append(dat[0]) # keep as string
        lon.append(dat[1]) # keep as string
        mmi.append(float(dat[2]))
    
    # make earthquake header
    smtxt = '<shakemap-data code_version="4.0 GSM" map_version="1">\n'
    earthquake = ''.join(('<earthquake id="', str(yyyymmddHHMMSS), '" lat="', str(eqla), '" lon="', \
                           str(eqlo), '" mag="', str(mag), '" year="', str(yyyymmddHHMMSS)[0:4], '" month="', str(yyyymmddHHMMSS)[4:6], \
                           '" day="', str(yyyymmddHHMMSS)[6:8], '" hour="', str(yyyymmddHHMMSS)[8:10], '" minute="', str(yyyymmddHHMMSS)[10:12], \
                           '" second="', str(yyyymmddHHMMSS)[12:], '" timezone="GMT" depth="', str(dep), \
                           '" locstring="', locstring, '" created="', str(int(time.time())), '"/>\n'))
    
    endtxt = '</stationlist>\n</shakemap-data>'
    
    smtxt += earthquake + '<stationlist created="' + str(int(time.time())) + '">\n'
    
    # now loop thru station data and get GM values
    i = 0
    # get stn loc
    for la, lo, mi in zip(lat, lon, mmi):
        code = 'OBS_'+str(i)
        source = 'AU'
        netid  = 'Intensity'
        name = 'Unknown (Intensity ' + toRoman(round(mi)) +')'
        
        station = '"'.join(('<station code=', code, ' name=', name, ' insttype="Observed"', \
                            ' lat=', la, ' lon=', lo, ' source=', source, \
                            ' netid=', netid, ' commtype="UNK">\n'))

        smtxt += station

        # get distance
        rngkm = distance(float(la), float(lo), float(eqla), float)        
        
        # make mmi comp
        pga = cmpsps2g(mmi2pgm_worden12(mi, 'pga', float(mag), rngkm)[0])
        pgv = (mmi2pgm_worden12(mi, 'pgv', mag, rngkm)[0])
        
        acc = str('%0.4f' % pga)
        vel = str('%0.4f' % pgv)
        smtxt += '"'.join(('<    comp name="DERIVED">\n        <acc value=', acc, '/>\n        <vel value=', \
                          vel, '/>\n    </comp>\n</station>\n'))

           
        i += 1

    # end text
    smtxt += endtxt
    
    # write to text
    f = open(str(yyyymmddHHMMSS)+'_dat.xml', 'wb')
    f.write(smtxt)
    f.close()
    
            
                                      
                          