def fix_stream_channels(mseedfile, accel=False):
    from obspy import read
    
    st = read(mseedfile)
    
    # loop thru traces
    for tr in st:
        if accel == False:
            
            if tr.stats['channel'].startswith('EL') or tr.stats['channel'].startswith('EY') \
            or tr.stats['channel'].startswith('DH') or tr.stats['channel'].startswith('DD') \
            or tr.stats['channel'].startswith('HY') or tr.stats['channel'].startswith('SY'):
            
                if tr.stats['station'].startswith('SWN1') or tr.stats['station'].startswith('SWN0') \
                    or tr.stats['station'].startswith('SWN2'):
                    tr.stats['channel'] = 'HH' + tr.stats['channel'][-1]
                else:
                    tr.stats['channel'] = 'EH' + tr.stats['channel'][-1]
        
        
        elif accel == True:
            tr.stats['channel'] = 'HN' + tr.stats['channel'][-1]
        '''
        elif tr.stats['channel'].startswith('EN') or tr.stats['channel'].startswith('DN'):
            tr.stats['channel'] = 'HN' + tr.stats['channel'][-1]
        '''
    
    # overwrite mseed
    st.write(mseedfile, format="MSEED")
    
def fix_swan_stream_channels(mseedfile):
    from obspy import read
    
    st = read(mseedfile)
    
    # loop thru traces
    for tr in st:
    
        if tr.stats['channel'].startswith('CH'):
            tr.stats['channel'] = 'HH' + tr.stats['channel'][-1]

    # overwrite mseed
    st.write(mseedfile, format="MSEED")
    
def fix_stream_channels_bb2sp(mseedfile):
    from obspy import read
    
    st = read(mseedfile)
    
    # loop thru traces
    for tr in st:
    
        if tr.stats['channel'].startswith('HH'):
            tr.stats['channel'] = 'EH' + tr.stats['channel'][-1]

    # overwrite mseed
    st.write(mseedfile, format="MSEED")
    
def fix_src_stream_channels(mseedfile):
    from obspy import read
    
    st = read(mseedfile)
    
    # loop thru traces
    st[0].stats['channel'] = 'EHE'
    st[1].stats['channel'] = 'EHN'
    st[2].stats['channel'] = 'EHZ'

    # overwrite mseed
    st.write(mseedfile, format="MSEED")


def fix_stream_network(mseedfile, newnet):
    '''
    newnet = new two letter code, e.g. 'AU'
    '''
    from obspy import read
    
    st = read(mseedfile)
    
    # loop thru traces
    for tr in st:
        tr.stats['network'] = newnet
        #print(tr.stats.network)
    
    # overwrite mseed
    st.write(mseedfile, format="MSEED")
# converts css to mseed
def css2mseed(cssfile, mseedpath):
    '''
    cssfile = wfdisc file in
    mseedfile = path to mseed out file 
    '''
    
    #from obspy.css.core import readCSS
    from obspy import read
    from data_fmt_tools import readGACSS
    from os import path
    
    #cssfile = 'moe_4.4/TOO/TOO_2012202.wfdisc'
    
    try:
        #st = readCSS(cssfile)
        st = read(cssfile)
    except:
        print('Using GA format CSS')
        st = readGACSS(cssfile)
    
    # make outfile
    outfile = '.'.join((st[0].stats['starttime'].strftime('%Y-%m-%dT%H.%M'), 'AU', st[0].stats['station'],'mseed'))
    outpath = path.join(mseedpath, outfile)
    
    # make sure data type ok
    for s in st:
        s.data = 1. * s.data
        
    st.write(outpath, format="MSEED")
    
    return st, outpath    
    
import numpy as np
DTYPE = {
    # Big-endian integers
    b's4': b'>i',
    b's2': b'>h',
    # Little-endian integers
    b'i4': b'<i',
    b'i2': b'<h',
    # ASCII integers
    b'c0': (b'S12', np.int),
    b'c#': (b'S12', np.int),
    # Big-endian floating point
    b't4': b'>f',
    b't8': b'>d',
    # Little-endian floating point
    b'f4': b'<f',
    b'f8': b'<d',
    # ASCII floating point
    b'a0': (b'S15', np.float32),
    b'a#': (b'S15', np.float32),
    b'b0': (b'S24', np.float64),
    b'b#': (b'S24', np.float64),
}

def get_ga_channel(chin, sensor):
    chin = chin.strip().decode()
    if chin == 'ez':
        chin_return = 'HHZ'
    elif chin == 'en':
        chin_return = 'HHN'
    elif chin == 'ee':
        chin_return = 'HHE'
    elif chin == 'sz':
        chin_return = 'BHZ'
    elif chin == 'sn':
        chin_return = 'BHN'
    elif chin == 'se':
        chin_return = 'BHE'
    elif chin == 'gz':
        chin_return = 'HNZ'
    elif chin == 'gn':
        chin_return = 'HNN'
    elif chin == 'ge':
        chin_return = 'HNE'
    else:
        chin_return = chin
        
    if sensor == 'CMG40T' and chin_return.startswith('HH'):
        chin_return.replace('HH', 'EH')
    elif sensor == 'CMG40T' and chin_return.startswith('BH'):
        chin_return.replace('BH', 'SH')
        
    return chin_return

# hacked version of obspy core to read old GA files
def readGACSS(filename, **kwargs):
    import os

    from obspy import Stream, Trace, UTCDateTime
    from obspy.core.compatibility import from_buffer

    """
    Reads a CSS waveform file and returns a Stream object.

    .. warning::
        This function should NOT be called directly, it registers via the
        ObsPy :func:`~obspy.core.stream.read` function, call this instead.

    :type filename: str
    :param filename: CSS file to be read.
    :rtype: :class:`~obspy.core.stream.Stream`
    :returns: Stream with Traces specified by given file.
    """
    # read metafile with info on single traces
    with open(filename, "rb") as fh:
        lines = fh.readlines()
        
    basedir = os.path.dirname(filename)
    traces = []
    # read single traces
    for line in lines:
        dat = line.strip().split()
        
        # format 1
        try:
            offset = int(dat[16])
            npts = int(dat[4])
            dirname = dat[14].strip().decode()
            filename = dat[15].strip().decode()
            dtype = DTYPE[dat[10]]
            fmt1 = True
        
        # format 2
        except:
            offset = int(dat[17])
            npts = int(dat[7])
            dirname = dat[15].strip().decode()
            filename = dat[16].strip().decode()
            dtype = DTYPE[dat[13]]
            fmt1 = False
        
        filename = os.path.join(basedir, dirname, filename)
            
        if isinstance(dtype, tuple):
            read_fmt = np.dtype(dtype[0])
            fmt = dtype[1]
        else:
            read_fmt = np.dtype(dtype)
            fmt = read_fmt
        with open(filename, "rb") as fh:
            fh.seek(offset)
            data = fh.read(read_fmt.itemsize * npts)
            data = from_buffer(data, dtype=read_fmt)
            data = np.require(data, dtype=fmt)
        header = {}
        
        if fmt1 == True:
            header['station'] = dat[2].strip().decode()
            header['channel'] = get_ga_channel(dat[3], dat[8])
            header['starttime'] = UTCDateTime(float(dat[1]))
            header['sampling_rate'] = float(dat[5])
            header['calib'] = float(dat[12])
            header['calper'] = float(dat[13])
        else:
            header['station'] = dat[0].strip().decode()
            header['channel'] = get_ga_channel(dat[1], '')
            header['starttime'] = UTCDateTime(float(dat[2])) # ok?
            header['sampling_rate'] = float(dat[8])
            header['calib'] = float(dat[2])
            header['calper'] = float(dat[10])
        
        tr = Trace(data, header=header)
        traces.append(tr)
    return Stream(traces=traces)


# merges seed files for multichannel data    
def append_seed(mseedfiles, outfile):

    '''
    mseedfiles = tuple of filenames
    '''
    
    from obspy import read
    
    #mseedfiles = ('2012171_105400_0a903_2_1.seed', '2012171_105600_0a903_2_1.seed')
    #outfile = '201210191054.GEES.EHE.seed'
    
    st = read(mseedfiles[0])
    
    # fixes for jump data
    if st[0].stats['network'] == 'XX':
        st[0].stats['network'] = 'AU'
    
    if st[0].stats['channel'] == '001':
        st[0].stats['channel'] = 'EHZ'
    elif st[0].stats['channel'] == '002':
        st[0].stats['channel'] = 'EHN'
    elif st[0].stats['channel'] == '003':
        st[0].stats['channel'] = 'EHE'
    elif st[0].stats['channel'] == '004':
        st[0].stats['channel'] = 'HNZ'
    elif st[0].stats['channel'] == '005':
        st[0].stats['channel'] = 'HNN'
    elif st[0].stats['channel'] == '006':
        st[0].stats['channel'] = 'HNE'
    
    if len(mseedfiles) > 1:
        for i in range(1,len(mseedfiles)):
            st += read(mseedfiles[i])
            
            # fixes for jump data
            if st[i].stats['network'] == 'XX':
                st[i].stats['network'] = 'AU'
    
            if st[i].stats['channel'] == '001':
                st[i].stats['channel'] = 'EHZ'
            elif st[i].stats['channel'] == '002':
                st[i].stats['channel'] = 'EHN'
            elif st[i].stats['channel'] == '003':
                st[i].stats['channel'] = 'EHE'
            elif st[i].stats['channel'] == '004':
                st[i].stats['channel'] = 'HNZ'
            elif st[i].stats['channel'] == '005':
                st[i].stats['channel'] = 'HNN'
            elif st[i].stats['channel'] == '006':
                st[i].stats['channel'] = 'HNE'
    
    # now write to file
    st.write(outfile, format='MSEED')

# function to merge data files from the ausarray deployment
def merge_ausarray_channels(folder):       
    from obspy import read
    from misc_tools import listdir_extension
    from sys import argv
    from os import path
    
    files = listdir_extension(folder, 'HHE')
    
    for f in files:
        print(f)
        try:
            st = read(path.join(folder, f))
            
            stn = read(path.join(folder, f.replace('HHE','HHN')))
            st += stn
            
            stz = read(path.join(folder, f.replace('HHE','HHZ')))
            st += stz
            
            # write to file
            tr = st[0]
            
            mseed = path.join('mseed', \
                               '.'.join((tr.stats.starttime.strftime('%Y-%m-%dT%H.%M'), \
                               tr.stats['network'], tr.stats['station'], 'mseed')))
            #mseed = path.join('mseed', f.replace('HHE','mseed'))
            print('Writing: '+mseed)
            st.write(mseed, format="MSEED")
            
        except:
            print('  Cannot read file')
        
# merge mseed files
# mseedfiles = tuple of files
def merge_seed(mseedfiles):
    from obspy.core import read, Stream
    from misc_tools import doy2ymd
    
    st = Stream()
    for i, wave in enumerate(mseedfiles):
        st += read(wave)
        
        # now merge
        st.merge(method=0, fill_value=0)
        
        # set starttime for file name
        if i == 0:
            starttime = st[0].stats['starttime'].strftime('%Y-%m-%dT%H.%M')
    
    #tmpname = mseedfiles[0].strip('.mseed').split('_')
    #print(tmpname
    #ymdhmd = doy2ymd(tmpname[0][0:4], tmpname[0][4:]) + str('%04d' % float(tmpname[1]))
    outfile = '.'.join((starttime, 'AU', st[0].stats['station'],'mrg','mseed'))
    #outfile = '.'.join((ymdhmd,sta,'mrg.mseed'))
    print('Merged file:', outfile)
    st.write(outfile, format='MSEED')
    
# merge mseed files with a given file prefix
# mseedfiles = tuple of files
def merge_seed_prefix(folder, file_prefix, station):
    from obspy.core import read, Stream
    from misc_tools import listdir_file_prefix
    from os import path
    
    '''
    Usage, e.g.:
         merge_seed_prefix('../mseed_dump/', '2004-03-05T00.', 'MEEK')
    '''
    
    # get file list
    mseedfiles = listdir_file_prefix(folder, file_prefix)
    
    st = Stream()
    i = 0
    for wave in mseedfiles:
        st_tmp = read(path.join(folder, wave))
        if st_tmp[0].stats.station == station:
            st += st_tmp
            
            # now merge
            st.merge(method=0, fill_value=0)
            
            # set starttime for file name
            if i == 0:
                starttime = st[0].stats['starttime'].strftime('%Y-%m-%dT%H.%M')
            
            i += 1
    
    outfile = path.join(folder, '.'.join((starttime, 'AU', st[0].stats['station'],'mrg','mseed')))
    print('Merged file:', outfile)
    st.write(outfile, format='MSEED')

def merge_seed_extract_centaur(folder, file_prefix, station, eventDateTime):
    from obspy.core import read, Stream, UTCDateTime
    from misc_tools import listdir_file_prefix
    from os import path
    
    '''
    folder: folder in which mseed files are located
    file_prefix: text string common to all files we want to merge
    station: station code 
    eventDateTime: event origintime to nearest minute, e.g. datetime.datetime(2018,12,1,13,27)
    '''
    
    # get file list
    mseedfiles = listdir_file_prefix(folder, file_prefix)
    
    st = Stream()
    for wave in mseedfiles:
        st_tmp = read(path.join(folder, wave))
        if st_tmp[0].stats.station == station:
            st += st_tmp
            
            # now merge
            st.merge(method=0, fill_value=0)
    
    # now cut to starttime
    starttime = UTCDateTime(eventDateTime.year, eventDateTime.month, eventDateTime.day, eventDateTime.hour, eventDateTime.minute) - 120
    endtime = starttime + 1200
    
    # now trim
    st_trim = st.trim(starttime, endtime)
    
    outfile = '.'.join((starttime.strftime('%Y-%m-%dT%H.%M'), 'AU', st_trim[0].stats['station'],'mseed'))
    print('Merged file:', outfile)
    st.write(outfile, format='MSEED')

# merge jump seed files to one
def merge_jump_seed(seedfolder, hhmm):
    
    from misc_tools import listdir_extension
    from data_fmt_tools import append_seed, merge_seed
    from os import sep, path, remove
    
    if not isinstance(hhmm, float):
        hhmm = float(hhmm)
    #seedfolder = '/Users/tallen/Documents/Earthquake_Data/GA_Network_Data/GHSS/2012171'
    #hhmm = 0246
    #print(seedfolder, hhmm
        
    if seedfolder.endswith(sep):
        seedfolder = seedfolder[0:-1]
    
    # get yyyydoy
    yyyydoy = seedfolder.split(sep)[-1]
    
    # get station
    stn = seedfolder.split(sep)[-2]
    
    # look for seed files
    seedfiles = listdir_extension(seedfolder, 'seed')
    #print(seedfiles
    
    # first, append files with similar timestamp
    quitAppend = False
    inc = 0
    mrgfiles = []
    while quitAppend == False:
        appendfiles = []
        for file in seedfiles:
            #print('_'.join((yyyydoy, str('%04d' % (hhmm+inc)))), hhmm
            if file.startswith('_'.join((yyyydoy, str('%04d' % (hhmm+inc))))):
            
                appendfiles.append(path.join(seedfolder,file))
        
        if not appendfiles: # quit append
            quitAppend = True
        
        else: # write to seed
            outseed = '_'.join((yyyydoy, str(hhmm+inc), stn.upper())) + '.mseed'
            append_seed(appendfiles, outseed)
            inc += 2
            
            # add for merging
            mrgfiles.append(outseed)
            
    # now merge into one seed file
    print(mrgfiles)
    merge_seed(mrgfiles)
    
    # remove tmp files
    for tmpfile in mrgfiles:
        remove(tmpfile)

# splits a mseed file with multiple stations to a single file per station
# with multiple channels
def split_mseed_stations(mseedfile, out_prefix):
    from obspy import read, Stream
    from numpy import array, unique
    
    '''
    mseedfile = path to mseed file, e.g. 'mseed_dump/2018-04-07 0547 58.ms'
    out_prefix = prefix for output file name, e.g. '2018-04-07T05.48.AU'
    '''
    
    # read streams
    st = read(mseedfile)
    
    # get unique stations
    stas = []
    for tr in st:
        stas.append(tr.stats.station)
    ustas = unique(array(stas))
    
    # now loop through unique stations and dump data
    for sta in ustas:
        new_trs = []
        for tr in st:
            if tr.stats.station == sta:
                new_trs.append(tr)
        
        # make new stream for one station
        new_st = Stream(traces=new_trs)
        
        # set filename
        newfile = '.'.join((out_prefix, sta, 'mseed'))
        
        # write mseed file
        print('    '+newfile)
        new_st.write(newfile, format="MSEED")
        
# trims segment from larger mseed file
def trim_seed(eventDateTuple, mseedFile):
    from obspy import read, UTCDateTime
    from os import path
    
    # format UTC datetime
    dtsplit = [str(x) for x in eventDateTuple] # convert dt to string  
    #print( dtsplit
    utcdt = '-'.join((dtsplit[0], dtsplit[1].zfill(2), dtsplit[2].zfill(2))) \
            + 'T' + ':'.join((dtsplit[3].zfill(2), dtsplit[4].zfill(2), '00.000'))
    
    satrttime = UTCDateTime(utcdt) - 120
    endtime = satrttime + 1500
    
    # read trace
    st = read(mseedFile)
    tr = st[0]
    
    # check end times
    if endtime > tr.stats['endtime']:
        endtime = tr.stats['endtime']
    
    # copy trace
    st_trim = st.copy()
    
    # trim trace
    st_trim = st_trim.trim(starttime=satrttime, endtime=endtime)
    trt = st_trim[0]
    # make output file
    
    outfile = '.'.join((trt.stats.starttime.strftime('%Y-%m-%dT%H.%M'), \
                        trt.stats['network'], trt.stats['station'], 'mseed'))
    
    print('Writing', outfile)
    st.write(outfile, format='MSEED')
    
    return st_trim

# trim big mseed file to something more manageable
def trim_mseed(mseedIn, dateStr, durn=600):
    from obspy import read, UTCDateTime, Stream, Trace
    from misc_tools import ymd2doy
    import datetime as dt
    from sys import argv
    from os import path
    
    starttime = UTCDateTime(dateStr) # e.g. 2017-01-01T04:11:33
    
    # read big file
    st = read(mseedFile)
    
    endtime = peakTime + durn
    
    # now trim
    st_trim = st.trim(starttime, endtime)
    
    # set filename
    mseedOut = '.'.join((st_trim[0].stats.starttime.strftime('%Y-%m-%dT%H.%M'), \
    
                          st_trim[0].stats['network'], st_trim[0].stats['station'], \
                          st[0].stats.channel, 'mseed'))
                          
    # write file
    print('Writing file:', mseedOut)
    outpath = path.join(mseedOut)           
    st_trim.write(outpath, format="MSEED")


# REMOVE DC OFFSET
def remove_dc_offset(data):
    from numpy import round, median
    data = round(data - median(data))

    return data
    
# mseed to ascii
def mseed2asc(mseedfile, ascout):
    from obspy.core import read
    from numpy import array, savetxt
    
    # read mseed
    st = read(mseedfile)
    
    # make data array
    data = []
    chans = ''
    for tr in st:
        data.append(tr.data)
        chans += tr.stats['channel'] + '\t'
    data = array(data)
    chans = chans[0:-1]
    
    # make header
    header  = 'station: ' + tr.stats['station'] + '\n'
    header += 'start time: ' + tr.stats['starttime'].strftime("%Y-%m-%d %H:%M:%S.%f")[0:-4] + ' UTC\n'
    header += 'sampling rate: ' + str(tr.stats['sampling_rate']) + ' Hz\n'
    header += 'units: g\n'
    header += 'npts: ' + str(tr.stats['npts']) + '\n'
    header += chans
    
    savetxt(ascout, data.T, fmt='%.6e', delimiter='\t', header=header)
    
def eqwave2mseed(eqwfile):
    from readwaves import readeqwave 
    from write_data import export_SAC
    from obspy.core import Trace, Stream, UTCDateTime, AttribDict
    from obspy.core.util import calcVincentyInverse
    from os import path, mkdir
    from numpy import array
    
    wavfile = 'waves//2012-07-20_0910_52_MOE4.txt'
    allsta, comps, allrecdate, allsec, allsps, alldata, allnsamp = readeqwave(wavfile)
    
    # set file and path names
    outfile = wavfile.strip('.txt') + '.mseed'
    #outdir = 'mseed'
    #outfile = path.join(outdir,filename)
    
    # try saving to sac directory
    try:
        f = open(outfile,'w')
    except:
        mkdir('mseed')
    
    st = Stream()
    for i, comp in enumerate(comps):
        data = array(alldata[i])
        stats = {'network': 'AU', 'station': allsta[i], \
                 'channel': comp, 'npts': allnsamp[i], 'sampling_rate': allsps[i]}
        
        ymd = allrecdate[i][0:8]
        hhmm = allrecdate[i][8:]
        stats['starttime'] = UTCDateTime(int(allrecdate[i][0:4]), int(allrecdate[i][4:6]), int(allrecdate[i][6:8]), \
                                           int(allrecdate[i][8:10]), int(allrecdate[i][10:12]), float(allsec[i]))
        print(stats['starttime'])
        tr = Trace(data=data, header=stats)
            
        # fill stream
        st.append(tr)
        
    # calculate distance and azimuth
    #stats['sac'] = sac
    
    st.write(outfile, format='MSEED')


def get_iris_data(dateTuple, sta, net, durn=600, client='IRIS'):     
    from obspy.core.utcdatetime import UTCDateTime
    from obspy.clients.fdsn.client import Client
    from os import path, makedirs
    
    '''
    Code to extract IRIS data, one station at a time.  Exports mseed file to 
    working directory
    
    datetime tuple fmt = (Y,m,d,H,M)
    sta = station
    durn = duration in secs
    '''
    
    #try:
    # get inputs
    #datetime = argv[1] # fmt = (Y,m,d,H,M)
    #sta = argv[2].upper() # station code
    
    # format UTC datetime
    dtsplit = [str(x) for x in dateTuple] # convert dt to string  
    #print( dtsplit
    utcdt = '-'.join((dtsplit[0], dtsplit[1].zfill(2), dtsplit[2].zfill(2))) \
            + 'T' + ':'.join((dtsplit[3].zfill(2), dtsplit[4].zfill(2), '00.000'))
    
    data_client = Client(client)
    t1 = UTCDateTime(utcdt) - 120
    t2 = t1 + durn
    t3 = t1 + 3
    
    bulk = [(net.upper(), sta.upper(), "*", "*", t1, t2),
            ("AU", "AFI", "1?", "BHE", t1, t3)]
    
    try:        
        st = data_client.get_waveforms_bulk(bulk)
        
        # save out to file
        tr = st[0]
        
        trname = path.join('iris_dump', \
                           '.'.join((tr.stats.starttime.strftime('%Y-%m-%dT%H.%M'), \
                           tr.stats['network'], tr.stats['station'], 'mseed')))
        
        # check if waves folder exists
        if not path.isdir('iris_dump'):
            makedirs('iris_dump')
            
        print('Writing file:', trname)
        st.write(trname, format="MSEED")
    except:
        print('Data not available:', sta.upper())
        # dummy data returned
        st = 0
        trname='null'
    '''
    except:
        print('\nUsage: \n    python get_iris_data.py <datetime tuple> <station code>\n'
        
        print('    e.g.: python get_iris_data.py (2010,4,20,0,16) kmbl\n'
    '''
    
    return st, trname

def get_iris_inventory(net, stalist):
    '''
    net = network code (e.g. "II")
    stalist = comma separated string (e.g., "TAU,WRAB")
    '''
    from obspy.clients.fdsn.client import Client
    
    client = Client("IRIS")
    
    inventory = client.get_stations(network=net, station=stalist, level="response")
    
    xmlname = net.lower()+'-inventory.xml'
    inventory.write(xmlname,'stationxml')

def get_auspass_data(dateTupple, durn=600, network='S1', station='*', channel='*'):
    '''
    datetime tuple fmt = (Y,m,d,H,M)
    '''
    from obspy import UTCDateTime, Stream
    from obspy.clients.fdsn import Client as FDSN_Client
    from numpy import unique
    from os import path, makedirs
    
    auspass = FDSN_Client('http://auspass.edu.au:8080',user_agent='trevor.allen@ga.gov.au') 
    
    '''
    Networks:
        2P = SWAN
        1K = ALFREX
        M8 = semi perm
        WG = WA Array
    '''
    location = "*" #n.b. AusPass data typically has a univeral blank (e.g. '') value for ALL station locations. "*" works just as well. 
    time0 = UTCDateTime(dateTupple[0],dateTupple[1],dateTupple[2],dateTupple[3],dateTupple[4]) - 120
    time1 = time0 + durn 
    
    #try:
    st = auspass.get_waveforms(network=network,station=station,location=location,channel=channel,starttime=time0,endtime=time1)
    
    if len(st) > 0:
        # get array of stations
        stalist = []
        for tr in st:
            stalist.append(tr.stats.station)
            
        stalist = unique(stalist)
        
        # check if waves folder exists
        if not path.isdir('auspass_dump'):
            makedirs('auspass_dump')
        
        # now loop thru unique stalist and output streams
        for us in stalist:
            new_st = Stream()
            for tr in st:
                if tr.stats.station == us:
                    new_st += tr
                    
            # save out to file
            tr = new_st[0]
            
            trname = path.join('auspass_dump', \
                               '.'.join((tr.stats.starttime.strftime('%Y-%m-%dT%H.%M'), \
                               tr.stats['network'], tr.stats['station'], 'mseed')))
            
            print('Writing file:', trname)
            new_st.write(trname, format="MSEED")
            
    #except:
    #    print('No AusPass data available...')
        
def get_swan_data(dateTupple, durn=600, network='2P', station='*', channel='*'):
    '''
    datetime tuple fmt = (Y,m,d,H,M)
    '''
    from obspy import UTCDateTime, Stream
    from obspy.clients.fdsn import Client as FDSN_Client
    from numpy import unique
    from os import path, makedirs
    
    if network == '2P':
        auspass = FDSN_Client('http://auspass.edu.au:8080',user_agent='trevor.allen@ga.gov.au',user='2p',password='kiwi30') 
    elif network == 'WG':
        auspass = FDSN_Client('http://auspass.edu.au:8080',user_agent='trevor.allen@ga.gov.au',user='WG',password='bunglebungle') 
    
    '''
    Networks:
        2P = SWAN
        1K = ALFREX
    '''
    location = "*" #n.b. AusPass data typically has a univeral blank (e.g. '') value for ALL station locations. "*" works just as well. 
    time0 = UTCDateTime(dateTupple[0],dateTupple[1],dateTupple[2],dateTupple[3],dateTupple[4]) - 120
    time1 = time0 + durn 
    
    try:
        st = auspass.get_waveforms(network=network,station=station,location=location,channel=channel,starttime=time0,endtime=time1)
        
        if len(st) > 0:
            # get array of stations
            stalist = []
            for tr in st:
                stalist.append(tr.stats.station)
                
            stalist = unique(stalist)
            
            # check if waves folder exists
            if not path.isdir('auspass_dump'):
                makedirs('auspass_dump')
            
            # now loop thru unique stalist and output streams
            for us in stalist:
                new_st = Stream()
                for tr in st:
                    if tr.stats.station == us:
                        new_st += tr
                        
                # save out to file
                tr = new_st[0]
                
                trname = path.join('auspass_dump', \
                                   '.'.join((tr.stats.starttime.strftime('%Y-%m-%dT%H.%M'), \
                                   tr.stats['network'], tr.stats['station'], 'mseed')))
                
                print('Writing file:', trname)
                new_st.write(trname, format="MSEED")
                
    except:
        print('No AusPass data available...')
    
def get_arclink_data(datetupple, sta, net):
    from obspy import UTCDateTime
    from obspy.clients.arclink.client import Client
    from os import path
    
    client = Client(user='test@obspy.org')
    
    # format UTC datetime
    dtsplit = [str(x) for x in datetupple] # convert dt to string  
    #print( dtsplit
    utcdt = '-'.join((dtsplit[0], dtsplit[1].zfill(2), dtsplit[2].zfill(2))) \
            + 'T' + ':'.join((dtsplit[3].zfill(2), dtsplit[4].zfill(2), '00.000'))
    
    t1 = UTCDateTime(utcdt) - 120
    t2 = t1 + 1500
    
    # save out to file
    #try:
    st = client.get_waveforms(net, sta, "", "*", t1, t2)
    tr = st[0]
    
    trname = path.join('iris_dump', \
                       '.'.join((tr.stats.starttime.strftime('%Y-%m-%dT%H.%M'), \
                       tr.stats['network'], tr.stats['station'], 'mseed')))
    
    # check if waves folder exists
    if not path.isdir('iris_dump'):
        makedirs('iris_dump')
        
    print('Writing file:', trname)
    st.write(trname, format="MSEED")
    '''
    except:
        print('Data not available:', sta.upper())
        # dummy data returned
        st = 0
        trname='null'
    '''
    return st, trname

def get_nat_cwb_data(Y,m,d,H,M,td_start, td_end):
    '''
    origintime = tuple fmt (Y,m,d,H,M)
    td_start, td_end: time deltas in seconds
    '''
    import datetime
    from os import path, makedirs
    from obspy.core import utcdatetime
    #from obspy.core.event import Event, Magnitude, Origin
    from obspy.clients.neic.client import Client
    
    ##########################################################################
    # set constants and funcs
    ##########################################################################
    
    # initialize the cwb port
    client=Client(host='10.7.161.60',port=2061,debug=False, timeout=60)
    
    ##########################################################################
    # set event time
    ##########################################################################
    dt = datetime.datetime(Y,m,d,H,M)
    
    print(dt)
    
    # convert datetime object to UTCdatetime
    dt = utcdatetime.UTCDateTime(dt)
    
    ''' the time window to request the data will be 20 minutes, check maximum travel time and increase this value accordingly '''
    #end_time=start_time+960 # 16 minutes
    start_time = dt + datetime.timedelta(seconds=td_start)
    end_time   = dt + datetime.timedelta(seconds=td_end) # 16 minutes
    #end_time   = dt + datetime.timedelta(seconds=600) # 5 minutes
    
    
    ''' get all waveform data available, use wildcards to reduce the data volume and speedup the process,
    unfortunately we need to request few times for every number of characters that forms the station name '''
    # kluge to fix non-retrieval of data  - loop through alphabet integers
    for ch in range(ord('A'), ord('Z')+1):
        print('Stations beginning with ', chr(ch))
    #    st_3 = client.get_waveforms("AU", chr(ch)+"??", "", "[BSEH]?[ENZ]", start_time,end_time)
    #    st_4 = client.get_waveforms("AU", chr(ch)+"???", "", "[BSEH]?[ENZ]", start_time,end_time)
        st_3 = client.get_waveforms("AU", chr(ch)+"??", "", "[BSEH][HN][ENZ]", start_time,end_time)
        st_4 = client.get_waveforms("AU", chr(ch)+"???", "", "[BSEH][HN][ENZ]", start_time,end_time)
        if ch == ord('A'):
            if len(st_4) > 0:
                st=st_3+st_4
            else:
                st=st_3
        
        else:
            if len(st_4) > 0:
                st+=st_3+st_4
            else:
                st+=st_3
    
    # Cleanup duplicate traces returned by server
     #       st.merge(-1) #-1 method only merges overlapping or adjacent traces with same i            
    # Now sort the streams by station and channel
    st.sort()
    # Cleanup duplicate traces returned by server
    for tr in st:
        if tr.stats['sampling_rate'] < 20:
            st.remove(tr)
        
    st.merge(0, fill_value='interpolate') #1 method only merges overlapping or adjacent traces with same id
    # Now sort the streams by station and channel
    st.sort()
    
    # check if waves folder exists
    if not path.isdir('cwb_dump'):
        makedirs('cwb_dump')
        
    # set mseed filename
    msfile = path.join('cwb_dump', dt.strftime('%Y%m%d%H%M')+'.mseed')
    
    # now write streams for each event to mseed
    st.write(msfile, format="MSEED") 
    print(st)
    
def get_sta_cwb_data(Y,m,d,H,M, td_start, td_end, sta):
    '''
    origintime = tuple fmt (Y,m,d,H,M)
    td_start, td_end: time deltas in seconds
    '''
    import datetime
    from os import path, makedirs
    from obspy.core import utcdatetime
    #from obspy.core.event import Event, Magnitude, Origin
    from obspy.clients.neic.client import Client
    
    ##########################################################################
    # set constants and funcs
    ##########################################################################
    
    # initialize the cwb port
    client=Client(host='10.7.161.60',port=2061,debug=False, timeout=60)
    
    ##########################################################################
    # set event time
    ##########################################################################
    dt = datetime.datetime(Y,m,d,H,M)
    
    # convert datetime object to UTCdatetime
    dt = utcdatetime.UTCDateTime(dt)
    
    ''' the time window to request the data will be 20 minutes, check maximum travel time and increase this value accordingly '''
    #end_time=start_time+960 # 16 minutes
    start_time = dt + datetime.timedelta(seconds=td_start)
    end_time   = dt + datetime.timedelta(seconds=td_end) # 16 minutes
    #end_time   = dt + datetime.timedelta(seconds=600) # 5 minutes
    
    
    ''' get all waveform data available, use wildcards to reduce the data volume and speedup the process,
    unfortunately we need to request few times for every number of characters that forms the station name '''
    st = client.get_waveforms("AU", sta, "", "[BSEH][HN][ENZ]", start_time,end_time)
    
    # Cleanup duplicate traces returned by server
     #       st.merge(-1) #-1 method only merges overlapping or adjacent traces with same i            
    # Now sort the streams by station and channel
    st.sort()
    
    # Cleanup duplicate traces returned by server
    for tr in st:
        if tr.stats['sampling_rate'] < 20:
            st.remove(tr)
        
    # remove low-sample recs
    st = remove_low_sample_data(st)
    
    st.merge(0, fill_value='interpolate') #1 method only merges overlapping or adjacent traces with same id
    
    # Now sort the streams by station and channel
    st.sort()
    
    # check if waves folder exists
    if not path.isdir('cwb_dump'):
        makedirs('cwb_dump')
        
    # set mseed filename
    msfile = path.join('cwb_dump', start_time.strftime('%Y-%m-%dT%H.%M')+ '.AU.' + sta + '.mseed')
    
    # now write streams for each event to mseed
    if len(st) > 0:
        st.write(msfile, format="MSEED") 
        print(st)
        
    return st, msfile
        
# parse station dataless seed
def get_stn_dataless_seed(network):
    from obspy.io.xseed import Parser
    from os import path, getcwd
        
    # set dataless path
    if getcwd().startswith('/nas'):
        try:
            dataless = Parser(path.join('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks', network, network+'.dataless'))
        except:
            dataless = Parser('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/AU/AU_dataless.286690.dataless')
        
    return dataless

# parse station dataless seed
def get_iris_stn_dataless_seed(network):
    from obspy.io.xseed import Parser
    from os import path, getcwd
        
    # set dataless path
    if getcwd().startswith('/nas'):
        dataless = Parser(path.join('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks', network, network+'.IRIS.dataless'))
    elif getcwd().startswith('C:'):
        print(path.join('C:\\Users\\u56903\\OneDrive - Geoscience Australia\\Networks', network, network+'.IRIS.dataless'))
        dataless = Parser(path.join('C:\\Users\\u56903\\OneDrive - Geoscience Australia\\Networks', network, network+'.IRIS.dataless'))
    else:    
        dataless = Parser(path.join('..\\..\\..\\Networks', network, network+'.IRIS.dataless'))
        
    return dataless

# gets station distance from dataless seed volume
def get_station_distance(st, dataless, eqlo, eqla):
    '''
    st = station miniseed stream
    dataless = dataless seed volume    
    '''
    
    from mapping_tools import distance
    from numpy import array
    
    # set arrays to return
    channel = []
    repi = []
    azim = []
    stations = []
    stalocs = []
    
    for tr in st:
        seedid=tr.get_id()
        channel.append(tr.stats.channel)
        stations.append(tr.stats.station)
        
        start_time = tr.stats.starttime
        
        #paz = dataless.get_paz(seedid,start_time)
        staloc = dataless.get_coordinates(seedid,start_time)
        stalocs.append(staloc)
        
        repi.append(distance(eqla, eqlo, staloc['latitude'], staloc['longitude'])[0])
        azim.append(distance(eqla, eqlo, staloc['latitude'], staloc['longitude'])[1])
               
    return stalocs, array(stations), array(channel), array(repi), array(azim)

def get_station_distance_stadat(sta, eqlo, eqla): 
    from data_fmt_tools import return_sta_data
    from mapping_tools import distance
    
    sta_data = return_sta_data(sta)
    
    distkm = distance(eqla, eqlo, sta_data['stla'], sta_data['stlo'])[0]
    azim = distance(eqla, eqlo, sta_data['stla'], sta_data['stlo'])[1]
    
    return distkm, azim
   
def return_all_au_station_data():
    from datetime import datetime
    from os import getcwd
    
    if getcwd().startswith('C:') or getcwd().startswith('Z:'):
        #au_station_file = '/nas/users/u56903/unix/Code/my_codes/au_station_data.dat'
        au_station_file = 'C:\\Code\\my_codes\\au_station_data.dat' # windows
    else:
        au_station_file = '/Users/trev/Documents/Code/my_codes/au_station_data.dat'
    
    lines = open(au_station_file).readlines()[1:]
    
    sta_dict = []
    
    for line in lines:
        dat = line.strip().split('\t')
        #print(line)
        
        if int(dat[5]) < 1:
            dat[5] = 1
        if int(dat[7]) < 1:
            dat[7] = 1
            
        tmp = {'sta':dat[0], 'stlo':float(dat[1]), 'stla':float(dat[2]), 
               'startdate':datetime(int(dat[4]),int(dat[5]),1), 
               'enddate':datetime(int(dat[6]),int(dat[7]),1)}
        
        # append to sta_dict
        sta_dict.append(tmp)
        
    return sta_dict

# get single station data
def return_sta_data(sta):
    
    #sta = sta.encode('ascii', 'ignore')
    sta_dict = return_all_au_station_data()
    
    for sd in sta_dict:
        if sd['sta'] == sta:
            sta_data = sd
            
    return sta_data

# parses txt files downloaded from e.g.: http://ds.iris.edu/gmap/#network=AU&planet=earth
def parse_iris_stationlist(stationlist):
    '''
    look at response.iris_gmap2stationlist for export formats
    '''
    
    from obspy import UTCDateTime
    from numpy import nan
    
    lines = open(stationlist).readlines() #[:]
    
    if lines[0].startswith('#Network'):
       startidx = 3
    else:
        startidx = 0
    
    staDict = []
    
    for line in lines[startidx:]:
        dat = line.strip().split('|')
        #print(dat)
        try:
            start = UTCDateTime(dat[6]).datetime
            stop  = UTCDateTime(dat[7]).datetime
            tmp = {'sta': dat[1], 'lat': float(dat[2]), 'lon': float(dat[3]), \
                   'elev': float(dat[4]), 'place': dat[5], 'net': dat[0], \
                   'starttime':start, 'stoptime':stop, 'sensitivity':nan}
        except:
            start = UTCDateTime(dat[-2]).datetime
            stop  = UTCDateTime(dat[-1]).datetime
            tmp = {'sta': dat[1], 'lat': float(dat[4]), 'lon': float(dat[5]), \
                   'elev': float(dat[6]), 'place': dat[1], 'net': dat[0], \
                   'sensitivity': float(dat[11]), 'starttime':start, 'stoptime':stop}
        
        staDict.append(tmp)
        
    return staDict

def iris_gmap2stationlist(iris_gmap):
    from data_fmt_tools import parse_iris_stationlist
    from misc_tools import dictlist2array
    from numpy import unique
    
    staDict = parse_iris_stationlist(iris_gmap)
    
    # remove channels
    stations = dictlist2array(staDict, 'sta')
    stations = unique(stations)
    
    staDict2 = []
    
    for s in stations:
        for sta in staDict:
            if sta['sta'] == s:
                s2 = sta
                
        staDict2.append(s2)
    
    # now write txt
    txt = ''
    txt2 = ''
    digsen = str('%0.1f' % (sta['sensitivity'] / 754.3)) # assume trillium 120
    for sta in staDict2:
    	txt += '\t'.join((sta['sta'],'B',sta['starttime'].strftime('%Y%m%d'),sta['stoptime'].strftime('%Y%m%d'),str(sta['lon']),str(sta['lat']),sta['net'],'-12345','-12345','754.3',digsen,'1','HHE','trillium-compact-120.paz')) + '\n'
    	txt += '\t'.join((sta['sta'],'B',sta['starttime'].strftime('%Y%m%d'),sta['stoptime'].strftime('%Y%m%d'),str(sta['lon']),str(sta['lat']),sta['net'],'-12345','-12345','754.3',digsen,'1','HHN','trillium-compact-120.paz')) + '\n'
    	txt += '\t'.join((sta['sta'],'B',sta['starttime'].strftime('%Y%m%d'),sta['stoptime'].strftime('%Y%m%d'),str(sta['lon']),str(sta['lat']),sta['net'],'-12345','-12345','754.3',digsen,'1','HHZ','trillium-compact-120.paz')) + '\n'
    	
    	txt2 += '\t'.join((sta['sta'],str(sta['lon']),str(sta['lat']),'1',str(sta['starttime'].year), str(sta['starttime'].month), \
    	                   str(sta['stoptime'].year), str(sta['stoptime'].month))) + '\n'
    # now write
    f = open(sta['net']+'_stationlist.txt', 'w')
    f.write(txt)
    f.close()
    
    f = open(sta['net']+'_station_data.dat', 'w')
    f.write(txt2)
    f.close()


# this has not been tested!
def eqwave2stationlist(eqwave_txtfile):
    from readwaves import readeqwave
    from sys import argv
    
    wavfile = argv[1]
    sta, comps, recdate, sec, sps, data, nsamp, cntpvolt, sen, gain = readeqwave(wavfile)
    
    lines = ''
    for i in range(0, len(sta)):
        line = '\t'.join((sta[i], comps[i][0], recdate[i][0:8], '25990101', 'LON', 'LAT', 'MEL', '-12345', '-12345', \
                          sen[i], cntpvolt[i], gain[i], comps[i], 'NULL'))
        lines += line + '\n'
    
    # make filename
    slfile = wavfile.split('.')[0].split('/')[1]+'.stationlist'
    slfile = slfile.replace(' ', '_')
    f = open(slfile, 'wb')
    f.write(lines)
    f.close()




def remove_acceleration_data(st):
    # first split the stream and remerge
    st = st.split().split().merge(method=1, fill_value=0)
    
    for tr in st:
        if tr.stats.channel[1] == 'N':
            print('    Removing '+tr.get_id())
            st = st.remove(tr)
            
    st.merge()
    return st

def remove_low_sample_data(st):
    from numpy import array, unique, zeros_like, where
    
    for tr in st:
        if tr.stats.sampling_rate < 1.0:
            st = st.remove(tr)
            
    # first split the stream and remerge
    st = st.split().split().merge(method=1, fill_value=0)
    
    rmlsr = False
    for tr in st:
        if tr.stats.sampling_rate >= 80:
            rmlsr = True
        
            
        # remove junk stations
        try:
            if tr.stats.channel.encode('ascii','ignore').startswith('L') \
                or tr.stats.channel.encode('ascii','ignore').startswith('V') \
                or tr.stats.channel.encode('ascii','ignore').startswith('U') \
                or tr.stats.sampling_rate < 1.0:
                st = st.remove(tr)
                
                #or tr.stats.channel.encode('ascii','ignore').startswith('C') \
        except:
            if tr.stats.channel.startswith('L') \
                or tr.stats.channel.startswith('V') \
                or tr.stats.channel.startswith('U') \
                or tr.stats.sampling_rate < 1.0:
                st = st.remove(tr)
                
                #or tr.stats.channel.startswith('C') \
                
    # now strip low sample data from stream
    if rmlsr == True:
        for tr in st:
            if tr.stats.sampling_rate < 80:
                st = st.remove(tr)
                
    # now look for dupicate station codes with different sampling rates
    channels = []
    sampling_rates = []
    for tr in st:
        try:
            channels.append(tr.stats.channel.encode('ascii','ignore'))
        except:
            channels.append(tr.stats.channel)
        sampling_rates.append(tr.stats.sampling_rate)
    
    channels = array(channels)
    sampling_rates = array(sampling_rates)
    delidx = zeros_like(sampling_rates)
    
    unique_channels = unique(array(channels))
    
    # loop thru unique channels and remove lsr duplicates
    for uc in unique_channels:
        
        idx = where(channels == uc)[0]
        maxSR = max(sampling_rates[idx])
        #print(uc, idx, maxSR
        
        for i in idx:
            if not sampling_rates[i] == maxSR:
                delidx[i] = 1.
    
    #print(delidx            
    # now purge channels
    for tr, di in zip(st, delidx):
        if di == 1.:
            st = st.remove(tr)
            
    # finally, if duplicates
    st.merge()
                       
    return st
 
def return_trace_datetime_array(tr):
    from datetime import timedelta
    
    # get times array
    times = tr.times()
    
    dt_times = []
    for time in times:
        dt_times.append(tr.stats.starttime.datetime + timedelta(seconds=time))
        
    return dt_times

# reformates date strings from obspy plotting
def reformat_trace_mpl_labels(ax):
    import matplotlib as mpl 
    
    ticks = ax.get_xticks()

    ticklabels = [mpl.dates.num2date(t).strftime('%H:%M:%S') for t in ticks]
    ax.set_xticklabels(ticklabels)
    
    return ticklabels
    
    

# parse gmprocess flatfile
def parse_gmprocess_flatfile(flatfile):
    from numpy import array
    
    lines = open(flatfile).readlines()
    
    # get header
    header = lines[0].split(',')
    
    # get periods
    per = []
    for key in header:
        if key.startswith('SA('):
            per.append(float(key.strip('SA(').strip(')')))
    
    per = array(per)
    
    recs = []
    # now loop thru recs
    for line in lines[1:]:
        rec = dict.fromkeys(header, 0)
        
        keys = rec.keys()
        
        dat = line.split(',')
        sa = []
        for i, key in enumerate(keys):
            if key == 'PGA' or key == 'PGV':
                rec[key] = float(dat[i])
            else:
                rec[key] = dat[i]
            
            if key.startswith('SA('):
                sa.append(float(dat[i]))
        
        sa = array(sa)
        rec['per'] = per
        rec['sa'] = sa
        
        # make list
        recs.append(rec)    
        
    return recs
