# converts css to mseed
def css2mseed(cssfile):
    from obspy.css.core import readCSS
    from data_fmt_tools import readGACSS
    
    #cssfile = 'moe_4.4/TOO/TOO_2012202.wfdisc'
    
    try:
        st = readCSS(cssfile)
    except:
        st = readGACSS(cssfile)
    
    # make outfile
    outfile = '.'.join((st[0].stats['starttime'].strftime('%Y%m%d%H%M'), st[0].stats['station'],'mseed'))
    
    # make sure data type ok
    for s in st:
        s.data = 1. * s.data
        
    st.write(outfile, format="MSEED")
    
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

def get_ga_channel(chin):
    if chin == 'ez':
        return 'EHZ'
    elif chin == 'en':
        return 'EHN'
    elif chin == 'ee':
        return 'EHE'
    elif chin == 'sz':
        return 'SHZ'
    elif chin == 'sn':
        return 'SHN'
    elif chin == 'se':
        return 'SHE'
    elif chin == 'gz':
        return 'HNZ'
    elif chin == 'gn':
        return 'HNN'
    elif chin == 'ge':
        return 'HNE'

# hacked version of obspy core to read old GA files
def readGACSS(filename, **kwargs):
    import os

    from obspy import Stream, Trace, UTCDateTime
    from obspy.core.compatibility import frombuffer

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
        print dat
        npts = int(dat[4])
        dirname = dat[14].strip().decode()
        filename = dat[15].strip().decode()
        filename = os.path.join(basedir, dirname, filename)
        offset = int(dat[16])
        dtype = DTYPE[dat[10]]
        if isinstance(dtype, tuple):
            read_fmt = np.dtype(dtype[0])
            fmt = dtype[1]
        else:
            read_fmt = np.dtype(dtype)
            fmt = read_fmt
        with open(filename, "rb") as fh:
            fh.seek(offset)
            data = fh.read(read_fmt.itemsize * npts)
            data = frombuffer(data, dtype=read_fmt)
            data = np.require(data, dtype=fmt)
        header = {}
        header['station'] = dat[2].strip().decode()
        header['channel'] = get_ga_channel(dat[3]).strip().decode()
        header['starttime'] = UTCDateTime(float(dat[1]))
        header['sampling_rate'] = float(dat[5])
        header['calib'] = float(dat[12])
        header['calper'] = float(dat[13])
        tr = Trace(data, header=header)
        traces.append(tr)
    return Stream(traces=traces)


# merges seed files for multichannel data    
def append_seed(mseedfiles, outfile):

    '''
    mseedfiles = tuple of filenames
    '''
    
    from obspy.core import read
    
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
    
# merge mseed files
# mseedfiles = tuple of files
def merge_seed(mseedfiles):
    from obspy.core import read, Stream
    from misc_tools import doy2ymd
    
    st = Stream()
    for wave in mseedfiles:
        st += read(wave)
        
        # now merge
        st.merge(method=0, fill_value=0)
    
    sta = st[0].stats['station']
    tmpname = mseedfiles[0].strip('.mseed').split('_')
    print tmpname
    ymdhmd = doy2ymd(tmpname[0][0:4], tmpname[0][4:]) + str('%04d' % float(tmpname[1]))
    outfile = '.'.join((ymdhmd,sta,'mrg.mseed'))
    print outfile
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
    #print seedfolder, hhmm
        
    if seedfolder.endswith(sep):
        seedfolder = seedfolder[0:-1]
    
    # get yyyydoy
    yyyydoy = seedfolder.split(sep)[-1]
    
    # get station
    stn = seedfolder.split(sep)[-2]
    
    # look for seed files
    seedfiles = listdir_extension(seedfolder, 'seed')
    #print seedfiles
    
    # first, append files with similar timestamp
    quitAppend = False
    inc = 0
    mrgfiles = []
    while quitAppend == False:
        appendfiles = []
        for file in seedfiles:
            #print '_'.join((yyyydoy, str('%04d' % (hhmm+inc)))), hhmm
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
    print mrgfiles
    merge_seed(mrgfiles)
    
    # remove tmp files
    for tmpfile in mrgfiles:
        remove(tmpfile)
    


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
        print stats['starttime']
        tr = Trace(data=data, header=stats)
            
        # fill stream
        st.append(tr)
        
    # calculate distance and azimuth
    #stats['sac'] = sac
    
    st.write(outfile, format='MSEED')


def get_iris_data(dateTuple, sta, net):     
    from obspy.core.utcdatetime import UTCDateTime
    from obspy.clients.fdsn.client import Client
    
    '''
    Code to extract IRIS data, one station at a time.  Exports mseed file to 
    working directory
    
    datetime tuple fmt = (Y,m,d,H,M)
    sta = station
    '''
    
    #try:
    # get inputs
    #datetime = argv[1] # fmt = (Y,m,d,H,M)
    #sta = argv[2].upper() # station code
    
    # format UTC datetime
    dtsplit = [str(x) for x in dateTuple] # convert dt to string  
    #print  dtsplit
    utcdt = '-'.join((dtsplit[0], dtsplit[1].zfill(2), dtsplit[2].zfill(2))) \
            + 'T' + ':'.join((dtsplit[3].zfill(2), dtsplit[4].zfill(2), '00.000'))
    
    client = Client("IRIS")
    t1 = UTCDateTime(utcdt)
    t2 = t1 + 900
    t3 = t1 + 3
    
    bulk = [(net.upper(), sta.upper(), "*", "*", t1, t2),
            ("AU", "AFI", "1?", "BHE", t1, t3)]
            
    st = client.get_waveforms_bulk(bulk)
    
    # save out to file
    tr = st[0]
    #for tr in st:
    trname = '.'.join((tr.stats.starttime.strftime('%Y%m%d%H%M'), \
                       tr.stats['station'], \
                       tr.stats['network'],'mseed'))
    
    print 'Writing file:', trname                   
    st.write(trname, format="MSEED")
    '''
    except:
        print '\nUsage: \n    python get_iris_data.py <datetime tuple> <station code>\n'
        
        print '    e.g.: python get_iris_data.py (2010,4,20,0,16) kmbl\n'
    '''

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
    
    print dt
    
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
        print 'Stations beginning with ', chr(ch)
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
    print st
    
def get_sta_cwb_data(Y,m,d,H,M,td_start, td_end, sta):
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
    
    print dt
    
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
        
    st.merge(0, fill_value='interpolate') #1 method only merges overlapping or adjacent traces with same id
    
    # Now sort the streams by station and channel
    st.sort()
    
    # check if waves folder exists
    if not path.isdir('cwb_dump'):
        makedirs('cwb_dump')
        
    # set mseed filename
    msfile = path.join('cwb_dump', dt.strftime('%Y-%m-%dT%H.%M')+ '.AU.' + sta + '.mseed')
    
    # now write streams for each event to mseed
    if len(st) > 0:
        st.write(msfile, format="MSEED") 
        print st
        