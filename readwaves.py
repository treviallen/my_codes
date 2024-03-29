# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 11:15:56 2012

@author: u56903

Code to read various ASCII wavefile formats, including:
    1) SRC SeisMac
    2) ES&S eqWave
    3) NMX ASCII
    4) TSPAIR (from: obspy.core.ascii.writeTSPAIR)

"""
import numpy as np

def check_file_fmt(wavfile):
    # check to see if mseed first
    from obspy.core import read

    print('\nObsPy Installed')

    # try testing if SAC or miniSEED
    try:
        st = read(wavfile)
        sacseed = True
        fmt = 'mseed'

    except:
        sacseed = False
        print(wavfile)
        f = open(wavfile, 'r')
        firstline = f.readline()
        if firstline.find('Nanometrics') == 0:
            fmt = 'nmx'
        elif firstline.find('ARRIVALS') == 0:
            fmt = 'eqw'
        elif firstline.find('Event:') == 0:
            fmt = 'bkn'
        elif firstline.find('TIMESERIES') == 0:
            fmt = 'tspair'
        elif firstline.find('FileType=TimeSeries') == 0:
            fmt = 'sm'

    return fmt
    
# return data
def return_data(wavfile):  
    fmt = check_file_fmt(wavfile)
    datadict = {}
    if fmt == 'eqw':
        # read the eqWave text file
        #allsta, comps, allrecdate, allsec, allsps, alldata, allnsamp = readeqwave(wavfile)
        allsta, comps, allrecdate, allsec, allsps, alldata, allnsamp, cntpvolt, sen, gain, datadict = readeqwave(wavfile)
    elif fmt == 'nmx':
        allsta, comps, allrecdate, allsec, allsps, alldata, allnsamp = readnmx(wavfile)
    elif fmt == 'tspair':
        allsta, comps, allrecdate, allsec, allsps, alldata, allnsamp = readtspair(wavfile)
    elif fmt == 'sm':
        allsta, comps, allrecdate, allsec, allsps, alldata, allnsamp, datadict = readseismac(wavfile)
    elif fmt == 'bkn':
        allsta, comps, allrecdate, allsec, allsps, alldata, allnsamp, cntpvolt, sen, gain, instrument, secoffset = readbkn(wavfile)
    elif fmt == 'mseed':
        from obspy.core import read
        st = read(wavfile)
        allsta, comps, allrecdate, allsec, allsps, alldata, allnsamp, allnetid = readseed(st)
    else:
        '\nFile format not recognised!'
        
    #print('return_data', alldata[:, 0]
    return allsta, comps, allrecdate, allsec, allsps, alldata, allnsamp, fmt, datadict

def readeqwave(wavfile):
    from datetime import datetime
    
    # READ HEADER INFO
    cntpvolt = 'null' 
    sen = 'null'
    gain = 'null'
    lats = -9999
    
    #print('\nReading header info...'
    header = open(wavfile).readlines()

    readdat = 0

    i = -1
    for k in range(0,len(header)):
        line = header[k]
        line = line.strip('\n')
        line = line.strip('\t')
        line = line.replace('\t\t','\t')
        
        # get station name
        ind = line.find('#sitename')
        if ind >= 0:
            sta = line.split('\t')
            sta = sta[0:-1]

        # get components
        ind = line.find('#component')
        if ind >= 0:
            comps = line.split('\t')
            comps = comps[0:-1]

        # get start time
        ind = line.find('#year month day')
        if ind >= 0:
            ymd = line.split('\t')
            ymd = ymd[0:-1]

        ind = line.find('#hour minute')
        if ind >= 0:
            hhmm = line.split('\t')
            hhmm = hhmm[0:-1]

        ind = line.find('#second')
        if ind >= 0:
            sec = line.split('\t')
            if len(sec) >= 2:
                sec = sec[0:-1]
                if sec[-1] == '': # This seems to be something odd with eqWave3!
                    sec = sec[0:-1]
                intsec = []
                for j in range(0,len(sec)):
                    tmpsec = float(sec[j])
                    intsec.append(str(int(np.floor(tmpsec))))
        
        
        ind = line.find('#samples per second')
        if ind >= 0:
            sps = line.split('\t')
            sps = sps[0:-1]
            sps = [round(float(x)) for x in sps]
            #print(sps)
        
        ind = line.find('#Counts/Volt')
        if ind >= 0:
            cntpvolt = line.split('\t')
            cntpvolt = cntpvolt[0:-1]
        
        ind = line.find('#Volts/Unit')
        if ind >= 0:
            sen = line.split('\t')
            sen = sen[0:-1]

        ind = line.find('#Flat Gain')
        if ind >= 0:
            gain = line.split('\t')
            gain = gain[0:-1]
            	
        ind = line.find('#latitude')
        if ind >= 0:
            lats = line.split('\t')
            lats = lats[0:-1]
        
        ind = line.find('#longitude')
        if ind >= 0:
            lons = line.split('\t')
            lons = lons[0:-1]
            	
        ind = line.find('#number of samples')
        if ind >= 0:
            nsampstr = line.split('\t')
            nsampstr = nsampstr[0:-1]
            if nsampstr[-1] == '': # This seems to be something odd with eqWave3!
                    nsampstr = nsampstr[0:-1]
            nsamp = []
            for j in range(0,len(nsampstr)):
                nsamp.append(int(nsampstr[j]))

            # get max nsamp
            maxsamp = max(nsamp)

            # make data array
            data = np.zeros((len(comps), maxsamp))

        # get data for new file
#        ind = line.find('       	       	')
#        if ind >= 0:
#            readdat = 0
        
        # this part reads in the data in raw counts
        if readdat == 1:
            i += 1
            
            if i < maxsamp:
                dat = (line.split('\t'))
                for j in range(0,len(nsamp)):
                    if i >= nsamp[j]:
                        data[j, i] = np.nan
                    else:
                        data[j, i] = int(round(float(dat[j].strip('\n'))))
                        

        ind = line.find('--------')
        if ind >= 0:
            readdat = 1
            #print('Reading data...'

    # make date array
    alldatestr = []
    for j in range(0,len(sec)):
        #alldatestr.append(ymd[j]+hhmm[j]+intsec[j])
        alldatestr.append(ymd[j]+hhmm[j])
    
    # set additional data
    tmptime = str(ymd[0])+str(hhmm[0])+str(tmpsec)
    tmptime = datetime.strptime(tmptime,"%Y%m%d%H%M%S.%f")
    if lats == -9999:
        from data_fmt_tools import return_sta_data
        austa = return_sta_data(sta[0])
        lats = austa['stla']
        lons = austa['stlo']
        
    datadict = {'net':'MEL', 'datetime':tmptime, 'lats': lats, 'lons':lons}
#
    
    # reformat kelunji components
    # assume no BB
    i = 0
    #print(comps)
    #print(sps)
    for comp in comps:
        if comp.startswith('EH') or comp.startswith('v') or comp.startswith('EL') or comp.endswith(' t') \
           or comp.endswith(' T') or comp.endswith(' T 5') or comp.startswith('SP') \
           or comp.startswith('c01') or comp.startswith('c02') or comp.startswith('c03') \
           or comp.endswith('vel') or comp.endswith('v') or comp.endswith(' 1') \
           or comp.endswith(' 2') or comp.endswith(' 3') or comp.startswith('HH') \
           or comp.startswith('U Tran'):
            it = 'H'
            g = 'E'
            #print('vel', comp)
        
        if comp.find('Acc') >= 0 or comp.endswith('a') or comp.find(' A ') >= 0 \
           or comp.find('c04') >= 0 or comp.find('c05') >= 0 or comp.find('c06') >= 0 \
           or comp.startswith('HN') or comp.startswith('e 2') or comp.startswith('n 2') \
           or comp.startswith('x DIG1') or comp.startswith('y DIG1') or comp.startswith('z DIG1') \
           or comp.startswith('x 1') or comp.startswith('y 1') or comp.startswith('z 1') \
           or comp.startswith('v 2') or comp.endswith('acc') or comp.startswith('acc'):
            it = 'N'
            g = 'H'
            #print('accel', comp)
            
        else:
            it = 'H'
        
        if it == 'H' and int(sps[i]) > 80:
            g = 'E'
        else:
            g = 'S'

        if comp[0] == 'H' or comp[0] == 'B':
            g = comp[0]

        if comp[1] == 'N':
            it = comp[1]

        # now get orientation
        if comp.find('z ') >= 0 or comp.find('v ') >= 0 or comp.find('Up') >= 0 \
            or comp.find('u ') >= 0 or comp.find('U Tran') >= 0 or comp.find('Z ') >= 0 \
            or comp.find('HHZ') >= 0 or comp.find('BHZ') >= 0 or comp.find('HNZ') >= 0 \
            or comp.find('EHZ') >= 0 or comp.startswith('c03') or comp.startswith('BNZ'):
             o = 'Z'
        elif comp.find('x ') >= 0 or comp.find('e ') >= 0 or comp.find('East') >= 0 \
            or comp.find('E Tran') >= 0 or comp.find('E ') >= 0 or comp.find('X ') >= 0 \
            or comp.find('HHE') >= 0 or comp.find('BHE') >= 0 or comp.find('HNE') >= 0 \
            or comp.find('EHE') >= 0 or comp.startswith('c01') or comp.startswith('BNE') \
            or comp.startswith('HN1'):
             o = 'E'
        elif comp.find('y ') >= 0 or comp.find('n ') >= 0 or comp.find('North') >= 0 \
            or comp.find('N Tran') >= 0 or comp.find('N ') >= 0 or comp.find('Y ') >= 0 \
            or comp.find('HHN') >= 0 or comp.find('BHN') >= 0 or comp.find('HNN') >= 0 \
            or comp.find('EHN') >= 0 or comp.startswith('c02') or comp.startswith('BNN') \
            or comp.startswith('HN2'):
             o = 'N'
        else:
             o = 'U' # Unknown

        comps[i] = g + it + o
        #print(comps
        
        i += 1
    
    return sta, comps, alldatestr, sec, sps, data, nsamp, cntpvolt, sen, gain, datadict

# this function reads nmx ascii format
def readnmx(wavfile):
    # READ HEADER INFO
    print('Reading header info...')
    header = open(wavfile).readlines()

    readdat = 0
    i = 0
    sta = []
    comps = []
    sps = []
    datestr = []
    nsamp = []
    for line in header:
#        line = line.rstrip(' \t\n\r')
        # get station name and get components
        ind = line.find('StnLocChn:')
        if ind >= 0:
            stacomp = line.strip().split(' ')
            sta.append(stacomp[1])
            comps.append(stacomp[-1])

        # get start time
        ind = line.find('Start Time:')
        if ind >= 0:
            y = line[line.find(':') + 2:line.find(':') + 6].strip()
            m = line[line.find(':') + 7:line.find(':') + 9].strip()
            d = line[line.find(':') + 10:line.find(':') + 12].strip()
            hh =line[line.find(':') + 13:line.find(':') + 15].strip()
            mm =line[line.find(':') + 16:line.find(':') + 18].strip()
            sec =line[line.find(':') + 19:].strip()
            intsec = str(int(np.floor(float(sec))))
            ymd = y+m+d
            hhmm = hh+mm
            datestr.append(ymd+hhmm+intsec)
        ind = line.find('Sample Rate:')
        if ind >= 0:
            sps.append(int(float(line[line.find(':') + 2:].strip())))

        ind = line.find('Number of Samples:')
        if ind >= 0:
            nsamp.append(int(float(line[line.find(':') + 2:].strip())))
            # make data array
            data = np.zeros((nsamp[0],1))
#        ind = line.find('DC Offset:')
#        if ind >= 0:
#            dc_Offset = int(float(line[line.find(':') + 2:].strip()))

        # get data for new file
        # this part reads in the data in raw counts
        if readdat == 1:
            if i < nsamp[0]:
                dat = line.split(',')
                col = len(dat) - 1

                for j in range(0,col):
                    data[i+j] = int(dat[j].strip())
                i += 5

        ind = line.find('Format Version: 5.0')
        if ind >= 0:
            readdat = 1
            #print('Reading data...'

    return sta, comps, datestr, sec, sps, data, nsamp

"""this function reads simple TSPAIR format
TSPAIR can be output using obspy functions for any stream data, eg:
st.write("tspair.ascii", format="TSPAIR")
http://docs.obspy.org/packages/autogen/obspy.core.ascii.writeTSPAIR.html"""
def readtspair(wavfile):
    from datetime import datetime
    print('Reading header info...')
    header = open(wavfile).readlines()

    readdat = 0
    i = 0
    nsamp = []
    sta = []
    comps = []
    sps = []
    datestr = []
    for line in header:
        # get header info
        if i == 0:
            headinfo = line.strip('\n').split(',')

            # get station name
            tmpsta = headinfo[0].split('_')
            sta.append(tmpsta[1])

            # get component
            comps.append(tmpsta[-2])

            # get nsamples
            tmpsamp = headinfo[1].strip().split()
            nsamp.append(int(tmpsamp[0]))

            # make data array
            data = np.zeros((nsamp[0],1))

            # get sample rate
            tmpsps = headinfo[2].strip().split()
            sps.append(int(tmpsps[0]))

            # get record time
            tmptime = headinfo[3].strip()
            tmptime = datetime.strptime(tmptime,"%Y-%m-%dT%H:%M:%S.%f")
            sec = tmptime.strftime("%S")
            datestr.append(tmptime.strftime("%Y%m%d%H%M%S"))

        # now read data
        else:
            dat = line.strip('\n').split()
            data[i-1] = int(float(dat[1]))

        i += 1


    return sta, comps, datestr, sec, sps, data, nsamp

def readseismac(wavfile):
    # READ HEADER INFO
    from datetime import datetime

    #print('\nReading header info...'
    header = open(wavfile, 'rU')
#    header = header[0]

    readdat = 0

    i = -1
    for line in header:
#        line = line.rstrip(' \t\n\r')
#        print(line
        # get station name
        ind = line.find('SiteName=')
        if ind >= 0:
            sta = line.strip('\n').split('=')
            sta = sta[1]
        
        ind = line.find('Operator=')
        if ind >= 0:
            net = line.strip('\n').split('=')
            net = net[1]
            
        # get components
        ind = line.find('TIME')
        if ind >= 0:
            comps = line.strip('\n').split('\t')
            comps = comps[1:]

        # get start time
        ind = line.find('StartTime=')
        if ind >= 0:
            ymd = line.strip('\n').split('=')
            ymd = ymd[1]
            print(ymd)

        ind = line.find('sps')
        if ind >= 0:
            sps = line.strip('\n').split('sps\t')

        ind = line.find('RecordLength=')
        if ind >= 0:
            nsampstr = line.split('=')
            nsamp = int(nsampstr[1])

        # this part reads in the data in raw counts
        if readdat == 1:
            i += 1
            if i < nsamp:
                dat = (line.split('\t'))
                for j in range(0,len(comps)):
                    if i >= nsamp:
                        data[j, i] = np.nan
                    else:
                        data[j, i] = int(dat[j+1].strip('\n'))

        ind = line.find('************')
        if ind >= 0:
            readdat = 1

            # make data array
            data = np.zeros((len(comps), nsamp))

            #print('Reading data...'

    # make date array
    allnsamp = []
    alldatestr = []
    allsta = []
    allsec = []
    tmptime = datetime.strptime(ymd,"%Y-%m-%d %H%M %S.%f")
    #print(comps)
    sec = tmptime.strftime("%S")
    for j in range(0,len(comps)):
        alldatestr.append(tmptime.strftime("%Y%m%d%H%M%S")) # 2022-02-21 edit
        alldatestr.append(tmptime)
        allnsamp.append(nsamp)
        allsta.append(sta)
        allsec.append(sec)
    
    # set additional data
    if net == 'ASC':
        net = 'AU'
    elif net == 'SRC':
        net = 'MEL'
    else:
        net = 'MEL'
    if sta[0] == 'NLD' or sta[0] == 'RIV' or sta[0] == 'KIM':
       net = 'AU'
    datadict = {'net':net, 'datetime':tmptime}
    
    # reformat kelunji components
    # assume no BB
    i = 0
    for comp in comps:
        ind = comp.find('Acc')
        if ind >= 0:
            it = 'N'
            g = 'H'
        elif comp.endswith('4') or comp.endswith('5') or comp.endswith('6'):
            it = 'N'
            g = 'H'
        else:
            it = 'H'
            if int(sps[i]) > 80:
                g = 'E'
            else:
                g = 'S'

        if comp[0] == 'H' or comp[0] == 'B':
            g = comp[0]

        if comp[1] == 'N':
            it = comp[1]

        # now get orientation
        if comp.find('z ') >= 0 or comp.find('v ') >= 0 or comp.startswith('Up') \
            or comp.find('u ') >= 0 or comp.find('U Tran') >= 0 or comp.find('Z ') >= 0 \
            or comp.find('HHZ') >= 0 or comp.find('BHZ') >= 0 or comp.find('HNZ') >= 0 \
            or comp.startswith('vertical') or comp.endswith('3') or comp.endswith('6'):
             o = 'Z'
             #print('up?')
        elif comp.find('x ') >= 0 or comp.find('e ') >= 0 or comp.startswith('East') \
            or comp.find('E Tran') >= 0 or comp.find('E ') >= 0 or comp.find('X ') >= 0 \
            or comp.find('HHE') >= 0 or comp.find('BHE') >= 0 or comp.find('HNE') >= 0 \
            or comp.startswith('east') or comp.endswith('1') or comp.endswith('4'):
             o = 'E'
             #print('east?')
        elif comp.find('y ') >= 0 or comp.find('n ') >= 0 or comp.startswith('North') \
            or comp.find('N Tran') >= 0 or comp.find('N ') >= 0 or comp.find('Y ') >= 0 \
            or comp.find('HHN') >= 0 or comp.find('BHN') >= 0 or comp.find('HNN') >= 0 \
            or comp.startswith('north') or comp.endswith('2') or comp.endswith('5'):
             o = 'N'
             #print('north?')
        else:
             o = 'U' # Unknown

        comps[i] = g + it + o
        #print(comps)

        i += 1

    return allsta, comps, alldatestr, allsec, sps, data, allnsamp, datadict

# here we use obspy to read miniseed files
#def readseed(allst):
#    # make data array
#    data = np.zeros((maxsamp,len(comps)))
#    i = 0
#    for st in allst:
#        data[0,??] = st[0].data

# read data from BKN deployement
def readbkn(wavfile):
    from numpy import ones_like, zeros, array
    
    # for testing
    # bknfile = 'U:\\Geoscience_Australia\\eqSource\\Waves\\2002\\03\\05_0147\\BK1_02064_0147.G'
    
    lines = open(wavfile).readlines()
    
    # get date str
    dat = lines[0].strip().split()
    datestr = dat[2].replace('/','') + dat[3][0:6].replace(':','')
    sec = float(lines[0].strip().split(':')[-1])
    
    # get station
    dat = lines[2].strip().split()
    sta = dat[1]
    
    # get n samples and sps
    try:
        dat = lines[12].strip().split()
        nsamp = int(dat[0])
        sps = int(dat[3])
        
        # get comps
        dat = lines[16].strip().split()
        comps = dat[2:]
    except:
        dat = lines[7].strip().split()
        nsamp = int(dat[0])
        sps = int(dat[3])
        
        # get comps
        dat = lines[11].strip().split()
        comps = dat[2:]    

    # finish comps    
    for i, comp in enumerate(comps):
        if comp[0] == 'G' or sta == 'GOK' or sta == 'CMC' or sta == 'KOO1' \
           or sta == 'KOO2':
            comps[i] = 'HN' + comp[-1]
            #print(sta, 'Acc', comp, comps[i], wavfile
        elif comp[0] == 'S' and sps > 80:
            comps[i] = 'EH' + comp[-1]
            #print(sta, 'Vel', comp, comps[i], wavfile
            
    # get sensitivity
    gain = 1.0
    cntpvolt = 1.0
    instrument = lines[6].split()[-1]
    try:
        try:
            zsen = 1.0 / (float(lines[7].split()[1]) * 1E-9)
            esen = 1.0 / (float(lines[11].split()[1]) * 1E-9)
            nsen = 1.0 / (float(lines[9].split()[1]) * 1E-9)
            sen = array([esen, nsen, zsen])
            
        except:
            sen = array([1.0, 1.0, 1.0])
    
    # for 2004 GOK and CMC
    except:
        sen = array([1.0, 1.0, 1.0])
        
    
    # now get data
    readdat = False
    i = 0
    data = zeros((len(comps), nsamp))
    
    for line in lines:
        if readdat == True:
            dat = line.strip().split()
            data[0,i] = float(dat[1])
            data[1,i] = float(dat[2])
            data[2,i] = float(dat[3])
            
            # get second offset
            if i == 0:
                secoffset = float(dat[0])

            i += 1
        else:
            #print(line
            if line.startswith('--------------'):
                readdat = True
        
    # now make lists
    allsta = []
    alldatestr = []
    allnsamp = []
    for comp in comps:
        allsta.append(sta)
        alldatestr.append(datestr)
        allnsamp.append(nsamp)
        
    sps = ones_like(allnsamp) * sps
    allsec = ones_like(allnsamp) * sec
    
    return allsta, comps, alldatestr, allsec, sps, data, allnsamp, cntpvolt, sen, gain, instrument, secoffset

# here we use obspy to read miniseed files
def readseed(st):
    from sys import argv 
    from numpy import hstack, zeros_like, shape
    from data_fmt_tools import remove_low_sample_data
    
    st = remove_low_sample_data(st)
        
    # make data arrays
    allsta = []
    comps = []
    alldatestr = []
    allsec = []
    allsps = []
    alldata = []
    allnsamp = []
    allnetid = []
    
    # first merge data gaps
    st.merge()
    
    # get min samples
    min_npts = 10**20
    for tr in st:
        if tr.stats['npts'] < min_npts:
            min_npts = tr.stats['npts']
            
    # get max samples
    max_npts = 0
    for tr in st:
        if tr.stats['npts'] > max_npts:
            max_npts = tr.stats['npts']
            maxarray = tr.data
    
    # populate data arrays
    for i, tr in enumerate(st):
        allsta.append(tr.stats['station'])
        alldatestr.append(tr.stats['starttime'].strftime("%Y%m%d%H%M%S"))
        allsec.append(tr.stats['starttime'].strftime("%S"))
        allsps.append(tr.stats['sampling_rate'])
        allnsamp.append(tr.stats['npts'])
        allnetid.append(tr.stats['network'])
        #allnsamp.append(min_npts) # not sure why this line was added!
              
        # get channel 
        #print(tr
        if tr.stats['channel'] == 'HHE' or tr.stats['channel'] == 'HNE' or \
           tr.stats['channel'] == 'BHE' or tr.stats['channel'] == 'ENE' or \
           tr.stats['channel'] == 'EHE' or tr.stats['channel'] == 'SHE' or \
           tr.stats['channel'] == 'SNE' or tr.stats['channel'] == 'BNE' or \
           tr.stats['channel'] == 'HN1' or tr.stats['channel'] == 'HLE' or \
           tr.stats['channel'] == 'CPE' or tr.stats['channel'] == 'EYE' or \
           tr.stats['channel'] == 'CHE':
            comp = tr.stats['channel']
        elif tr.stats['channel'] == 'HHN' or tr.stats['channel'] == 'HNN' or \
             tr.stats['channel'] == 'BHN' or tr.stats['channel'] == 'ENN' or \
             tr.stats['channel'] == 'EHN' or tr.stats['channel'] == 'SHN' or \
             tr.stats['channel'] == 'SNN' or tr.stats['channel'] == 'BNN' or \
             tr.stats['channel'] == 'HN2' or tr.stats['channel'] == 'HLN' or \
             tr.stats['channel'] == 'CPN' or tr.stats['channel'] == 'EYN' or \
             tr.stats['channel'] == 'CHN':
            comp = tr.stats['channel']
        elif tr.stats['channel'] == 'HHZ' or tr.stats['channel'] == 'HNZ' or \
             tr.stats['channel'] == 'BHZ' or tr.stats['channel'] == 'ENZ' or \
             tr.stats['channel'] == 'EHZ' or tr.stats['channel'] == 'SHZ' or \
             tr.stats['channel'] == 'SNZ' or tr.stats['channel'] == 'BNZ' or \
             tr.stats['channel'] == 'HLZ' or tr.stats['channel'] == 'CPZ' or \
             tr.stats['channel'] == 'EYZ'  or tr.stats['channel'] == 'CHZ':
            comp = tr.stats['channel']
        else:
            numtrue = False
            print(tr.stats['channel'])
            if tr.stats['channel'] == '001' or tr.stats['channel'] == '004' or tr.stats['channel'] == 'HH1' or tr.stats['channel'] == 'BH1' or tr.stats['channel'] == 'HLE':
                o = 'E'
                it = 'H'
                numtrue = True
            elif tr.stats['channel'] == '002' or tr.stats['channel'] == '005' or tr.stats['channel'] == 'HH2'or tr.stats['channel'] == 'BH2' or tr.stats['channel'] == 'HLN':
                o = 'N'
                it = 'H'
                numtrue = True
            if tr.stats['channel'] == 'HN1':
                o = 'H'
                it = 'N'
                numtrue = True
            if tr.stats['channel'] == 'HN2':
                o = 'H'
                it = 'N'
                numtrue = True
            elif tr.stats['channel'] == '003' or tr.stats['channel'] == '006' or tr.stats['channel'] == 'HLZ':
                o = 'Z'
                numtrue = True
                
            if tr.stats['sampling_rate'] >= 80:
                g = 'E'
                if tr.stats['channel'].startswith('H'):
                    g = 'H'
            else:
                if tr.stats['channel'].startswith('B'):
                    g = 'B'
                else:
                    g = 'S'
                
            if numtrue == True:
                if int(tr.stats['channel'][2]) > 3 or tr.stats['channel'][1] == 'N':
                    it = 'N'
                else:
                    it = 'H'
            comp = g + it + o
        
        # append comp
        comps.append(comp) 
        
        tmpdat = zeros_like(maxarray)
        #print(tr.stats['station'], len(maxarray), tr.stats['npts']
        tmpdat[0:tr.stats['npts']] = tr.data
        #print(len(tmpdat)
        if i == 0:
            alldata = tmpdat.reshape(max_npts,1)
        else:
            alldata = hstack((alldata, tmpdat.reshape(max_npts,1)))
            
    #print(shape(alldata)
    return allsta, comps, alldatestr, allsec, allsps, alldata, allnsamp, allnetid

# make text to select channel
def select_channel(sta, comps):
    # if only one channel, ignore
    if len(comps) == 1:
        selchan = 0
        i = 0
    else:
        chantxt = '\n'+'Select channel: \n'
        for i in range(0,len(comps)):
            chantxt = chantxt + str(i+1) + ') ' + sta[i] + ' ' + comps[i] + '\n'
        chantxt = chantxt + '\n' + 'Channel number > '
        try:
            selchan = int(raw_input(chantxt)) - 1
        except:
            selchan = int(input(chantxt)) - 1
    return comps[selchan], selchan

# REMOVE DC OFFSET
def remove_dc_offset(data):
    data[0,:] = np.round(data[0,:] - np.median(data[0,:]))

    return data
    
