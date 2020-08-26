# -*- coding: utf-8 -*-
"""
Created on Mon Aug 06 17:14:54 2012

@author: u56903
"""
import os
import numpy as np

# use info collected to get systematic output filename
def get_filename(sta, recdate, chan, inst_ty, dt, sps):
    import datetime

    # get event time string
    evdate = recdate + datetime.timedelta(minutes=int(dt/60))
    datestr = evdate.strftime("%Y%m%d%H%M")

    # get inst code
    if inst_ty == 'N':
        inst_ty = 'H'
        inst_code = 'N'
    else:
        inst_code = 'H'

    # get orientation
    if chan.find('z ') >= 0 or chan.find('v ') >= 0 or chan.find('Up') >= 0 \
       or chan.find('u ') >= 0 or chan.find('BHZ') >= 0 or chan.find('SHZ') >= 0  \
       or chan.find('EHZ') >= 0 or chan.find('HNZ') >= 0 or chan.find('HHZ') >= 0 \
       or chan.find('U Tran') >= 0 or chan.endswith('Z'):
        orient = 'Z'
    elif chan.find('x ') >= 0 or chan.find('e ') >= 0 or chan.find('East') >= 0 \
       or chan.find('BHE') >= 0 or chan.find('SHE') >= 0 or chan.endswith('E')  \
       or chan.find('EHE') >= 0 or chan.find('HNE') >= 0 or chan.find('HHE') >= 0:
        orient = 'E'
    elif chan.find('y ') >= 0 or chan.find('n ') >= 0 or chan.find('North') >= 0 \
       or chan.find('BHN') >= 0 or chan.find('SHN') >= 0 or chan.endswith('N') \
       or chan.find('EHN') >= 0 or chan.find('HNN') >= 0 or chan.find('HHN') >= 0:
        orient = 'N'
    else:
        orient = 'U' # Unknown

    # get final filename
    filename = ''
    joinstr = (datestr, '.', sta, '.', inst_ty,  inst_code, orient)
    filename = filename.join(joinstr)

    return filename, evdate

# this function gets text for the file header
def get_header_text(filename, sta, evdate, sps, stla, stlo, eqla, eqlo, eqdep, \
                    eqmag, rhyp, azim, pga, pgv, lofreq, hifreq, filetype):
    # get processing time
    import datetime as dt
    proctime = dt.datetime.now()

    tfile =  'FILE:\t' + filename
    tdate =  'DATE:\t' + evdate.strftime("%Y-%m-%d %H:%M")
    tsta =   'STA:\t' + sta
    tstla =  'STLA:\t' + str(stla)
    tstlo =  'STLO:\t' + str(stlo)
    tsps =   'SR:\t' + str(sps) + ' Hz'
    teqla =  'EQLA:\t' + str(eqla)
    teqlo =  'EQLO:\t' + str(eqlo)
    teqdep = 'EQDEP:\t' + str(eqdep)  + ' km'
    teqmag = 'EQMAG:\t' + str(eqmag)
    trhyp =  'RHYP:\t' + str("%0.1f" % rhyp) + ' km'
    tazim =  'AZIM:\t' + str("%0.1f" % azim) + ' degrees'
    if filetype == 'psa' or filetype == 'cor':
        tpga =   'PGA:\t' + str("%0.3f" % pga) + ' mm/s/s'
        tpgv =   'PGV:\t' + str("%0.3f" % pgv) + ' mm/s\n'
    proc =   'Record Processed: ' + proctime.strftime("%Y-%m-%d %H:%M") + '\n'
    buff =   '-----------------------------------------------------------------'
    if filetype == 'cor':
        comm = 'Instrument corrected time histories filtered using a 4th order\n' + \
               'Butterworth bandpass between ' + str(lofreq) + '-' \
               + str(hifreq) + ' Hz\n'
        units = 'Disp (m)\tVel (m/s)\tAcc (m/s/s)'

    elif filetype == 'fds':
        comm = '\nDisplacement Fourier spectra claculated using time histories \n' + \
               'filtered using a 4th order Butterworth bandpass between \n' + \
                str(lofreq) + '-' + str(hifreq) + ' Hz\n'
        units = 'Frequency (Hz)\tFDS (m-s)'

    elif filetype == 'psa':
        comm = '5% Damped Pseudo Response Spectral Acceleration (PSA) claculated\n' + \
               'using time histories filtered using a 4th order Butterworth\n' + \
               'bandpass between ' + str(lofreq) + '-' + str(hifreq) + ' Hz\n'
        units = 'Period (s)\tPSA (m/s/s)'

    elif filetype == 'wa':
        comm = '\nSynthesised Wood-Anerson displacement seismogram claculated\n' + \
               'using time histories filtered using a 4th order Butterworth\n' + \
               'bandpass between ' + str(lofreq) + '-' + str(hifreq) + ' Hz\n'
        units = 'W-A Disp (mm)'

    header = '\n'
    if filetype == 'psa' or filetype == 'cor':
        joinstr = (tfile, tdate, tsta, tstla, tstlo, tsps, teqla, teqlo, teqdep, \
                   teqmag, trhyp, tazim, tpga, tpgv, comm, proc, buff, units, buff)
    else:
        joinstr = (tfile, tdate, tsta, tstla, tstlo, tsps, teqla, teqlo, teqdep, \
                   teqmag, trhyp, tazim, comm, proc, buff, units, buff)
    header = header.join(joinstr) + '\n'

    return header

# write output time histories
def write_waves(sta, evdate, sps, idisp, ivel, iacc, filename, stla, stlo, \
                eqla, eqlo, eqdep, eqmag, rhyp, azim, lofreq, hifreq):

    # get PGA for file header in mm/s/s
    pga = max(abs(iacc))*1000

    # get PGV for file header in mm/s
    pgv = max(abs(ivel))*1000

    # set file and path names
    filename = filename + '.cor'
    outdir = 'cor'
    outfile = os.path.join(outdir,filename)

    # set header info
    header = get_header_text(filename, sta, evdate, sps, stla, stlo, eqla, eqlo, eqdep, \
                             eqmag, rhyp, azim, pga, pgv, lofreq, hifreq, 'cor')

    # set data array
    data = np.vstack((idisp, ivel, iacc))

    # try saving to cor directory
    try:
        f = open(outfile,'w')
    except:
        os.mkdir('cor')
        f = open(outfile,'w')

    # write header
    f.write(header)
    f.close()

    # now append data
    f = open(outfile, 'a')
    np.savetxt(f, data.T, delimiter='\t', fmt='%1.6e')
    f.close()

# this function writes out the FFT data
def write_fft(sta, evdate, sps, freq, corfftr, corffti, filename, stla, stlo, \
              eqla, eqlo, eqdep, eqmag, rhyp, azim, lofreq, hifreq, inst_ty):

    #n = len(freq[0]) / 2
    n = int(np.floor(len(freq) / 2))

    # set file and path names
    filename = filename + '.fds'
    outdir = 'fds'
    outfile = os.path.join(outdir,filename)

    # set header info
    header = get_header_text(filename, sta, evdate, sps, stla, stlo, eqla, eqlo, eqdep, \
                             eqmag, rhyp, azim, -12345, -12345, lofreq, hifreq, 'fds')
    
    '''
    (filename, sta, evdate, sps, stla, stlo, eqla, eqlo, eqdep, \
                    eqmag, rhyp, azim, pga, pgv, lofreq, hifreq, filetype)
    '''
    # calculate displacement spectra
    dispamp = np.zeros(np.shape(freq))
    dispfftr = np.zeros(np.shape(freq))
    dispffti = np.zeros(np.shape(freq))

    if inst_ty == 'N': # Assume accelerometer
        dispfftr[1:] = corfftr[1:]/(2.*(np.pi)*freq[1:])**2
        dispffti[1:] = corffti[1:]/(2.*(np.pi)*freq[1:])**2
        dispamp[1:] = np.sqrt((1./sps)**2 * (dispfftr[1:]**2 + dispffti[1:]**2))
    else: # assume seismometer
        '''
        dispfftr[0,1:] = corfftr[0,1:]/(2.*(np.pi)*freq[0,1:])
        dispffti[0,1:] = corffti[0,1:]/(2.*(np.pi)*freq[0,1:])
        dispamp[0,1:] = np.sqrt((1./sps)**2 * (dispfftr[0,1:]**2 + dispffti[0,1:]**2))
        '''
        dispfftr[1:] = corfftr[1:]/(2.*(np.pi)*freq[1:])
        dispffti[1:] = corffti[1:]/(2.*(np.pi)*freq[1:])
        dispamp[1:] = np.sqrt((1./sps)**2 * (dispfftr[1:]**2 + dispffti[1:]**2))

    # set data array
    #data = np.vstack((freq[0,1:n],dispamp[0,1:n]))
    data = np.vstack((freq[1:n],dispamp[1:n]))

    # try saving to fft directory
    try:
        f = open(outfile,'w')
    except:
        os.mkdir('fds')
        f = open(outfile,'w')

    # write header
    f.write(header)
    f.close()

    # now append data
    f = open(outfile, 'a')
    np.savetxt(f, data.T, delimiter='\t', fmt='%1.6e')
    f.close()


# write output response spectra
def write_response_spectra(sta, evdate, sps, T, psa, pga, pgv, filename, stla, stlo, \
                          eqla, eqlo, eqdep, eqmag, rhyp, azim, lofreq, hifreq):

    # get PGA for file header in mm/s/s
    pga *= 1000

    # get PGV for file header in mm/s
    pgv *= 1000

    # set file and path names
    filename = filename + '.psa'
    outdir = 'psa'
    outfile = os.path.join(outdir,filename)

    # set header info
    header = get_header_text(filename, sta, evdate, sps, stla, stlo, eqla, eqlo, eqdep, \
                             eqmag, rhyp, azim, pga, pgv, lofreq, hifreq, 'psa')

    # set data array
    data = np.vstack((T, psa))

    # try saving to spc directory
    try:
        f = open(outfile,'w')
    except:
        os.mkdir('psa')
        f = open(outfile,'w')

    # write header
    f.write(header)
    f.close()

    # now append data
    try:
        f = file(outfile, 'a')
    except:
        f = open(outfile, 'a')
    np.savetxt(f, data.T, delimiter='\t', fmt='%1.5e')
    f.close()

# this function writes Wood-Anderson displacement (in mm) to file
def write_WA_disp(sta, evdate, sps, wadisp, filename, stla, stlo, \
                  eqla, eqlo, eqdep, eqmag, rhyp, azim, lofreq, hifreq):

    # set file and path names
    filename = filename + '.wa'
    outdir = 'wa'
    outfile = os.path.join(outdir,filename)

    # set header info
    '''
    filename, sta, evdate, sps, stla, stlo, eqla, eqlo, eqdep, \
                             eqmag, rhyp, azim, pga, pgv, lofreq, hifreq,
    '''                             
    header = get_header_text(filename, sta, evdate, sps, stla, stlo, eqla, eqlo, eqdep, \
                             eqmag, rhyp, azim, 0, 0, lofreq, hifreq, 'wa')

    # try saving to wa directory
    try:
        f = open(outfile,'w')
    except:
        os.mkdir('wa')
        f = open(outfile,'w')

    # write header
    f.write(header)
    f.close()

    # now append data
    f = file(outfile, 'a')
    np.savetxt(f, wadisp, delimiter='\t', fmt='%1.5e')
    f.close()

# this function writes ML estimates to file
def write_ML_dat(filename, rhyp, repi, logA, R35, BJ84, HB87, GG91, MLM92, A10):
    evstacomp = filename.split('.')
    event = evstacomp[0]
    sta = evstacomp[1]
    comp = evstacomp[2]

    # set file and path names
    filename = event + '.ml'
    outdir = 'ml'
    outfile = os.path.join(outdir,filename)

    header = 'SITE\tCOMP\tRHYP\tREPI\tLOGA\tR35\tBJ84\tHB87\tGG91\tMLM92\tA10\n'
    newtxt = ''

    # try saving to ml directory
    try:
        f = open(outfile,'r')
        f.close()
    except:
        try:
            f = open(outfile,'w')
        except:
            os.mkdir('ml')
            f = open(outfile,'w')
        f.write(header)
        f.close()

    # now write mag estimates for each stn
    newtxt = ''
    newsta = True
    data = open(outfile).readlines()
    for line in data:
        tmp = line.strip('\t')
        # overwrite current station
        if tmp[0] == sta and tmp[1] == comp:
            replacesta = '\t'
            joinstr = (sta, comp,str("%0.1f" % rhyp), str("%0.1f" % repi), str("%0.3f" % logA), \
                       str("%0.2f" % R35), str("%0.2f" % BJ84), str("%0.2f" % HB87), \
                       str("%0.2f" % MLM92), str("%0.2f" % A10))
            replacesta = replacesta.join(joinstr) + '\n'
            newtxt = newtxt + replacesta
            newsta = False
        else:
            # copy existing lines
            newtxt = newtxt + line

    # if new station, append to file
    if newsta == True:
        newstastr = '\t'
        joinstr = (sta, comp,str("%0.1f" % rhyp), str("%0.1f" % repi), str("%0.3f" % logA), \
                   str("%0.2f" % R35), str("%0.2f" % BJ84), str("%0.2f" % HB87), \
                   str("%0.2f" % MLM92), str("%0.2f" % A10))
        newstastr = newstastr.join(joinstr) + '\n'
        newtxt = newtxt + newstastr

    # now overwrite file
    f = open(outfile,'w')
    f.write(newtxt)
    f.close()

# this function uses Obspy modules to export to SAC format
def export_SAC(filename, site, netid, sta, stla, stlo, eqla, eqlo, eqdep, \
               eqmag, sps, iwav):
    try:
        # try importing obspy modules
        from obspy.core import Trace, Stream, UTCDateTime, AttribDict
        from obspy.core.util import calcVincentyInverse

        # set file and path names
        filename = filename + '.sac'
        outdir = 'sac'
        outfile = os.path.join(outdir,filename)

        ymd = filename[0:4]
        hhmm = filename[4:8]
        starttime = ymd+'T'+hhmm+'00.0'

        # try saving to sac directory
        try:
            f = open(outfile,'w')
        except:
            os.mkdir('sac')

        nsamp = len(iwav)

        component = filename.split('.')
        component = component[1]
        stats = {'network': '', 'station': sta, 'location': '',
                 'channel': component, 'npts': nsamp, 'sampling_rate': sps}


        stats['starttime'] = UTCDateTime(starttime)
        sac = AttribDict()
        sac.__setattr__('stlon',stlo)
        sac.__setattr__('stla', stla)
        sac.__setattr__('idep', 1)
        sac.__setattr__('evla', eqla)
        sac.__setattr__('evlo', eqlo)
        sac.__setattr__('evdp', eqdep)

        # calculate distance and azimuth
        diazbaz = calcVincentyInverse(eqla,eqlo,stla,stlo)
        distkm = diazbaz[0] / 1000.
        R = 6371.
        gcarc = distkm * 360 / (2 * np.pi * R)
        az = diazbaz[1]
        baz = diazbaz[2]
        sac.__setattr__('gcarc', gcarc)
        sac.__setattr__('az', az)
        sac.__setattr__('baz', baz)
        sac.__setattr__('dist', distkm)
        stats['sac'] = sac

        st = Stream([Trace(data=iwav, header=stats)])
        st.write(outfile, format='SAC')

    except:
        # print statement if cannot access obspy
        print('\nCannot import Obspy modules!')
