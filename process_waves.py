# -*- coding: utf-8 -*-
"""
Created on Wed Jul 25 15:05:20 2012

@author: Trevor Allen

process_eqwave.py

This code reads ES&S eqWave format text files and performs the following tasks
depending on the users requirements:

    1) Do instrument correction and write output time history
    2) Calculate FFT and write output spectrum
    3) Calculate and output 5% damped response spectrum
    4) Convert to Wood-Anderson time history and calculate local magnitude
    5) Fit Brune spectrum and calculate moment magnitude
    6) Export to miniSEED format
    7) Export to SAC format
    8) Plot instrument response

"""
# IMPORT LIBRARIES
import numpy as np
#import matplotlib.pyplot as plt
import response, readwaves, spectral_analysis, plotting, write_data, spatial_tools
from sys import argv
from os import getcwd, path
from datetime import datetime

#*******************************************************************************
# GET RESPONSE INFORMATION
#def checkpaz():
#    var = raw_input('Use PAZ file (y/n)? > ')
#    if var == 'y':
#        return 1
#    else:
#        return 0

#*******************************************************************************
# SET TASK
wavfile = argv[1]

# get absolute path
cwd = getcwd()
wavfile = path.join(cwd, wavfile)

def task_options():
    tasks = '\nSelect task: \n' \
             + '1) Do instrument correction and write output time history\n' \
             + '2) Calculate FFT and write output spectrum\n' \
             + '3) Calculate and output 5% damped response spectrum\n' \
             + '4) Convert to Wood-Anderson time history and calculate local magnitude\n' \
             + '5) Fit Brune spectrum and calculate moment magnitude\n' \
             + '6) Export to miniSEED format\n' \
             + '7) Export to SAC format\n' \
             + '8) Plot instrument response\n\n' \
             + 'Task number [3] > '
    try:
        selection = raw_input(tasks)
    except:
        selection = input(tasks)
    return selection

"""****************************************************************************
Common functions to read data file and get instrument parameters
****************************************************************************"""

def common_read(allsta, comps, allrecdate, allsec, allsps, alldata, allnsamp, sacseed):
    
    #print(len(alldata[:,17])
    # select component
    chan, chan_no = readwaves.select_channel(allsta, comps)

    # unify kelunji channel names
#    tmpchan = chan.split(' Trans')
#    chan = ''
#    for i in range(0,len(tmpchan)):
#        chan = chan + tmpchan[i]

    # get channel data of interest
    sta = allsta[chan_no].strip()
    sps = int(allsps[chan_no])
    
    if sacseed == True:
        chan_dat = alldata[:,chan_no]
    else:
        chan_dat = alldata[chan_no] # edited as of 2017-09-26 to read new eqwave
    
    tmprecdate = allrecdate[chan_no]
    recdate = datetime.strptime(tmprecdate, "%Y%m%d%H%M%S")

    # remove trailing nans
    nantrue = True
    while nantrue == True:
        if chan_dat[-1] == np.nan:
            chan_dat = chan_dat[0:-1]
        else:
            nantrue = False

    chan_dat = chan_dat[0:allnsamp[chan_no]]
    chan_dat = chan_dat.reshape(1,len(chan_dat))

    # remove DC offset
    chan_dat = readwaves.remove_dc_offset(chan_dat)

    # check to see if response exists
    nat_freq, inst_ty, damping, sen, recsen, gain, pazfile, stlo, stla, netid \
        = response.get_response_info(sta,recdate,chan)

    # if info does not exist, ask for user input
    if nat_freq == -12345 and pazfile == 'NULL':
        inst_ty, mindate, maxdate, stlo, stla, netid, nat_freq, damping, \
            sen, recsen, gain, pazfile = response.enter_response_info(sta, chan, sps)

        # write new data to file
        response.write_response_info(sta, inst_ty, mindate, maxdate, stlo, \
            stla, netid, nat_freq, damping, sen, recsen, gain, chan, pazfile)
    
    return sta, inst_ty, sps, recdate, nat_freq, damping, sen, recsen, gain, \
           chan, chan_no, chan_dat, stlo, stla, pazfile, alldata, netid

"""****************************************************************************
Common FFT functions
****************************************************************************"""
def common_fft(chan_dat, inst_ty, sps, seltask):
    # view wave and trim if necessary
    print(chan_dat)
    start_index, stop_index = plotting.trim_wave(chan_dat, sps, inst_ty, False)
    #taper_dat = chan_dat[0, start_index:stop_index]
    taper_dat = chan_dat[start_index:stop_index]
    taper_dat = taper_dat.reshape(1, len(taper_dat))

    # if doing FFT analysis, trim again to specify exact window
#    if seltask == '2':
#        start_index, stop_index = plotting.trim_wave(taper_dat, sps, inst_ty, True)
#        taper_dat = taper_dat[0,start_index:stop_index]
#        taper_dat = taper_dat.reshape(1, len(taper_dat))

    # taper time-domain before fft
    taper_dat, costaper = spectral_analysis.cosine_taper(taper_dat)
    
    # do acausal filter in time domain using 4th order Butterworth filter
#    lofreq, hifreq, filtwave = spectral_analysis.acausal_butter_filter(inst_ty, sps, taper_dat, seltask)
#    filtwave = taper_dat.reshape(1, len(filtwave))

    # get FFT of data stream for full record
#    freq, wavfft = spectral_analysis.calc_fft(filtwave, sps)
    freq, wavfft = spectral_analysis.calc_fft(taper_dat[0], sps)

    # filter FFT using 4th order Butterworth filter
    lofreq, hifreq, wavfft = spectral_analysis.butter_filter_user(inst_ty, sps, freq, wavfft, seltask)

    # get time offset for filenames
    dt = start_index / float(sps)

    return freq, lofreq, hifreq, wavfft, dt

 #   return freq, wavfft

"""****************************************************************************
Common functions to remove FAP response
****************************************************************************"""

def common_resp(freq, nat_freq, damping, sen, recsen, gain, wavfft, inst_ty):

    # calculate FAP instrument response
    if pazfile == 'NULL':
        real_resp, imag_resp = response.fap_response(freq, nat_freq, damping, \
                                                     sen, recsen, gain, inst_ty)

    # or use PAZ file
    else:
        real_resp, imag_resp = response.paz_response(freq, pazfile, sen, recsen, \
                                                     gain, inst_ty)

    # deconvolve instrument response from record and return corrected FFT
    corfftr, corffti = response.deconvolve_instrument(real_resp, imag_resp, wavfft)

    return corfftr, corffti, real_resp, imag_resp

"""****************************************************************************
# Begin main code
****************************************************************************"""
# set boolean opperators
continue_loop = True

plot_outputs = True
'''
var = raw_input('\nPlot outputs ([y]/n)? > ')
if var == 'y' or var == '':
    plot_outputs = True
elif var == 'n':
    plot_outputs = False
'''
# first check to see if obspy installed
try:
    # try importing obspy modules
    from obspy.core import read

    print('\nObsPy Installed')

    # try testing if SAC or miniSEED
    try:
        st = read(wavfile)
        sacseed = True

    except:
        sacseed = False

except:
    # print(statement if cannot access obspy
    print('\nCannot import ObsPy modules!')

# if SAC or SEED read waves
if sacseed == True:
    allsta, comps, allrecdate, allsec, allsps, alldata, allnsamp = readwaves.readseed(st)
    fmt = 'obspy'

# esle read text file format
else:
    # check text file format
    fmt = readwaves.check_file_fmt(wavfile)
    if fmt == 'eqw':
        # read the eqWave text file
        try:
            allsta, comps, allrecdate, allsec, allsps, alldata, allnsamp, cntpvolt, allsen, allgain = readwaves.readeqwave(wavfile)
        except:
            allsta, comps, allrecdate, allsec, allsps, alldata, allnsamp = readwaves.readeqwave(wavfile)
    elif fmt == 'nmx':
        allsta, comps, allrecdate, allsec, allsps, alldata, allnsamp = readwaves.readnmx(wavfile)
    elif fmt == 'tspair':
        allsta, comps, allrecdate, allsec, allsps, alldata, allnsamp = readwaves.readtspair(wavfile)
    elif fmt == 'sm':
        allsta, comps, allrecdate, allsec, allsps, alldata, allnsamp = readwaves.readseismac(wavfile)
    else:
        '\nFile format not recognised!'

while continue_loop == True:
    # write options to screen for selection
    seltask = task_options()

    if seltask == '':
        seltask = '3' # response spectra

    if seltask != '8':
        # do common read functions
        sta, inst_ty, sps, recdate, nat_freq, damping, sen, recsen, gain, chan, \
                chan_no, chan_dat, stlo, stla, pazfile, alldata, netid = \
                common_read(allsta, comps, allrecdate, allsec, allsps, alldata, allnsamp, sacseed)

        print(chan)
        # do common fft functions
        chan_dat = chan_dat[0]
        freq, lofreq, hifreq, wavfft, dt = common_fft(chan_dat, inst_ty, sps, seltask)
        
        # do common instrument deconvolution functions
        corfftr, corffti, real_resp, imag_resp = common_resp(freq, nat_freq, damping, sen, \
                                                 recsen, gain, wavfft, inst_ty)
        
        # get filename
        filename, evdate = write_data.get_filename(sta, recdate, chan, inst_ty, dt, sps)

        # get source-to-site distance
        rhyp, azim, eqla, eqlo, eqdep, eqmag, evdate = spatial_tools.get_eq_distance(stlo, stla, evdate)

    if seltask == '1': # plot instrument corrected time histories
        if plot_outputs == True:
            # plot displacement, velocity and acceleration time history
            idisp, ivel, iacc = plotting.plot_dva(freq, sps, corfftr, corffti, \
                                                  filename, inst_ty, chan_no)

        # code for testing
    #    vfft= np.fft.fft(ivel)
    #    freq = np.fft.fftfreq(len(vfft), d=1./sps)
    #    velamp = np.sqrt((1./sps)**2 * (vfft.real**2 + vfft.imag**2))
    #    plt.loglog(freq[1:len(vfft)/2],velamp[1:len(vfft)/2])

        # write time histories to file
        write_data.write_waves(sta, evdate, sps, idisp, ivel, iacc, \
                               filename, stla, stlo, eqla, eqlo, eqdep, eqmag, \
                               rhyp, azim, lofreq, hifreq)

    elif seltask == '2': # calculate FFT

        # plot Fourier displacement spectra
        plotting.plot_fft(sps, freq, corfftr, corffti, filename, inst_ty)

        # write time histories to file
        write_data.write_fft(sta, evdate, sps, freq, corfftr, corffti, \
                               filename, stla, stlo, eqla, eqlo, eqdep, eqmag, \
                               rhyp, lofreq, hifreq, inst_ty)

    elif seltask == '3': # calculate response spectrum

        # prepare data for response spectrum calculation
        pgv, iacc, ivel = spectral_analysis.prep_psa(corfftr, corffti, freq, inst_ty)

        # calculate 5% damped response spectra
        T, psa, pga = spectral_analysis.calc_response_spectra(iacc.real, sps, 5.0, 1./hifreq, 1./lofreq)

        if plot_outputs == True:
            # plot response spectra
            plotting.plot_response_spectra(T, psa, pga, filename, chan_no)

        # write output to file
        write_data.write_response_spectra(sta, evdate, sps, T, psa, pga, pgv, \
                 filename, stla, stlo, eqla, eqlo, eqdep, eqmag, rhyp, azim, lofreq, hifreq)

    elif seltask == '4': # calculate Wood-Anderson seismogram
        import calculate_magnitudes

        # convolve with Wood-Anderson instrument and get displacement wave
        wadisp = response.convolve_WoodAnderson(freq, corfftr, corffti, inst_ty)

#        if plot_outputs == True:
        # plot W-A displacement time history and get windo for ML
        watrim = plotting.plot_WoodAnderson(wadisp, sps, filename, chan_no)

        # calculate magnitudes
        print(filename)
        logA = np.log10(max(abs(watrim)))
        calculate_magnitudes.main(filename, logA, rhyp, eqdep)

        # write output to file
        write_data.write_WA_disp(sta, evdate, sps, wadisp, \
                 filename, stla, stlo, eqla, eqlo, eqdep, eqmag, rhyp, lofreq, hifreq)

    elif seltask == '6': # export as miniSEED
        print('\nModule not yet functional...')

    elif seltask == '7': # export as SAC
        # get velocity time history
        pgv, iacc, ivel = spectral_analysis.prep_psa(corfftr, corffti, freq, inst_ty)

        # now export to SAC format
        if inst_ty == 'N': # if accelerometer
            write_data.export_SAC(filename, sta, netid, sta, stla, stlo, \
                                  eqla, eqlo, eqdep, eqmag, sps, iacc)

        else: # if velocity seismometer
            write_data.export_SAC(filename, sta, netid, sta, stla, stlo, \
                                  eqla, eqlo, eqdep, eqmag, sps, iacc)

    elif seltask == '8': # plot instrument response
        # do common read functions
        sta, inst_ty, sps, recdate, nat_freq, damping, sen, recsen, gain, chan, \
                chan_no, chan_dat, stlo, stla, pazfile, alldata = \
                common_read(allsta, comps, allrecdate, allsec, allsps, alldata, allnsamp)

        # make dummy dataset of size 131072
        if inst_ty == 'H':
            dummydat = np.random.rand(1,131072)
        else:
            dummydat = np.random.rand(1,32768)

        # do FFT to get dummy frequencies
        freq, wavfft = spectral_analysis.calc_fft(dummydat, sps)

        # calculate and plot instrument response
        plotting.plot_instrument_resp(sta, inst_ty, freq, nat_freq, damping, \
                                      sen, recsen, gain, pazfile, chan_no)


    # ask user if users would like to perform another task using current wavfile
    try:
        var = raw_input('\nPerform another task using wavefile: ' + wavfile + ' ([y]/n)? > ')
    except:
        var = input('\nPerform another task using wavefile: ' + wavfile + ' ([y]/n)? > ')
    if var == 'n':
        continue_loop = False
    else:
        continue_loop = True
