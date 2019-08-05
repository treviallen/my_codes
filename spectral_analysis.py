# -*- coding: utf-8 -*-
"""
Created on Sun Jul 29 11:01:37 2012

@author: u56903

This script does various FFT analysis

"""

# this function tapers time-domain window prior to FFT
# uses cosine taper that scales from 0 to 1 fot the first 15% of record length
# and another to scale from 1 to 0 from 85% of the record length
# between 10-90% of the record, taper = 1.0
def cosine_taper(data):
    import numpy as np
#    import matplotlib.pyplot as plt

    r = 0.2
    n = len(data[0])
    m = n*r / 2.0
    x = np.arange(0,n)
    y = np.ones(np.shape(data)) # this is the taper function

    for i in range(0,n):
        i1 = i + 1
        if i <= m-1:
            y[0,i] = 0.5 * i * (1-np.cos(2*np.pi*i / (2*m))) / i1
        elif i >= n-m:
            y[0,i] = 0.5 * i * (1-np.cos(2*np.pi*(n-i-1) / (2*m))) / i1

    # plot taper function
#    plt.figure(2,figsize=(16,9))
#    plt.plot(x,y[0],'-')
#    plt.show()

    return data * y, y

# this function asks for and applies a 4th order bandpass Butterworth filter
def butter_filter_user(inst_ty, sps, freq, wavfft, seltask):
    import numpy as np
    
    filtfft = wavfft
    #print(np.shape(filtfft), 'filtfft'
    
    # filter order
    ford = 4.0

    # set initial hi & lo pass frequencies
    if inst_ty == 'B' or inst_ty == 'N' or inst_ty == 'H':
        lofreq = 0.02
    else:
        lofreq = 0.2

    # override if plotting FFT
    if seltask == '2':
        lofreq = 0.02

    # set low pass at Nyquist frequency
    hifreq = sps / 2.0

    # ask for user input
    txt = '\n'+'High pass Butterworth corner ['+str(lofreq)+' Hz] > '
    var = raw_input(txt)
    if var != '':
        lofreq = float(var)

    txt = '\n'+'Low pass Butterworth corner ['+str(hifreq)+' Hz] > '
    var = raw_input(txt)
    if var != '':
        hifreq = float(var)

    # do high pass filter
    #print(freq, 'freq'
    #print(np.shape(filtfft[1:]), 'filtfft[1:]'
    #filtfft[1:] /= np.sqrt(1 + (abs(freq[0,1:]) / lofreq)**(-2*ford))
    filtfft[1:] /= np.sqrt(1 + (abs(freq[1:]) / lofreq)**(-2*ford))
    
    # do low pass filter
    #filtfft[1:] /= np.sqrt(1 + (abs(freq[0,1:]) / hifreq)**(2*ford))
    filtfft[1:] /= np.sqrt(1 + (abs(freq[1:]) / hifreq)**(2*ford))

    # remover zero freq
    filtfft[0] *= 0.0

    return lofreq, hifreq, filtfft

# this function automatically applies a 4th order bandpass Butterworth filter
# using predefined filter corners
def butter_filter_auto(inst_ty, sps, freq, wavfft):
    import numpy as np
    filtfft = wavfft
    # filter order
    ford = 4.0

    # set initial hi & lo pass frequencies
    if inst_ty == 'B' or inst_ty == 'N' or inst_ty == 'H':
        lofreq = 0.02
    else:
        lofreq = 0.2

#    # override if plotting FFT
#    if seltask == '2':
#        lofreq = 0.02

    # set low pass at Nyquist frequency
    hifreq = sps / 2.0

    # do high pass filter
    filtfft[1:] /= np.sqrt(1 + (abs(freq[0,1:]) / lofreq)**(-2*ford))

    # do low pass filter
    filtfft[1:] /= np.sqrt(1 + (abs(freq[0,1:]) / hifreq)**(2*ford))

    # remover zero freq
    filtfft[0] *= 0.0

    return lofreq, hifreq, filtfft

# this function asks for and applies an ACAUSAL 4th order bandpass Butterworth filter
def acausal_butter_filter(inst_ty, sps, data, seltask):
    import numpy as np
    from scipy.signal import butter, filtfilt, iirfilter , lfilter

    # filter order
    ford = 4.0

    # set initial hi & lo pass frequencies
    if inst_ty == 'B' or inst_ty == 'N' or inst_ty == 'H':
        lofreq = 0.02
    else:
        lofreq = 0.2

    # override if plotting FFT
    if seltask == '2':
        lofreq = 0.02

    # set low pass at Nyquist frequency
    hifreq = sps / 2.0

    # ask for user input
    txt = '\n'+'High pass Butterworth corner ['+str(lofreq)+' Hz] > '
    var = raw_input(txt)
    if var != '':
        lofreq = float(var)

    txt = '\n'+'Low pass Butterworth corner ['+str(hifreq)+' Hz] > '
    var = raw_input(txt)
    if var != '':
        hifreq = float(var)

    # do high pass filter
    angc = 2.0 * np.pi
    fe = 0.5 * sps
    f1 = lofreq * angc / fe
    f2 = hifreq * angc / fe
    f = np.array([f1, f2])
#    f = [f1, f2]
#    b, a = butter(ford, f, btype = 'high')
    [b, a] = iirfilter(ford, f1, btype='high', ftype='butter', output='ba')
    firstpass = lfilter(b, a, data[0])
    filtwave = lfilter(b, a, firstpass[::-1])[::-1]
#    filtwave = filtfilt(b, a, data)

#    [b, a] = iirfilter(ford, f2, btype='lowpass', ftype='butter', output='ba')
#    firstpass = lfilter(b, a, filtwave)
#    filtwave = lfilter(b, a, firstpass[::-1])[::-1]
#    filtwave = filtfilt(b, a, data)

#    # do low pass filter
#    f = hifreq / fe
#    b, a = butter(ford, f, btype = 'low')
#    filtwave = filtfilt(b, a, filtwave, padlen=150)

    return lofreq, hifreq, filtwave

# This fuction gets the FFT of the data stream between start
# index x1 and stop index x2
def calc_fft(data, sps):
#    import matplotlib.pyplot as plt
    import numpy as np
    try:
        n = len(data)
    
        # calc FFT
        wavfft = np.fft.fft(data, n)
        
        # get frequencies
        freq = np.fft.fftfreq(n, d=1./sps)
        
    except:
        n = len(data[0])
        
        # calc FFT
        wavfft = np.fft.fft(data[0],n)
        
        freq = freq.reshape(1,n)
        
    return freq, wavfft

# return instrument corrected velocity
def get_cor_velocity(corfftr, corffti, freq, inst_ty):
    import numpy as np
    
    # get velocity time history for PGV
    complex_array = corfftr + 1j*corffti
    
    if inst_ty == 'N': # if accelerometer
        complex_array[1:] = complex_array[1:] / (2 * np.pi * abs(freq[1:]))
        complex_array[0] = 0 + 1j*0
        
    n = len(corfftr)
    ivel = np.fft.ifft(complex_array,n)
    pgv = max(abs(ivel.real))
    
    return pgv, ivel


# this function prepares data for the response spectral calculation
def prep_psa(corfftr, corffti, freq, inst_ty):
    import numpy as np
    
    # get velocity time history for PGV
    #complex_array = corfftr[0] + 1j*corffti[0]
    complex_array = corfftr + 1j*corffti
    
    if inst_ty == 'N': # if accelerometer
        complex_array[1:] = complex_array[1:] / (2 * np.pi * abs(freq[1:]))
        complex_array[0] = 0 + 1j*0
        
    n = len(corfftr)
    ivel = np.fft.ifft(complex_array,n)
    pgv = max(abs(ivel.real))
        
    # get acceleration time history
    complex_array = corfftr + 1j*corffti
    
    if inst_ty != 'N': # if seismometer
        complex_array = complex_array * (2 * np.pi * abs(freq))

    n = len(corfftr)
    iacc = np.fft.ifft(complex_array,n)

    return pgv, iacc, ivel
    
# this function prepares data for the response spectral calculation
def prep_psa_simple(corfftr, corffti, freq, inst_ty):
    import numpy as np
    
    # get acceleration time history
    complex_array = corfftr + 1j*corffti
    
    if inst_ty != 'N': # if seismometer
        complex_array = complex_array * (2 * np.pi * abs(freq))

    n = len(corfftr)
    iacc = np.fft.ifft(complex_array,n)

    return iacc

    
# this function calculated response spectral acceleration of a SODF
# structure with a damping of h (usually 5%)
def calc_response_spectra(iacc, sps, h, minT, maxT):
    import numpy as np

    print('\nCalculating ' + str(h) + '% damped response spectrum. Please be patient.')

    damp = h * 0.01
    dt = 1. / sps

    # set size of 100 periods
    nT = 100    
    T = np.logspace(-2,1,num=nT)
    
    # strip periods that are out of range
    idx = np.where((T >=  minT) & (T <= maxT))[0]
    nTstrip = len(idx)
    Tstrip = T[idx]
    
    psa = np.zeros((nTstrip,1))
    sd = np.zeros((nTstrip,1))

    # set PGA value
    pga = np.max(np.abs(iacc))

    # now start first loop
    for i in range(0, nTstrip):
        K = 2 * np.pi / Tstrip[i]
        C = 2 * damp * K
        beta = 0.25
        gamma = 0.5
        K *= K

        B = 1.0 / (beta*dt*dt) + (gamma*C)/(beta*dt)
        A = B + K
        E = 1.0 / (beta*dt) + (gamma/beta-1)*C
        G = 1.0 / (2*beta) - 1.0

        x = 0;
        xp = 0;
        xpp = iacc[0]
        maxx = x

        # start second loop
        for j in range(0,len(iacc)):
            xn = (-iacc[j]+B*x+E*xp+G*xpp) / A
            xppn = (xn-x-dt*xp-dt*dt*xpp / 2+dt*dt*beta*xpp) / (beta*dt*dt)
            xpn = xp+dt*xpp+dt*gamma*(xppn-xpp)

            x = xn
            xpp = xppn
            xp = xpn
            xn = abs(x)

            # save max val
            if xn > maxx:
                maxx = xn

        sd[i] = maxx
        psa[i] = maxx*K

    return Tstrip, psa.flatten(), pga