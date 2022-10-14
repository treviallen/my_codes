"""
These are tools for plotting

"""

# a function to annotate max and min ground motions
def annotate_maxmin(plt, ax ,data):
    import numpy as np

    maxval = str("%0.3f" % max(data))
    minval = str("%0.3f" % min(data))
    ylim = ax.get_ylim()
    xlim = ax.get_xlim()
    ymin = -1*max(np.abs(ylim))
    ymax = max(np.abs(ylim))
    ax.set_ylim((ymin, ymax))
    plt.annotate('max = '+maxval, (0,0), fontsize=10, \
                 xytext=((xlim[1]-xlim[0])*0.02, ymax*0.80), xycoords='data')
    plt.annotate('min = '+minval, (0,0), fontsize=10, \
                 xytext=((xlim[1]-xlim[0])*0.02, ymin*0.85), xycoords='data')

# this function allows interactive trimming of time histories if required
def trim_wave(data, sps, inst_ty, fftanal):
    import numpy as np
    import matplotlib.pyplot as plt
    plt.ion()

    figure = plt.figure(6,figsize=(16,9))
    # do smart way
    #try:
    nsamp = len(data)
    
    # set time array
    tvect = np.linspace(0, nsamp/sps, num=nsamp)
    
    ax = figure.add_subplot(3,1,1)
    ax.plot(tvect, data,'-', lw=0.5,color='orange')
    '''
    # do stupid way!
    except:
        nsamp = len(data)
    
        # set time array
        tvect = np.reshape(np.linspace(0,nsamp/sps,num=nsamp),(1,nsamp))
    
        ax = figure.add_subplot(3,1,1)
        ax.plot(tvect[0],data[0],'-', lw=0.5, color='orange')
    '''
    if fftanal == False:
        plt.title('Select by clicking twice about time window of interest' + '\n' \
                   + '(hint: set t2 < t1 to return whole wave)')
    else:
        plt.title('Reselect wave for frequency domain analysis' + '\n' \
                   + '(hint: set t2 < t1 to return whole wave)')

    if inst_ty == 'N':
        plt.ylabel('Acceleration (Counts)')
    else:
        plt.ylabel('Velocity (Counts)')

    plt.xlabel('Time (s)')
    
    try:
        annotate_maxmin(plt, ax ,data)
    except:
        annotate_maxmin(plt, ax ,data[0])
    pts = np.asarray(plt.ginput(2,timeout=-1))

    # get x indicies to replot
    x1 = int(pts[0,0]*sps)
    x2 = int(pts[1,0]*sps)

    # if want whole window, make x2 < x1
    if x2 <= x1:
        x1 = 0
        x2 = nsamp - 1

    ax = figure.add_subplot(3,1,3)
    plt.title('Trimmed Time History')
    try:
        ax.plot(tvect[x1:x2],data[x1:x2],'-', lw=0.5,color='green')
        annotate_maxmin(plt, ax ,data[x1:x2])
    except:
        ax.plot(tvect[0,x1:x2],data[0,x1:x2],'-', lw=0.5,color='green')
        annotate_maxmin(plt, ax ,data[0,x1:x2])
    
    plt.ylabel('Velocity (Counts)')
    plt.xlabel('Time (s)')
    
#    figure.savefig('trim_example.png', format='png')
    plt.close(figure)
    #plt.show()

    return x1, x2

# plot displacement, velocity and acceleration time history
# assumes velocity input!
def plot_dva(freq, sps, corfftr, corffti, header, inst_ty, chan_no):
    import numpy as np
    import matplotlib.pyplot as plt
    import os
    '''
    n = len(corfftr[0])
    tvect = np.reshape(np.linspace(0,n/sps,num=n),(1,n))
    '''
    n = len(corfftr)
    #print(n
    tvect = np.linspace(0,n/sps,num=n)

    # if velocity seismometer
    if inst_ty == 'S' or inst_ty == 'B' or inst_ty == 'H' or inst_ty == 'E':

        # get displacement wave
        '''
        freq[0,0] = 1.0
        dispfftr = corfftr[0] / (2 * np.pi * abs(freq[0]))
        dispffti = corffti[0] / (2 * np.pi * abs(freq[0]))
        dispfftr[0] = 0
        dispffti[0] = 0
        freq[0,0] = 0
        '''
        freq[0] = 1.0
        dispfftr = corfftr / (2 * np.pi * abs(freq))
        dispffti = corffti / (2 * np.pi * abs(freq))
        dispfftr[0] = 0
        dispffti[0] = 0
        freq[0] = 0
        
        complex_array = dispfftr + 1j*dispffti

        idisp = np.fft.ifft(complex_array,n)

        # get velocity wave
        '''
        velfftr = corfftr[0]
        velffti = corffti[0]
        '''
        velfftr = corfftr
        velffti = corffti
        complex_array = velfftr + 1j*velffti
        ivel = np.fft.ifft(complex_array,n)

        # get acceleration wave
        '''
        accfftr = corfftr[0] * (2 * np.pi * abs(freq[0]))
        accffti = corffti[0] * (2 * np.pi * abs(freq[0]))
        '''
        accfftr = corfftr * (2 * np.pi * abs(freq))
        accffti = corffti * (2 * np.pi * abs(freq))
        complex_array = accfftr + 1j*accffti
        iacc = np.fft.ifft(complex_array,n)

    # if accelerometer
    elif inst_ty == 'N':

        # get displacement wave
        freq[0] = 1.0
        dispfftr = corfftr / (2 * np.pi * abs(freq))**2
        dispffti = corffti / (2 * np.pi * abs(freq))**2
        dispfftr[0] = 0.
        dispffti[0] = 0.

        complex_array = dispfftr + 1j*dispffti

        idisp = np.fft.ifft(complex_array,n)

        # get velocity wave
        velfftr = corfftr / (2 * np.pi * abs(freq))
        velffti = corffti / (2 * np.pi * abs(freq))
        freq[0] = 0.
        velfftr[0] = 0.
        velffti[0] = 0.
        complex_array = velfftr + 1j*velffti
        ivel = np.fft.ifft(complex_array,n)

        # get acceleration wave
        accfftr = corfftr
        accffti = corffti
        complex_array = accfftr + 1j*accffti
        iacc = np.fft.ifft(complex_array,n)

    # plot waves
    fignum = 10 + chan_no
    figure = plt.figure(fignum,figsize=(16,9))
    ax = figure.add_subplot(3,1,1)
    #ax.plot(tvect[0],idisp.real * 1000,'-',color='red') # converted to mm
    ax.plot(tvect,idisp.real * 1000,'-', lw=0.5,color='red') # converted to mm
    plt.ylabel('Displacement (mm)')
    plt.title(header)
    annotate_maxmin(plt, ax ,idisp.real * 1000)

    ax = figure.add_subplot(3,1,2)
    ax.plot(tvect,ivel.real*1000,'-', lw=0.5,color='green') # converted to mm/s
    #ax.plot(tvect[0],ivel.real*1000,'-',color='green') # converted to mm/s
    plt.ylabel('Velocity (mm/s)')
    annotate_maxmin(plt, ax ,ivel.real * 1000)

    ax = figure.add_subplot(3,1,3)
    #ax.plot(tvect[0],iacc.real * 100 / 9.81,'-',color='blue') # converted to %g
    ax.plot(tvect,iacc.real * 100 / 9.81,'-', lw=0.5,color='blue') # converted to %g
#    plt.ylabel(r'Acceleration (m/s$^2$)')
    plt.ylabel('Acceleration (%g)')
    plt.xlabel('Time (s)')
    annotate_maxmin(plt, ax ,iacc.real * 100 / 9.81)
    plt.show()

    filename = header + '.png'
    outdir = 'fig'
    outfile = os.path.join(outdir,filename)
    # try saving to fig directory
    try:
        dirList=os.listdir(outdir)
    except:
        os.mkdir('fig')

    figure.savefig(outfile, format='png',bbox_inches='tight')

    return idisp.real, ivel.real, iacc.real # only return real component

def plot_WoodAnderson(wadisp, sps, header, chan_no):
    import numpy as np
    import matplotlib.pyplot as plt
    plt.ion()

    n = len(wadisp)
    tvect = np.reshape(np.linspace(0,n/sps,num=n),(1,n))

    # plot
    fignum = 40 + chan_no
    figure = plt.figure(fignum,figsize=(16,9))
    ax = figure.add_subplot(3,1,1)
    ax.plot(tvect[0],wadisp,'-', lw=0.5,color='orange')
    plt.ylabel('WA Displacement (mm)')
    plt.xlabel('Time (s)')
    plt.title(header)
    annotate_maxmin(plt, ax ,wadisp)
    pts = np.asarray(plt.ginput(2,timeout=-1))

    # get x indicies to replot
    x1 = int(pts[0,0]*sps)
    x2 = int(pts[1,0]*sps)

    # if want whole window, make x2 < x1
    if x2 <= x1:
        x1 = 0
        x2 = n - 1

    # replot trimmed wave
    #print(np.shape(wadisp))
    ax = figure.add_subplot(3,1,3)
    ax.plot(tvect[0,x1:x2],wadisp[x1:x2],'-', lw=0.5,color='green')
    plt.ylabel('WA Displacement (mm)')
    plt.xlabel('Time (s)')
    plt.title(header)
    annotate_maxmin(plt, ax ,wadisp)
    #plt.savefig('wa_example.png', format='png', bbox_inches='tight')

    plt.close(figure)
    #plt.show()

    return wadisp[x1:x2], x1, x2

# this function plots the 5% damped response spectra
def plot_response_spectra(T, psa, pga, header, chan_no):
#    import numpy as np
    import matplotlib.pyplot as plt

    # plot response
    fignum = 30 + chan_no
    fig = plt.figure(fignum,figsize=(9,9))
    ax1 = fig.add_subplot(111)
    ax1.loglog(T, psa,'b')
    ax1.set_xlabel('Period (s)')
    ax1.set_ylabel('Pseudo Spectral Acceleration (m/s^2)')
    plt.savefig('psa_example.png', format='png')
    plt.title(header)

    plt.show()

# this function plots displacement Fourier sprectrum
def plot_fft(sps, freq, corfftr, corffti, header, inst_ty):
    import numpy as np
    import matplotlib.pyplot as plt
    
    if inst_ty == 'N': # for accelerometer
        corfftr = corfftr[1:]/(2.*(np.pi)*freq[1:])**2
        corffti = corffti[1:]/(2.*(np.pi)*freq[1:])**2
    else: # or seismometer
        #corfftr = corfftr[0,1:]/(2.*(np.pi)*freq[0,1:])
        #corffti = corffti[0,1:]/(2.*(np.pi)*freq[0,1:])
        corfftr = corfftr[1:]/(2.*(np.pi)*freq[1:])
        corffti = corffti[1:]/(2.*(np.pi)*freq[1:])
    dispamp = np.sqrt((1./sps)**2 * (corfftr**2 + corffti**2))

    # get plotting length
    pl = int(round(len(dispamp) * 0.45))
    #pltfreq = freq[0,2:]
    pltfreq = freq[2:]
    pltdisp = dispamp[1:]
    
    plt.figure(2,figsize=(9,9))
    plt.loglog(pltfreq[0:pl], pltdisp[0:pl],'b')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Displacement Fourier Spectrum (m-s)')
    plt.savefig('fft_example.png', format='png')
    plt.show()

# this function plots the instrument response
def plot_instrument_resp(sta, inst_ty, freq, nat_freq, damping, sen, recsen, \
                         gain, pazfile, chan_no):
    import numpy as np
    import matplotlib.pyplot as plt
    import response

    # get instrument response
    if pazfile == 'NULL':
        resp_real, resp_imag = response.fap_response(freq, nat_freq, damping, sen, \
                                                     recsen, gain, inst_ty)
    else:
        resp_real, resp_imag = response.paz_response(freq, pazfile, sen, recsen, \
                                                     gain, inst_ty)

    respa = np.sqrt(resp_real**2 + resp_imag**2)
#    respp = np.arctan2(resp_real,resp_imag) * 180.0 / np.pi
    respp = np.angle(resp_real + 1j*resp_imag, deg=True)

    # plot amplitude
    n = int(len(freq) / 2)
    fignum = 80 + chan_no
    fig = plt.figure(fignum,figsize=(16,7))
    ax1 = fig.add_subplot(111)
    p1 = ax1.loglog(freq[1:n], respa[1:n],'b')
    ax1.set_xlabel('Frequency (Hz)')
    if inst_ty == 'N':
        ax1.set_ylabel('Sensor Magnification (Counts/m/s**2)')
    else:
        ax1.set_ylabel('Sensor Magnification (Counts/m/s)')

    # plot phase
    ax2 = ax1.twinx()
    p2 = ax2.semilogx(freq[1:n], respp[1:n],'g')
    ax2.set_ylabel('Phase Angle (Degrees)')
    plt.legend( (p1[0], p2[0]), ['Amplitude','Phase'],loc=2)
    if pazfile == 'NULL':
        plt.title(sta)
    else:
        plt.title(sta + ' (' + pazfile.upper().strip('.paz') + ')')
#    fig.savefig('resp_example.png', format='png')
    plt.show()



