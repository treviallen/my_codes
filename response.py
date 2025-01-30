# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 11:01:37 2012

@author: u56903

This script checks to see if instrument FAP information exists for the current
wave file

"""

def get_response_info(sta,recdate,chan,netid):
    import datetime as dt
    from os import getcwd
    
    nat_freq = -12345
    inst_ty = ''
    damping = 0.707
    sen = 1.0
    recsen = 1.0
    gain = 1.0
    #netid = ''
    pazfile = 'NULL' 
    
    # make hack to network codes
    if netid.strip()=='ME' or netid.strip() == 'SR' or netid.strip() == 'VI' or netid.strip() == 'GM':
        netid = 'MEL'
    elif netid.strip() == 'AB':
        if sta == 'CCRM' or sta == 'GLEAD' or sta == 'MALDN' or sta == 'CAWDR' or sta == 'CORA' \
           or sta == 'FINNS' or sta == 'BOGA' or sta == 'MILD' or sta == 'BUCHN' or sta == 'GOOG':
            netid = 'MEL'
    elif netid.strip() == 'OZ':
        if sta == 'WILL' or sta == 'GLMM' or sta == 'JEER' or sta == 'DTMM' or sta == 'HOPM' \
           or sta == 'KORUM':
            netid = 'MEL'
    
    if sta == 'DROM' or sta == 'KORUM' or sta == 'JENM' or sta == 'EASTN' or sta == 'GLADS' or sta == 'TOMM' \
       or sta == 'CDNM' or sta == 'FSHM' or sta == 'AWOON' or sta == 'TPND' or sta == 'GLADS' or sta == 'TOMM':
        netid = 'MEL'   
        
    if netid == 'MEL':
        netid = 'OZ'
        
    if netid == '':
        if sta == 'MALDN' or sta == 'GLEAD':
            netid = 'OZ'
            
        if sta == 'FR27' or sta == 'NAP':
            netid = 'AD'
    
    '''
    elif netid.strip() == 'AB':
        netid = 'UM'
    '''    
    if sta == 'RIV' or sta == 'KIM' or sta == 'CHI' or sta == 'NLD' or sta == 'NPS' or sta == 'DOW' \
       or sta.startswith('NOR'):
        netid = 'AU'
        
    '''
    if sta == 'ARKL' or sta == 'HWK':
        netid = 'AD'
    '''    
    if sta == 'S88U' or sta == 'STGU' or sta == 'HODL':
        netid = 'UM'
    
    # check if sitename in file
    cwd = getcwd()
    if cwd.startswith('/nas'):
        stalist = '/nas/users/u56903/unix/Code/my_codes/stationlist.dat'        
    elif cwd.startswith('C:'):
        stalist = 'C:\Code\my_codes\stationlist.dat' 
    else:
        stalist = '//Users//trev//Documents//Code//my_codes//stationlist.dat'
        
    stlo = -12345.0
    stla = -12345.0
    try:
        stadat = open(stalist).readlines()
        
    # reset stalist to local drives on rhe-compute
    except:
        #stalist = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Ground_Motion/Data/stationlist.dat'   
        stadat = open(stalist).readlines()
        
    for line in stadat[1:]:
        # skip first row
        tmp = line.rstrip().split('\t')
        #print(line)
        mindate = dt.datetime.strptime(tmp[2], "%Y%m%d")
        maxdate = dt.datetime.strptime(tmp[3], "%Y%m%d")
        #print(mindate, maxdate, recdate)
        if tmp[0] == sta and tmp[12].strip() == chan and tmp[6] == netid:
            if recdate >= mindate and recdate <= maxdate:
                inst_ty = tmp[1]
                stlo = float(tmp[4])
                stla = float(tmp[5])
                #netid = tmp[6]
                nat_freq = float(tmp[7])
                damping = float(tmp[8])
                sen = float(tmp[9])
                recsen = float(tmp[10])
                gain = float(tmp[11])
                pazfile = tmp[13].rstrip()
                
        elif netid == 'AB':
            if tmp[0] == sta and tmp[12].strip() == chan:
                if recdate >= mindate and recdate <= maxdate:
                    
                    inst_ty = tmp[1]
                    stlo = float(tmp[4])
                    stla = float(tmp[5])
                    netid = tmp[6]
                    nat_freq = float(tmp[7])
                    damping = float(tmp[8])
                    sen = float(tmp[9])
                    recsen = float(tmp[10])
                    gain = float(tmp[11])
                    pazfile = tmp[13].rstrip()
                    
    if stlo == -12345:
        print(recdate,': Station', sta, chan, 'not found...')

    return nat_freq, inst_ty, damping, sen, recsen, gain, pazfile, stlo, stla, netid

# parses staionlist and outputs to list of dictionaries
def stationlist2dict():
    import datetime as dt
    from os import getcwd
    
    # check if sitename in file
    cwd = getcwd()
    if cwd.startswith('/nas'):
        stalist = '/nas/users/u56903/unix/Code/my_codes/stationlist.dat'
    else:
        stalist = '//Users//trev//Documents//Code//my_codes//stationlist.dat'
        
    stadat = open(stalist).readlines()
    
    stalistDict = []
    for line in stadat[1:]:
        if not line.startswith('#'):
            tmp = line.rstrip().split('\t')
            #print(line)
            mindate = dt.datetime.strptime(tmp[2], "%Y%m%d")
            maxdate = dt.datetime.strptime(tmp[3], "%Y%m%d")
            #print(mindate, maxdate, recdate)
            sta = tmp[0]
            inst_ty = tmp[1]
            stlo = float(tmp[4])
            stla = float(tmp[5])
            netid = tmp[6]
            nat_freq = float(tmp[7])
            damping = float(tmp[8])
            sen = float(tmp[9])
            recsen = float(tmp[10])
            gain = float(tmp[11])
            chan = tmp[12]
            pazfile = tmp[13].rstrip()
            
            sld = {'sta':sta, 'start':mindate, 'stop':maxdate, 'insttype':inst_ty, 'stlo':stlo, 'stla':stla, 'channel':chan, \
                   'netid':netid, 'sen_sensitivity':sen, 'rec_sensitivity':recsen, 'gain':gain, 'pazfile':pazfile}
                 	
            stalistDict.append(sld)
        
    return stalistDict

def iris_gmap2stationlist(gmap):

	lines = open(gmap).readlines()
	stas = []
	if lines[0].startswith('#Network'):
		skiprows = 1
	else:
		skiprows = 3 
	for line in lines[skiprows:]:
		line = line.replace('||', '|')
		dat = line.replace('"','').strip().split('|')
		#print(dat)
		tdict = {'net': dat[0], 'sta': dat[1].upper(), 'lat':float(dat[2]), 'lon':float(dat[3]), \
			       'start':dat[-2][0:10].replace('-',''), 'stop': dat[-1][0:10].replace('-','')}
		stas.append(tdict)
		
	# now write txt
	txt = ''
	txt2 = ''
	for sta in stas:
		txt += '\t'.join((sta['sta'],'B',sta['start'],sta['stop'],str('%0.4f' % sta['lon']),str('%0.4f' % sta['lat']),sta['net'],'-12345','-12345','-12345','-12345','1','HHE','inst.paz')) + '\n'
		txt += '\t'.join((sta['sta'],'B',sta['start'],sta['stop'],str('%0.4f' % sta['lon']),str('%0.4f' % sta['lat']),sta['net'],'-12345','-12345','-12345','-12345','1','HHN','inst.paz')) + '\n'
		txt += '\t'.join((sta['sta'],'B',sta['start'],sta['stop'],str('%0.4f' % sta['lon']),str('%0.4f' % sta['lat']),sta['net'],'-12345','-12345','-12345','-12345','1','HHZ','inst.paz')) + '\n'
	
		txt2 += '\t'.join((sta['sta'],str(sta['lon']),str(sta['lat']),'1',sta['start'][0:4],sta['start'][4:6],sta['stop'][0:4],sta['stop'][4:6])) + '\n'
	# now write
	f = open(gmap[0:-4]+'_stationlist.txt', 'w')
	f.write(txt)
	f.close()

	f = open(gmap[0:-4]+'_station_data.dat', 'w')
	f.write(txt2)
	f.close()
	
def fdsn_text2stationlist(gmap):

	lines = open(gmap).readlines()
	stas = []
	if lines[0].startswith('#Network'):
		skiprows = 1
	else:
		skiprows = 3 
	for line in lines[skiprows:]:
		line = line.replace('||', '|')
		dat = line.replace('"','').strip().split('|')
		#print(dat)
		tdict = {'net': dat[0], 'sta': dat[1].upper(), 'chan': dat[2].upper(), 'lat':float(dat[3]), 'lon':float(dat[4]), \
			       'sensor':dat[9], 'sensitivity':float(dat[10]), \
			       'start':dat[-2][0:10].replace('-',''), 'stop': dat[-1][0:10].replace('-','')}
		stas.append(tdict)
		
	# now write txt
	txt = ''
	for sta in stas:
		
		txt += '\t'.join((sta['sta'], sta['chan'][0], sta['start'],sta['stop'],str('%0.4f' % sta['lon']), \
		                  str('%0.4f' % sta['lat']),sta['net'],'-12345','-12345',str('%0.6e' % sta['sensitivity']),'1.0', \
		                  '1.0',sta['chan'],sta['sensor'])) + '\n'
	# now write
	f = open(gmap[0:-4]+'_stationlist.txt', 'w')
	f.write(txt)
	f.close()


# extracts data for stationlist from eqwave file
def eqwave2staionlist(wavfile):
    from readwaves import readeqwave

    sta, comps, recdate, sec, sps, data, nsamp, cntpvolt, sen, gain, datadict = readeqwave(wavfile)
    
    lines = ''
    for i in range(0, len(sta)):
        line = '\t'.join((sta[i], comps[i][0], recdate[i][0:8], '25990101', \
        	                datadict['lons'][i], datadict['lats'][i], 'MEL', '-12345', '-12345', \
                          sen[i], cntpvolt[i], gain[i], comps[i], 'NULL'))
        lines += line + '\n'
    
    # make filename
    slfile = wavfile.split('.')[0].split('/')[1]+'.stationlist'
    slfile = slfile.replace(' ', '_')
    f = open(slfile, 'wb')
    f.write(lines)
    f.close()

# this function lists available PAZ files for selection
def get_paz_list():
    import os

    # read paz file
    dirList=os.listdir('paz')
    i = 0
    for fname in dirList:
        #print(str(i+1) + ')\t' + fname.strip('.paz'))
        i += 1

    try:
        var = raw_input('\nSelect PAZ file > ')
    except:
        var = input('\nSelect PAZ file > ')
    pazfile = dirList[int(var)-1]

    return pazfile

# this function reads PAZ files
def read_pazfile(in_pazfile):
    from os import path
    import numpy as np
    from os import getcwd
    
    cwd = getcwd()
    # open paz file
    if cwd.startswith('/nas'):
        pazpath = '//nas//users//u56903//unix//paz' # for rhe-compute
    if cwd.startswith('C:'):
        pazpath = 'C:\\Code\\paz' # for rhe-compute        
    else:
        pazpath = '//Users//trev//Documents//Earthquake_Data//paz' # for bob
    
    try:
        pazfile = path.join(pazpath, in_pazfile)
        #print(pazfile)
        paztxt = open(pazfile).readlines()
    except:
        pazpath = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Ground_Motion/Data/paz'
        pazfile = path.join(pazpath, in_pazfile)
        paztxt = open(pazfile).readlines()
        
    # get zeros first
    zeros = []
    num = paztxt[0].split('\t')
    nzeros = int(num[1])
    for i in range(1,nzeros+1):
        tmpz = paztxt[i].strip('\n').split('\t')
#        print(tmpz
        zeros.append(complex(float(tmpz[0]),float(tmpz[1]))*2*np.pi) # convert to angular frequency
#    print(zeros

    # now get poles
    poles = []
    num = paztxt[nzeros+1].split('\t')
    npoles = int(num[1])
    for i in range(nzeros+2,nzeros+npoles+2):
        tmpp = paztxt[i].strip('\n').split('\t')
#        print(tmpp
        poles.append(complex(float(tmpp[0]),float(tmpp[1]))*2*np.pi) # convert to angular frequency
#    print(poles

    # get constant
    constant = paztxt[nzeros+npoles+2].split('\t')
    constant = float(constant[1])

    # get normalising frequency
    normf = paztxt[nzeros+npoles+3].split('\t')
    normf = float(normf[1])

    return poles, zeros, constant, normf

# this function reads PAZ files
def read_sac_respfile(in_pazfile):
    from os import path
    import numpy as np
    from os import getcwd
    
    paztxt = open(in_pazfile).readlines()
        
    # get zeros first
    zeros = []
    num = paztxt[0].split()
    nzeros = int(num[1])
    for i in range(1,nzeros+1):
        tmpz = paztxt[i].strip('\n').split()
        zeros.append(complex(float(tmpz[0]),float(tmpz[1]))) #*2*np.pi) # convert to angular frequency

    # now get poles
    poles = []
    num = paztxt[nzeros+1].split()
    npoles = int(num[1])
    for i in range(nzeros+2,nzeros+npoles+2):
        tmpp = paztxt[i].strip('\n').split()
        poles.append(complex(float(tmpp[0]),float(tmpp[1]))) #*2*np.pi) # convert to angular frequency

    # get constant
    constant = paztxt[nzeros+npoles+2].split()
    constant = float(constant[1])

#    # get normalising frequency
#    normf = paztxt[nzeros+npoles+3].split('\t')
#    normf = float(normf[1])

    return poles, zeros, constant

# returns normalisation factor from lists of poles & zeros
def normfact_from_paz(p, z):
    '''
    p = list of poles in rad-sec in real+imagj fmt, eg: [-0.01234+0.01234j, -0.01234-0.01234j, -39.1800+49.1200j, -39.1800-49.1200j]
    z = list of zeros in rad-sec, eg: [0.,0.,0.]
    '''
    
    import scipy.signal as signal
    
    #print(signal.ltisys.zpk2tf(z, p, 1.))
    return signal.ltisys.zpk2tf(z, p, 1.)[1][2]

# this function leads the user to input station information
def enter_response_info(sta, chan, sps):
    print('\nEnter site information for ' + sta + ' ' + chan)
    try:
        var = raw_input('\n'+'Velocity sensor or accelerometer ([v]/a)? > ')
    except:
        var = input('\n'+'Velocity sensor or accelerometer ([v]/a)? > ')
     
    if var == '' or var == 'v':
        inst_ty = 'V'
    elif var == 'a':
        inst_ty = 'N' # use IRIS code

    var = raw_input('\n'+'Open date of station (YYYYMMDD [00010101]) > ')
    if var == '':
        mindate = '00010101'
    else:
        mindate = var

    var = raw_input('\n'+'Close date of station (YYYYMMDD [25990101]) > ')
    if var == '':
        maxdate = '25990101'
    else:
        maxdate = var

    var = raw_input('\n'+sta+' network ID [PV] > ')
    if var == '':
        netid = 'PV'
    else:
        netid = var

    print('\nLook-up station details here: http://www.isc.ac.uk/registries/listing/')

    var = raw_input('\n'+sta+' station latitude (decimal degrees) > ')
    if var == '':
        stla = -12345
    else:
        stla = float(var)

    var = raw_input('\n'+sta+' station longitude (decimal degrees) > ')
    if var == '':
        stlo = -12345
    else:
        stlo = float(var)

    # ask for pole-zero file
    var = raw_input('\nUse poles and zeros file ([y]/n)?  > ')
    if var == 'y' or var == '':
        # get paz file list
        pazfile = get_paz_list()

        # enter null values for nat freq and damp
        nat_freq = -12345
        damping = -12345

    # ask for user input
    else:
        var = raw_input('\n'+'Enter sensor natural frequency [1.0 Hz] > ')
        if var == '':
            nat_freq = 1.0
        else:
            nat_freq = float(var)

        # get IRIS instrument type
        if 1.0 / nat_freq < 10 and inst_ty == 'V': #sec
            if sps >= 80:
                inst_ty = 'E'
            else:
                inst_ty = 'S'
        elif inst_ty == 'V':
            if sps >= 80:
                inst_ty = 'H'
            else:
                inst_ty = 'B'

        var = raw_input('\n'+'Enter seismometer damping [0.707] > ')
        if var == '':
            damping = 0.707
        else:
            damping = float(var)

        # set pazfile to null value
        pazfile = 'NULL'

    # if seismometer type not automatically assigned above, ask user
    if inst_ty == 'V':
        var = raw_input('\nShort Period or Broad Band instrument ([s]/b)? > ')
        if var == '' or var == 's':
            if sps >= 80:
                inst_ty = 'E'
            else:
                inst_ty = 'S'
        elif var == 'b':
            if sps >= 80:
                inst_ty = 'H'
            else:
                inst_ty = 'B'

    # if using seismometer
    if inst_ty == 'S' or inst_ty == 'H' or inst_ty == 'B' or inst_ty == 'E':
        var = raw_input('\n'+'Enter seismometer sensitivity [1.0 V/m/s] > ')
        if var == '':
            sen = 1.0
        else:
            sen = float(var)

    # if using accelerometer
    else:
        var = raw_input('\n'+'Enter accelerometer sensitivity [1.0 V/m/s**2] > ')
        if var == '':
            sen = 1.0
        else:
            sen = float(var)

    var = raw_input('\n'+'Enter recorder sensitivity [1.0 Count/V] > ')
    if var == '':
        recsen = 1.0
    else:
        recsen = float(var)

    var = raw_input('\n'+'Enter recorder gain [1.0] > ')
    if var == '':
        gain = 1.0
    else:
        gain = float(var)

    return inst_ty, mindate, maxdate, stlo, stla, netid, nat_freq, \
        damping, sen, recsen, gain, pazfile

# write new response info to file
def write_response_info(sta, inst_ty, mindate, maxdate, stlo, stla, netid, nat_freq, \
        damping, sen, recsen, gain, chan, pazfile):

    stalist = 'stationlist.dat'
    newtxt = '\t'
    joinstr = (sta, inst_ty, mindate, maxdate, str(stlo), str(stla), netid, \
              str(nat_freq), str(damping), str(sen), str(recsen), str(gain), \
              chan, pazfile)
    newtxt = newtxt.join(joinstr)
    newsta = open(stalist,'a')
    newsta.write('\n'+newtxt)
    newsta.close()

# generate real and imaginary components for FAP response
def fap_response(freq, nat_freq, damping, sen, recsen, gain, inst_ty):
    import numpy as np
    import scipy.signal as signal

    angc = 2.0 * np.pi
    ampfact = sen * recsen * gain

    poles = [-(damping + np.sqrt(1 - damping ** 2) * 1j) * angc * nat_freq]
    poles.append(-(damping - np.sqrt(1 - damping ** 2) * 1j) * angc * nat_freq)
    zeros = [0j, 0j]

    # if accelerometer, convert units from counts/g to counts/m/s**2
#    if inst_ty == 'N':
#        ampfact /= 9.81

    b, a = signal.ltisys.zpk2tf(zeros, poles, ampfact)
    w, resp = signal.freqs(b, a, freq * angc)

    return resp.real, resp.imag

# generate real and imaginary components for PAZ response
def paz_response(freq, pazfile, sen, recsen, gain, inst_ty):
    import numpy as np
    import scipy.signal as signal
    
    dispResp = False
    
    #print('freq', np.shape(freq)
    #print(freq

    # read PAZ file
    #print(pazfile)
    poles, zeros, constant, normf = read_pazfile(pazfile)
    #print(poles, zeros)
    angc = 2.0 * np.pi

    # if constant unknown, get nearest freq to normalising frequency
    
    if constant == -12345:
        # if normalizing frequency also unknown, assume normf = 5.0 Hz
        if normf == -12345:
            normf = 5.0 # maybe change this to 5 Hz to be safe
        elif normf == -99:
            dispResp = True
        
        # get amp at normf - updated 2019-01-04
        b, a = signal.ltisys.zpk2tf(zeros, poles, 1.0)
        w, resp = signal.freqs(b, a, freq * angc)
        
        cal_resp = abs(np.exp(np.interp(np.log(normf), np.log(freq), np.log(abs(resp)))))
        constant = 1. / cal_resp
        

    # use normalizing factor from file
    else:
        print('\nUsing normalisation constant =', constant, '\n')
        
        if normf == -99.: # displacement response
            #print(constant, sen
            cal_freq = 1.0 # assumed
            #constant *= angc #* cal_freq # convert from disp to vel?
            dispResp = True
            sen /= angc
        else:
            constant *= angc**(len(poles)-len(zeros))
        
    #print(constant, sen, dispResp
    if dispResp == True:
        sen /= angc
    #print(constant, sen
    # combine amp factors
    #print('Norm Const:', constant
    ampfact = sen * recsen * gain * constant

    # if accelerometer, convert units from counts/g to counts/m/s**2
#    if inst_ty == 'N':  # now defining in counts/m/s**2
#        ampfact /= 9.81

    # get shape of response
    b, a = signal.ltisys.zpk2tf(zeros, poles, ampfact)
    w, resp = signal.freqs(b, a, freq * angc)
    
    # check if disp response and convert to velocity
    '''
    print('disptrue', dispResp, normf
    if dispResp == True:
        print('disptrue'
        resp.real *= angc*freq
        #resp.imag *= angc*freq
    '''

    return resp.real, resp.imag

# remove FAP response from recorded FFT
def deconvolve_instrument(rresp, iresp, wavfft):
    idenom = rresp**2 + iresp**2
    idenom[0] = 0.0001

    corfftr = (wavfft.real * rresp + wavfft.imag * iresp) / idenom
    corfftr[0] = 0.
    corffti = (wavfft.imag * rresp - wavfft.real * iresp) / idenom
    corffti[0] = 0.

    return corfftr, corffti
    
def common_resp(freq, nat_freq, damping, sen, recsen, gain, wavfft, inst_ty, pazfile):

    # calculate FAP instrument response
    if pazfile == 'NULL':
        real_resp, imag_resp = fap_response(freq, nat_freq, damping, \
                                                     sen, recsen, gain, inst_ty)

    # or use PAZ file
    else:
        real_resp, imag_resp = paz_response(freq, pazfile, sen, recsen, \
                                                     gain, inst_ty)

    # deconvolve instrument response from record and return corrected FFT
    corfftr, corffti = deconvolve_instrument(real_resp, imag_resp, wavfft)

    return corfftr, corffti, real_resp, imag_resp


# get Wood-Andersion instrument response and convolve with input seismogram
def convolve_WoodAnderson(freq, corfftr, corffti, inst_ty, ampfact=2080, damping=0.8):
    import numpy as np
    import scipy.signal as signal

    # get Wood-Anderson response
    nat_freq = 1.25
    angc = 2.0 * np.pi

    poles = [-(damping + np.sqrt(1 - damping ** 2) * 1j) * angc * nat_freq]
    poles.append(-(damping - np.sqrt(1 - damping ** 2) * 1j) * angc * nat_freq)
    zeros = [0j, 0j]

    b, a = signal.ltisys.zpk2tf(zeros, poles, ampfact)
    w, resp = signal.freqs(b, a, freq * 2 * np.pi)

    # convolve Wood-Anderson response
    wafftr = (corfftr * resp.real - corffti * resp.imag)
    waffti = (corfftr * resp.imag + corffti * resp.real)

    # get Wood-Anderson Displacement time history
    '''
    n = len(wafftr[0])
    freq[0,0] = 1.0
    dispfftr = wafftr[0] / (2 * np.pi * abs(freq[0]))
    dispffti = waffti[0] / (2 * np.pi * abs(freq[0]))
    '''
    
    if inst_ty == 'N':
        n = len(wafftr)
        freq[0] = 1.0
        dispfftr = wafftr / (2 * np.pi * abs(freq))**2
        dispffti = waffti / (2 * np.pi * abs(freq))**2
    else:
        n = len(wafftr)
        freq[0] = 1.0
        dispfftr = wafftr / (2 * np.pi * abs(freq))
        dispffti = waffti / (2 * np.pi * abs(freq))
    
    # set zero freq to 0
    dispfftr[0] = 0
    dispffti[0] = 0

    complex_array = dispfftr + 1j*dispffti
    wadisp = np.fft.ifft(complex_array,n)

    # return time history and convert from m to mm
    return wadisp.real * 1000
