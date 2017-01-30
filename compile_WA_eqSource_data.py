# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 15:33:56 2015

@author: tallen
"""
from os import path, listdir

evdatafile = 'WA.event_pars.txt'

lines = open(evdatafile).readlines()

events = []
for line in lines:
    dat = line.strip().split()
    evdict = {'date': dat[0], 'hhmm': dat[1], 'lon': float(dat[2]), \
              'lat': float(dat[3]), 'dep': float(dat[4]), \
              'mw': float(dat[6])}
    evdict['evid'] = evdict['date'].replace('-','') + evdict['hhmm']
              
    events.append(evdict)

''' now look for stns used '''
# parse eqSource event dat files
recs = []
for event in events:
    evdatfile = path.join('..','..','eqSource','Events',event['evid'][0:4],event['evid'][4:6], \
                           event['evid'][6:8] + '_' + event['evid'][8:] + '.dat')
    # now read event data file
    evlines = open(evdatfile).readlines()
    
    # now get sites used
    readsites = False
    for el in evlines:
        if el.startswith('\n'):
            readsites = False
        
        # get sites here
        if readsites == True:
            stndat = el.strip().split()
            if stndat[0].endswith('A') or stndat[0].endswith('V'):
                rec = {'evid': event['evid'], 'stn': stndat[0][0:-1], \
                        'rawstn': stndat[0], 'rhyp': float(stndat[1])}
            else:
                rec = {'evid': event['evid'], 'stn': stndat[0], \
                       'rawstn': stndat[0], 'rhyp': float(stndat[1])}
            recs.append(rec)
            
        if el.startswith('Site'):
            readsites = True
    
''' now look for files in waves folder '''
for rec in recs:
    # look in path for file
    evpath = path.join('..','..','eqSource','Waves',rec['evid'][0:4],rec['evid'][4:6], \
                        rec['evid'][6:8] + '_' + rec['evid'][8:])
    
    files = listdir(evpath)
    for f in files:
        if f.find(rec['stn']) >= 0:
            if f.endswith('.txt') or f.endswith('.G') or f.endswith('.S'):
                rec['wavpath'] = path.join('..','..','eqSource','Waves',rec['evid'][0:4], \
                                      rec['evid'][4:6],rec['evid'][6:8] + '_' + rec['evid'][8:],f)
                
''' now parse data files '''
from readwaves import return_data 

for i, rec in enumerate(recs):
    # check file format
    allsta, recs[i]['comps'], allrecdate, allsec, recs[i]['sps'], recs[i]['data'], \
    recs[i]['nsamp'], fmt = return_data(rec['wavpath'])
    print fmt
