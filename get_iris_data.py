from obspy.core.utcdatetime import UTCDateTime
from obspy.clients.fdsn.client import Client
from sys import argv

'''
Code to extract IRIS data, one station at a time.  Exports mseed file to 
working directory
'''

try:
    # get inputs
    datetime = argv[1] # fmt = (Y,m,d,H,M)
    sta = argv[2].upper() # station code
    
    # format UTC datetime
    dtsplit = datetime[1:-1].split(',')
    utcdt = '-'.join((dtsplit[0], dtsplit[1].zfill(2), dtsplit[2].zfill(2))) \
            + 'T' + ':'.join((dtsplit[3].zfill(2), dtsplit[4].zfill(2), '00.000'))
    
    client = Client("IRIS")
    t1 = UTCDateTime(utcdt)
    t2 = t1 + 900
    t3 = t1 + 3
    
    bulk = [("AU", sta, "*", "*", t1, t2),
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
except:
    print '\nUsage: \n    python get_iris_data.py <datetime tuple> <station code>\n'
    
    print '    e.g.: python get_iris_data.py (2010,4,20,0,16) kmbl\n'
