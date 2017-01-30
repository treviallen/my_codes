from optparse import OptionParser
import os
import fnmatch
from datetime import datetime
import struct
from obspy.core import Trace, Stream, UTCDateTime
import numpy as np
from progressbar import Bar, ETA, FileTransferSpeed,  Percentage, ProgressBar, RotatingMarker


parser = OptionParser(usage="%prog -f ")
parser.add_option("-f","--rootfolder", type="string", dest="folder", help="folder")

(options, args) = parser.parse_args()

if options.folder == None:
        parser.error("i need folder exception")


def readFolderStructur(rootfolder):

    FileList = []

    for root,dirs,files in os.walk(rootfolder):

            for i in files:
                if fnmatch.fnmatch(i,'*.w*'):
                    fname = os.path.join(root,i)
                    FileList.append(fname)

    return FileList

def unpackFile(fobjpath):

    try:
        os.system("gunzip " + fobjpath)
    except:
        'uncompressing exception'

def convert(FileList):

    ca = {"BZ" : "BHZ", "BN":"BHN", "BE":"BHE","HZ" : "HHZ", "HN":"HHN", "HE":"HHE","LZ" : "LHZ", "LN":"LHN", "LE":"LHE"}

    #widgets = ['Something: ', Percentage(), ' ', Bar(marker=RotatingMarker()),' ', ETA(), ' ', FileTransferSpeed()]
    #pbar = ProgressBar(widgets=widgets,maxval=len(FileList))

    print 'Uncompressing Data\n'
    c=1;
    for f in FileList:

        if fnmatch.fnmatch(f,'*.gz'):
            unpackFile(f)

        #pbar.update(c)
        c+=1

    #pbar.finish()

    counter = 1;


    print 'Convert Data\n'
    for i in FileList:


        if fnmatch.fnmatch(i, '*wfdisc'):
            print "\n----------------------------------------------------------"
            print "\033[31m   WFDISCFILE: "+str(counter)+' of '+str(len(FileList)/2)+' : '+i+'\033[0m'
            path = os.path.dirname(i)

            fobj = open(i,'r')
            for line in fobj:
                print 'LAENGE: ',len(line)

                if len(line) == 284:

                    print line

                    head = map(str.strip, line.split())
                    inputf = os.path.join(path,head[16]) # input file name

                    if fnmatch.fnmatch(inputf,'*.w'):

                            print "INPUTFILE: "+inputf
                            timedata = map(int, head[2].split(".")) # 0 -> timestamp, 1-> millisecs

                                #headers for SH ascii file
                            SH = {
                                   "STATION": head[0],
                                   "COMP": head[1][-1].upper(), # guess from channel naming
                                  "START": ".".join((datetime.fromtimestamp(timedata[0]).strftime("%d-%b-%Y_%H:%M:%S"), "%03d" % timedata[1])),
                                  "DELTA": 1/float(head[8]),
                                  "LENGTH": int(head[7]),
                                  "CALIB": float(head[9])
                                  }

                            # binary format (big endian integers)
                            fmt = ">"+"i"*SH["LENGTH"]

                            # convert binary data
                            try:
                                bf=open(inputf, "rb")
                                shift = int(head[17])
                                bf.seek(shift)
                                data = struct.unpack(fmt, bf.read(struct.calcsize(fmt)))

                                # write miniseed
                                data=np.array(data,dtype=np.int32)
                                network=" "
                                station = head[0]
                                location=" "

                                channel=" "
                                if len(head[1]) == 2:
                                    for key in ca:
                                        if key == head[1]:
                                            channel =ca[key]
                                else:
                                    channel = head[1]

                                #channel=head[1]
                                npts=len(data)
                                sampling_rate=head[8]

                                outputname = inputf[:-2]+'.'+channel
                                stats = {'network': network, 'station': station, 'location': location,
                                         'channel': channel, 'npts': npts, 'sampling_rate': sampling_rate,
                                         'mseed' : {'dataquality' : 'D'}}
                                stats['starttime'] = UTCDateTime(timedata[0])
                                st = Stream([Trace(data=data, header=stats)])
                                outfile = os.path.join(path,outputname+".mseed")
                                print "OUTPUTFILE |-"+outfile
                                st.write(outfile, format='MSEED')
                            except:
                                print 'File not readable ----> ', inputf

                    else:
                        print '\033[31m  File %s not found \033[0m' %(inputf)


            counter+=int(1)



if __name__ == '__main__':

    L = readFolderStructur(options.folder)

    convert (L)