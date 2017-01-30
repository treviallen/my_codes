# finds files and gets geometric mean of two horizonal componets
def get_site_geomean(stn, folder):
    from fnmatch import filter
    from os import path, walk, system

    for root, dirnames, filenames in walk(folder):
        for filename in filenames:
            if filename.find(stn+'_HNE') >= 0:
                efile = path.join(root, filename)
            elif filename.find(stn+'_SLE') >= 0:
                efile = path.join(root, filename)

            if filename.find(stn+'_HNN') >= 0:
                nfile = path.join(root, filename)
            elif filename.find(stn+'_SLN') >= 0:
                nfile = path.join(root, filename)

    # read data
    T, SAe = read_psa(efile)
    T, SAn = read_psa(nfile)

    # get geometric mean and convert to g
    geomean = exp((log(SAe) + log(SAn)) / 2.) / 9.81 # convert from m/s**2 to g

    return T, geomean

# reads psa data files and returns period (T) and acceleration (SA) vectors
def read_psa(psafile):
    from numpy import array

    lines = open(psafile).readlines()[4:]  # ignore first 4 lines

    SA = []
    T = []
    for line in lines:
        dat = line.strip().split('\t')
        T.append(float(dat[0]))
        SA.append(float(dat[1]))

    return array(T), array(SA)

'''
start main
'''
# start of plotting routine
from numpy import array, arange, sqrt, exp, log
from sys import argv
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42

from calc_oq_gmpes import crustal_gsims

folder = argv[1]
folder = 'Spectra_Seiscomp3_Bener_Meriah' # example
colTrue = argv[2]

# set event details
mag  = 6.1
dep = 10.
ztor = 7. # guess
rake = 30. # USGS CMT
dip  = 90.

# set site details
vs30 = [560., 560., 560., 560.]
rrup = array([54.2, 151.3, 164.3, 170.8])
rjb = sqrt(rrup**2 - dep**2) # assume point source; i.e. repi = rjb
stn  = ['MLSI', 'LASI', 'LHMI', 'CEMA']

fig = plt.figure(1, figsize=(10, 10))
cmap = plt.cm.get_cmap('Spectral', 7)
cs = (cmap(arange(7)))

for i, s in enumerate(stn):
    # get ground motion estimates from GMPEs
    Zea06imt, CB08imt, CY08imt, BA11imt, Aea13imt, AA13imt \
        = crustal_gsims(mag, dep, ztor, dip, rake, rrup[i], rjb[i], vs30)

    ax = plt.subplot(2, 2, i+1)
    if colTrue == 'True':
        plt.loglog(Zea06imt['per'], exp(Zea06imt['sa']), lw=2., color=[cs[0][0],cs[0][1],cs[0][2]])
        plt.loglog(CB08imt['per'], exp(CB08imt['sa']),   lw=2., color=[cs[1][0],cs[1][1],cs[1][2]])
        plt.loglog(CY08imt['per'], exp(CY08imt['sa']),   lw=2., color=[cs[2][0],cs[2][1],cs[2][2]])
        plt.loglog(BA11imt['per'], exp(BA11imt['sa']),   lw=2., color=[cs[4][0],cs[4][1],cs[4][2]])
        plt.loglog(AA13imt['per'], exp(AA13imt['sa']), lw=2., color=[cs[5][0],cs[5][1],cs[5][2]])
        plt.loglog(Aea13imt['per'], exp(Aea13imt['sa']), lw=2., color=[cs[6][0],cs[6][1],cs[6][2]])
    else:
        plt.loglog(Zea06imt['per'], exp(Zea06imt['sa']), '-', lw=2.,  color='0.35')
        plt.loglog(CB08imt['per'], exp(CB08imt['sa']),   '--', lw=2., color='0.15')
        plt.loglog(CY08imt['per'], exp(CY08imt['sa']),   '-.', lw=2., color='0.15')
        plt.loglog(BA11imt['per'], exp(BA11imt['sa']),   '-', lw=2.,  color='0.65')
        plt.loglog(Aea13imt['per'], exp(Aea13imt['sa']), '--', lw=2.,  color='0.65')

    # get recorded SeisComP3 data
    T, geomean = get_site_geomean(s, folder)
    plt.loglog(T, geomean, lw=2., color='k')


    plt.xlabel('Period (s)')
    plt.ylabel('Spectral Acceleration (g)')
    plt.xlim([0.025, 10])
    plt.title(' '.join((s+'; MW =',str("%0.1f" % mag)+'; Rrup =',str(rrup[i]),'km')))
    plt.grid(which='both', color='0.5')

    if i == 0:
        plt.legend(['Zea06', 'CB08','CY08','BA11','Aea14','GeoMean'],loc=3)
        '''
        leg = plt.gca().get_legend()
        ltext  = leg.get_texts()
        plt.setp(ltext, fontsize='small')
        '''

plt.savefig(folder+'.pdf', format='pdf', dpi=150)
plt.show()

