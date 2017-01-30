# finds files and gets geometric mean of two horizonal componets
def get_site_geomean(stn, folder, t):
    from fnmatch import filter
    from os import path, walk, system
    from numpy import amax

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

    '''
    # get max of H component
    try:
        # read data
        T, SAe = read_psa(efile)
        T, SAn = read_psa(nfile)

        maxsa = []
        for i in range(0,len(SAe)):
            maxsa.append(max([SAe[i], SAn[i]]) / 100.)
    '''
    # get geomean
    try:
        T, SAe = read_psa(efile)
        T, SAn = read_psa(nfile)

        # get geometric mean and convert to g
        if efile.find('.psa5') >= 0:
            print efile
            geomean = exp((log(SAe) + log(SAn)) / 2.) / 9.81 # convert from m/s**2 to g
        else:
            geomean = exp((log(SAe) + log(SAn)) / 2.)

        # now get max val
        maxsa = []
        for i in range(0,len(SAe)):
            maxsa.append(max([SAe[i], SAn[i]]) / 9.81)

    except:
        T, geomean = read_psa(efile)

    if t == 0.0:
        geomean = maxsa

    return T, geomean

# reads psa data files and returns period (T) and acceleration (SA) vectors
def read_psa(psafile):
    from numpy import array

    lines = open(psafile).readlines()[3:]  # only use first 4 lines

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
from numpy import array, arange, sqrt, exp, log, log10, logspace, argwhere, interp
from sys import argv
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42

from calc_oq_gmpes import crustal_gsims #, get_T_index

folder = argv[1]
folder = 'Spectra_Seiscomp3_Bener_Meriah' # example
colTrue = argv[2]

# read stn & distance
sitefile = 'bener_meriah_lookup.txt'
lines = open(sitefile).readlines()

stn = []
dist = []
for line in lines:
    dat = line.strip().split('\t')
    stn.append(dat[0])
    dist.append(float(dat[1]))

# set event details
mag  = 6.1
dep = 10.
ztor = 7. # guess
rake = 30. # USGS CMT
dip  = 90.

# set site details
vs30 = [560.]
rjb = logspace(0,2.5,50)
rrup = sqrt(rjb**2 + dep**2) # assume point source; i.e. repi = rjb
#stn  = ['MLSI', 'LASI', 'LHMI', 'CEMA']

fig = plt.figure(1, figsize=(10, 10))
cmap = plt.cm.get_cmap('Spectral', 6)
cs = (cmap(arange(6)))

titles = ['PGA','SA(0.3)','SA(1.0)','SA(3.0)']

Tplot = [0.0, 0.3, 1.0, 3.0]
#Tplot = [0.0]
# loop thru periods
for j, t in enumerate(Tplot):
    ax = plt.subplot(2, 2, j+1)
    Zea06r = []
    CB08r = []
    CY08r = []
    BA11r = []
    Aea13r = []
    for i,r in enumerate(rrup):

        # get ground motion estimates from GMPEs
        Zea06imt, CB08imt, CY08imt, BA11imt, Aea13imt \
            = crustal_gsims(mag, dep, ztor, dip, rake, rrup[i], rjb[i], vs30)

        if t == 0.0:
            Zea06r.append(Zea06imt['pga'])
            CB08r.append(CB08imt['pga'])
            CY08r.append(CY08imt['pga'])
            BA11r.append(BA11imt['pga'])
            Aea13r.append(Aea13imt['pga'])
        else:
            #ti = get_T_index(Zea06imt, t)
            # interpolate log values to correct period
            Zea06r.append(interp(t, Zea06imt['per'], Zea06imt['sa']))
            CB08r.append(interp(t, CB08imt['per'], CB08imt['sa']))
            CY08r.append(interp(t, CY08imt['per'], CY08imt['sa']))
            BA11r.append(interp(t, BA11imt['per'], BA11imt['sa']))
            Aea13r.append(interp(t, Aea13imt['per'], Aea13imt['sa']))

    if colTrue == 'True':
        plt.loglog(rjb, exp(Zea06r), lw=2., color=[cs[0][0],cs[0][1],cs[0][2]])
        plt.loglog(rjb, exp(CB08r),  lw=2., color=[cs[1][0],cs[1][1],cs[1][2]])
        plt.loglog(rjb, exp(CY08r),  lw=2., color=[cs[2][0],cs[2][1],cs[2][2]])
        plt.loglog(rjb, exp(BA11r),  lw=2., color=[cs[4][0],cs[4][1],cs[4][2]])
        plt.loglog(rjb, exp(Aea13r), lw=2., color=[cs[5][0],cs[5][1],cs[5][2]])
    else:
        plt.loglog(rjb, exp(Zea06r), '-', lw=2.,  color='0.15')
        plt.loglog(rjb, exp(CB08r),  '--', lw=2., color='0.15')
        plt.loglog(rjb, exp(CY08r),  '-.', lw=2., color='0.15')
        plt.loglog(rjb, exp(BA11r),  '-', lw=2.,  color='0.55')
        plt.loglog(rjb, exp(Aea13r), '.', lw=2.,  color='0.55')

    # get recorded SeisComP3 data and plot
    for k, s in enumerate(stn):
        try:
            T, geomean = get_site_geomean(s, folder, t)
            # get interpolated value
            specval = interp(t, T, geomean)
            plt.loglog(dist[k], specval, 'k+', markersize=10, markeredgewidth=1.75)

        except:
            print 'Cannot find data for site:', s


    plt.xlabel('Rjb (km)')
    plt.ylabel('Spectral Acceleration (g)')
    plt.xlim([3, 300])
    plt.ylim([1E-4, 1])
    plt.title(titles[j])
    plt.grid(which='both', color='0.5')

    if j == 0:
        plt.legend(['Zea06', 'CB08','CY08','BA11','Aea14','Data'],loc=3)
        '''
        leg = plt.gca().get_legend()
        ltext  = leg.get_texts()
        plt.setp(ltext, fontsize='small')
        '''

plt.savefig(folder+'_rjb.pdf', format='pdf', dpi=150)
plt.show()

