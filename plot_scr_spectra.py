# finds files and gets geometric mean of two horizonal componets
def get_site_geomean(stn, folder):
    from fnmatch import filter
    from os import path, walk, system

    for root, dirnames, filenames in walk(folder):
        for filename in filenames:
            # check strong-motion
            if filename.find(stn) >= 0 and filename.find('NE.psa') >= 0:
                efile = path.join(root, filename)
                tfile = efile.split('NE.psa')
                nfile = ''.join((tfile[0],'NN.psa',tfile[1]))
            # check velocity sensors
            elif filename.find(stn) >= 0 and filename.find('HE.psa') >= 0:
                efile = path.join(root, filename)
                tfile = efile.split('HE.psa')
                nfile = ''.join((tfile[0],'HN.psa',tfile[1]))
                
            # get Z file
            elif filename.find(stn) >= 0 and filename.find('NZ.psa') >= 0:
                zfile = path.join(root, filename)
            elif filename.find(stn) >= 0 and filename.find('HZ.psa') >= 0:
                zfile = path.join(root, filename)
            

    try:
        try:
            # read data 
            T, SAe = read_psa(efile)
            T, SAn = read_psa(nfile)
            print efile, nfile
            
            # read psa deatails
            esta, esps, erhyp, epga, epgv, mag, dep = read_psa_details(efile)
            nsta, nsps, nrhyp, npga, npgv, mag, dep = read_psa_details(nfile)
            
            pga = max([epga, npga]) / (1000. * 9.81)
            rhyp = erhyp
            
            # get geometric mean and convert to g
            geomean = exp((log(SAe) + log(SAn)) / 2.) / 9.81 # convert from m/s**2 to g
        # just E-comp
        except:
            # read data
            print efile
            T, geomean = read_psa(efile)
            geomean = geomean / 9.81
            
            # read psa deatails
            sta, sps, rhyp, pga, pgv, mag, dep = read_psa_details(efile)
            pga = pga / (1000. * 9.81)
        
    except:
        # read data
        print zfile
        T, geomean = read_psa(zfile)
        geomean = geomean / 9.81
        
        # read psa deatails
        sta, sps, rhyp, pga, pgv, mag, dep = read_psa_details(zfile)
        pga = pga / (1000. * 9.81)
        
    return T, geomean, pga, rhyp

# reads psa data files and returns period (T) and acceleration (SA) vectors
def read_psa(psafile):
    from numpy import array

    lines = open(psafile).readlines()[24:]  # ignore first 24 lines

    SA = []
    T = []
    for line in lines:
        dat = line.strip().split('\t')
        T.append(float(dat[0]))
        SA.append(float(dat[1]))

    return array(T), array(SA)
    
# get other details for plotting
def read_psa_details(psafile):
    lines = open(psafile).readlines()[0:23]

    for line in lines:
        dat = line.strip().split('\t')
        if line.find('STA') >= 0:
            sta = dat[1]
        if line.find('SR') >= 0:
            sps = float(dat[1].split()[0])
        if line.find('RHYP') >= 0:
            rhyp = float(dat[1].split()[0])
        if line.find('PGA') >= 0:
            pga = float(dat[1].split()[0]) / 9.81
        if line.find('PGV') >= 0:
            pgv = float(dat[1].split()[0])
        if line.find('EQMAG') >= 0:
            mag = float(dat[1].split()[0])
        if line.find('EQDEP') >= 0:
            dep = float(dat[1].split()[0])
            
    return sta, sps, rhyp, pga, pgv, mag, dep

# calculate GMPEs and make subplots        
def makesubplt(i, fig, plt, sta, sps, mag, dep, ztor, dip, rake, rhyp):
    import matplotlib 
    matplotlib.rc('xtick', labelsize=8) 
    matplotlib.rc('ytick', labelsize=8) 
    from misc_tools import get_log_xy_locs
    import matplotlib.patheffects as PathEffects
    
    rrup = rhyp
    rjb = sqrt(rrup**2 - dep**2) # assume point source; i.e. repi = rjb
    
    # get station vs30
    vs30, isproxy = get_station_vs30(sta.split('.')[0])
    
    # get ground motion estimates from GMPEs
    Tea02imt, C03imt, AB06imt, Sea09imt, Sea09YCimt, Pea11imt, A12imt, Bea14imt, YA15imt, SP16imt \
        = scr_gsims(mag, dep, ztor, dip, rake, rrup, rjb, vs30)
        
    # lets adjust site to appropriate VS30
    Sea09YCimt = adjust_gmm_with_SS14(Sea09YCimt, 865., vs30)
    A12imt_hold = A12imt # for testing GM adjustment
    A12imt = adjust_gmm_with_SS14(A12imt, 820., vs30)
    
    ax = plt.subplot(3, 3, i)
    if colTrue == True:
        plt.loglog(Tea02imt['per'], exp(Tea02imt['sa']), lw=1.5, color=cs[0])
        #plt.loglog(C03imt['per'], exp(C03imt['sa']),     lw=1.5, color=cs[1])
        plt.loglog(AB06imt['per'], exp(AB06imt['sa']),   lw=1.5, color=cs[1])
        #plt.loglog(CY08imt['per'], exp(CY08imt['sa']),   lw=1.5, color=cs[3])
        plt.loglog(Sea09YCimt['per'], exp(Sea09YCimt['sa']), lw=1.5, color=cs[2])
        plt.loglog(Pea11imt['per'], exp(Pea11imt['sa']), lw=1.5, color=cs[3])
        plt.loglog(A12imt['per'], exp(A12imt['sa']),     lw=1.5, color=cs[4])
        #plt.loglog(A12imt_hold['per'], exp(A12imt_hold['sa']),  ls='--',   lw=1.5, color='k') # for testing GM adjustment
        plt.loglog(Bea14imt['per'], exp(Bea14imt['sa']), lw=1.5, color=cs[5])
        plt.loglog(YA15imt['per'], exp(YA15imt['sa']), lw=1.5, color=cs[6])
        plt.loglog(SP16imt['per'], exp(SP16imt['sa']), lw=1.5, color=cs[7])
    else:
        plt.loglog(Tea02imt['per'], exp(Tea02imt['sa']), '-', lw=1.5,  color='0.35')
        plt.loglog(C03imt['per'], exp(C03imt['sa']),     '--', lw=1.5, color='0.15')
        plt.loglog(AB06imt['per'], exp(AB06imt['sa']),   '-.', lw=1.5, color='0.15')
        plt.loglog(CY08imt['per'], exp(CY08imt['sa']),   '-', lw=1.5,  color='0.65')
        plt.loglog(Sea09Cimt['per'], exp(Sea09Cimt['sa']), '--', lw=1.5,  color='0.65')
        plt.loglog(Pea11imt['per'], exp(Pea11imt['sa']), '-.', lw=1.5,  color='0.65')
        plt.loglog(A12imt['per'], exp(A12imt['sa']),     '--', lw=1.5,  color='0.35')
        
    # get recorded process_waves.py psa data
    T, geomean, pga, rhyp = get_site_geomean(sta[0:-1], folder)
    plt.loglog(T, geomean, lw=1.5, color='k')
    
    #plt.xlabel('Period (s)', fontsize=9)
    if i >= 7:
        plt.xlabel('Period (s)', fontsize=9)
    if i == 1 or i == 4 or i == 7:
        plt.ylabel('Spectral Acceleration (g)', fontsize=9)
        
    plt.xlim([0.02, 5])
    plt.xlim([0.02, 5])
    #plt.title(' '.join((sta+'; MW =',str("%0.1f" % mag)+'; Rrup =',str(rrup),'km')), fontsize=9)
    plt.title(sta, fontsize=9)
    plt.grid(which='both', color='0.5')
    
    # annotate
    lims = ax.get_ylim()
    yloc = get_log_xy_locs(lims, 0.04)
    lims = ax.get_xlim()
    xloc = get_log_xy_locs(lims, 0.96)
    
    txt = 'Rhyp = '+str(rrup)+' km\n'+'Vs30 = '+str('%0.0f' % vs30)+' m/s'
    if isproxy == True:
        txt += '*'
    txth = plt.text(xloc, yloc, txt, ha='right', va='bottom', fontsize=8)
    txth.set_path_effects([PathEffects.withStroke(linewidth=4, foreground='w')])

    if i == 1:
        #plt.legend(['T02', 'C03','AB06','Sea09YC','Pea11','A12','Bea14', 'YA15', 'GeoMean'],loc=3, fontsize=7)
        plt.legend(['T02', 'AB06','Sea09YC','Pea11','A12','Bea14', 'YA15', 'SP16', 'GeoMean'],loc=3, fontsize=6.5)
        '''
        leg = plt.gca().get_legend()
        ltext  = leg.get_texts()
        plt.setp(ltext, fontsize='small')
        '''

'''
start main
'''
# start of plotting routine
from numpy import array, arange, sqrt, exp, log, unique, argsort
from sys import argv
from os import getcwd
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
import matplotlib as mpl
mpl.style.use('classic')

from calc_oq_gmpes import scr_gsims, adjust_gmm_with_SS14, get_station_vs30
from gmt_tools import cpt2colormap

folder = argv[1]
prefix = argv[2]
#colTrue = argv[3]
colTrue = True

# set event details
'''
mag  = 5.0
dep = 11.
'''
ztor = 0. # guess
rake = 90. # USGS CMT
dip  = 30.

# set site details
#vs30 = 760

ii = 1
fig = plt.figure(ii, figsize=(10, 10))
#cmap = plt.cm.get_cmap('Spectral', 7)
ncols = 9
if getcwd().startswith('/nas'):
    cptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/cpt/gay-flag-1978.cpt'
else:
    cptfile = '/Users/tallen/Documents/DATA/GMT/cpt/gay-flag-1978.cpt'
#cptfile = 'U:\\DATA\\GMT\\cpt\\gay-flag-1978.cpt'
cmap, zvals = cpt2colormap(cptfile, ncols)
cs = (cmap(arange(ncols)))


from fnmatch import filter
from os import path, walk, system

# build sites
sites = []
for root, dirnames, filenames in walk(folder):
    for filename in filenames:
        if filename.find('psa') >= 0:
            tmpsite = filename.split('.')[1]
            tmpcomp = filename.split('.')[2]
            sites.append('.'.join((tmpsite,tmpcomp)))

# if two H components, rename channel
#usites = unique(sites)
usites = []
for site1 in sites:
    cnt = 0
    for site2 in sites:
        if site1[0:-1] == site2[0:-1]:
            cnt += 1
            print site1
            
    if cnt == 2 or cnt == 3:
        usites.append(site1[0:-1]+'H')
    elif cnt == 1:
        usites.append(site1)
        
usites = unique(usites)

udists = []
for stn in usites:
    for root, dirnames, filenames in walk(folder):
        for filename in filenames:
            if filename.find(stn[0:-1]) >= 0:
                if filename.find('NE') >= 0:
                    psafile = path.join(root, filename)
                elif filename.find('HE') >= 0:
                    psafile = path.join(root, filename)
                elif filename.find('NZ') >= 0:
                    psafile = path.join(root, filename)
                elif filename.find('HZ') >= 0:
                    psafile = path.join(root, filename)
                elif filename.find('NH') >= 0:
                    psafile = path.join(root, filename)
                elif filename.find('HHH') >= 0:
                    psafile = path.join(root, filename)
                
    # get record details
    sta, sps, rhyp, pga, pgv, mag, dep = read_psa_details(psafile)
    udists.append(rhyp)
    
    
# now sort by distance
udists = array(udists)
idx=argsort(array(udists))
udists = udists[idx]
usites = usites[idx]

# loop thru sites ordered by distance
i = 0
for stn in usites:
    for root, dirnames, filenames in walk(folder):
        for filename in filenames:
            if filename.find(stn[0:-1]) >= 0:
                if filename.find('NE') >= 0:
                    psafile = path.join(root, filename)
                elif filename.find('HE') >= 0:
                    psafile = path.join(root, filename)
                elif filename.find('NZ') >= 0:
                    psafile = path.join(root, filename)
                elif filename.find('HZ') >= 0:
                    psafile = path.join(root, filename)
                elif filename.find('NH') >= 0:
                    psafile = path.join(root, filename)
                elif filename.find('HHH') >= 0:
                    psafile = path.join(root, filename)
                
    # get record details
    sta, sps, rhyp, pga, pgv, mag, dep = read_psa_details(psafile)
    
    #mag = 4.
    
    # now plot
    # make sub plot
    if stn != 'CDNM.HNH':
        i += 1
        makesubplt(i, fig, plt, stn, sps, mag, dep, ztor, dip, rake, rhyp)
        
        if i == 9:
            i = 0
            plt.savefig(prefix+'.Lake_Muir_'+str(ii)+'.png', format='png', dpi=300, bbox_inches='tight')
            ii += 1
            fig = plt.figure(ii, figsize=(10, 10))

plt.savefig(prefix+'.Lake_Muir_'+str(ii)+'.png', format='png', dpi=300, bbox_inches='tight')
plt.show()

