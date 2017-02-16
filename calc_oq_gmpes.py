# does the IM calculation given input GMPE, source and site info
def get_pga_sa(gmpe, sites, rup, dists, crust_ty):
    from openquake.hazardlib.imt import PGA, PGV, SA
    from openquake.hazardlib.const import StdDev
    from numpy import array, isnan, nan
    '''
    ###########################################################################
    Note - can get periods from imtls!!!
    
    imtl.keys() =  ['pga', 'per', 'sig', 'sa', 'inter', 'intra']
    imtl['per'] = periods!
    ###########################################################################
    '''
        
    
    # get gmpe periods
    if crust_ty == 'wcrust':
        try:
            pkeys = gmpe.COEFFS.sa_coeffs.keys()
        except:
            pkeys = gmpe.COEFFS_ASC.sa_coeffs.keys()
            
    elif crust_ty == 'inslab':
        try:
            pkeys = gmpe.COEFFS.sa_coeffs.keys()
        except:
            try:
                pkeys = gmpe.COEFFS_ASC.sa_coeffs.keys()
            except:
                try:
                    pkeys = gmpe.COEFFS_SSLAB.sa_coeffs.keys()
                except:
                    try:
                        pkeys = gmpe.COEFFS_ROCK.sa_coeffs.keys()
                    except:
                        pkeys = gmpe.COEFFS_SOIL.sa_coeffs.keys()
                
    elif crust_ty == 'interface':
        try:
            pkeys = gmpe.COEFFS.sa_coeffs.keys()
        except:
            try:
                pkeys = gmpe.COEFFS_ASC.sa_coeffs.keys()
            except:
                try:
                    pkeys = gmpe.COEFFS_SINTER.sa_coeffs.keys()
                except:
                    try:
                        pkeys = gmpe.COEFFS_ROCK.sa_coeffs.keys()
                    except:
                        pkeys = gmpe.COEFFS_SOIL.sa_coeffs.keys()
                    
    elif crust_ty == 'intraplate':
        try:          
            pkeys = gmpe.COEFFS_HARD_ROCK.sa_coeffs.keys()
        except:
            try:
                pkeys = gmpe.COEFFS.sa_coeffs.keys()
            except:
                if rup.hypo_depth > 10.:
                    pkeys = gmpe.COEFFS_DEEP.sa_coeffs.keys()
                else:
                     pkeys = gmpe.COEFFS_SHALLOW.sa_coeffs.keys()

    imt = {}
    periods = []
    for k in pkeys:
       periods.append(k[1])
    periods.sort()
    imt['per'] = periods

    imt['pga'] = gmpe.get_mean_and_stddevs(sites, rup, dists, PGA(), [StdDev.TOTAL])
    try:
        imt['pgv'] = gmpe.get_mean_and_stddevs(sites, rup, dists, PGV(), [StdDev.TOTAL])
    except:
        imt['pgv'] = nan
    imt['sa'] = []
    imt['sig'] = []
    imt['inter'] = []
    imt['intra'] = []
    for p in imt['per']:
        satmp = array([gmpe.get_mean_and_stddevs(sites, rup, dists, SA(p), [StdDev.TOTAL])[0][0]])
        if isnan(satmp[0]):
            imt['sa'].append(satmp[1])
        else:
            imt['sa'].append(satmp[0])
        imt['sig'].append(gmpe.get_mean_and_stddevs(sites, rup, dists, SA(p), [StdDev.TOTAL])[1][0][0])
        try:
            imt['inter'].append(gmpe.get_mean_and_stddevs(sites, rup, dists, SA(p), [StdDev.INTER_EVENT])[1][0][0])
            imt['intra'].append(gmpe.get_mean_and_stddevs(sites, rup, dists, SA(p), [StdDev.INTRA_EVENT])[1][0][0])
        except:
            imt['inter'].append(nan)
            imt['intra'].append(nan) 

    return imt
    
def get_pga_sa_gmpe(gmpe, sites, rup, dists, crust_ty):
    from openquake.hazardlib.imt import PGA, SA
    from openquake.hazardlib.const import StdDev
    from numpy import array, isnan, nan, rollaxis, swapaxes
    
    imt = {}
    periods = []
    periods = gmpe.imls[u'T']
    imt['per'] = periods
    
    imt['pga'] =  gmpe.get_mean_and_stddevs(sites, rup, dists, PGA(), [StdDev.TOTAL])
    imt['sa'] = []
    for p in imt['per']:
        satmp = gmpe.get_mean_and_stddevs(sites, rup, dists, SA(p), [StdDev.TOTAL])
        if isnan(satmp[0]):
            imt['sa'].append(satmp[1])
        else:
            imt['sa'].append(satmp[0])

    return imt

# calls and calculates candidate active crustal GMPEs - values returned in ln(g)
def crustal_gsims(mag, dep, ztor, dip, rake, rrup, rjb, vs30):
    from openquake.hazardlib.gsim.boore_1997 import BooreEtAl1997GeometricMean   
    from openquake.hazardlib.gsim.zhao_2006 import ZhaoEtAl2006Asc
    from openquake.hazardlib.gsim.campbell_bozorgnia_2008 import CampbellBozorgnia2008
    from openquake.hazardlib.gsim.chiou_youngs_2008 import ChiouYoungs2008
    from openquake.hazardlib.gsim.bindi_2011 import BindiEtAl2011
    from openquake.hazardlib.gsim.boore_atkinson_2011 import BooreAtkinson2011
    from openquake.hazardlib.gsim.akkar_2014 import AkkarEtAlRjb2014
    from openquake.hazardlib.gsim.boore_2014 import BooreEtAl2014
    from openquake.hazardlib.gsim.campbell_bozorgnia_2014 import CampbellBozorgnia2014
    from openquake.hazardlib.gsim.chiou_youngs_2014 import ChiouYoungs2014
    from openquake.hazardlib.gsim.base import RuptureContext, SitesContext, DistancesContext
    from atkinson_adams_2013 import atkinson_adams_2013
    from fault_tools import mag2rupwid_WC94
    from numpy import array, sqrt, log, exp
    
    crust_ty = 'wcrust'

    sites = SitesContext()
    sites.vs30 = array([float(vs30)])
    sites.vs30measured = 0
    #sites.z1pt0 = exp(28.5 - (3.82/8.)*log(sites.vs30**8 + 378.7**8)) # in m; from ChiouYoungs2008
    sites.z1pt0 = exp((-7.15 / 4.)*log((sites.vs30**4 + 571.**4) / (1360.**4 + 571.**4))) # in m; from ChiouYoungs2014
    sites.z2pt5 = (519 + 3.595 * sites.z1pt0) / 1000. #in km; from Kaklamanos etal 2011
    

    rup = RuptureContext()
    rup.mag = mag
    rup.hypo_depth = dep
    rup.dip = dip
    rup.rake = rake
    rup.ztor = rup.hypo_depth
    rup.width = mag2rupwid_WC94(mag, 'all') # should be checked

    dists = DistancesContext()
    dists.rrup = array([rrup])
    dists.rjb = array([rjb])    
    dists.rx = sqrt(dists.rrup**2 - rup.hypo_depth**2) # this is not correct, but good enough for now

    gmpe = BooreEtAl1997GeometricMean()
    Bea97imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)    
    
    gmpe = ZhaoEtAl2006Asc()
    Zea06imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)

    gmpe = CampbellBozorgnia2008()
    CB08imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)

    gmpe = ChiouYoungs2008()
    CY08imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)

    gmpe = BindiEtAl2011()
    Bea11imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)  
    
    gmpe = BooreAtkinson2011()
    BA11imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)

    gmpe = AkkarEtAlRjb2014()
    Aea14imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)
    
    gmpe = BooreEtAl2014()
    Bea14imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)
    
    gmpe = CampbellBozorgnia2014()
    CB14imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)
    
    gmpe = ChiouYoungs2014()
    CY14imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)
    
    # prepare Atkinson & Adams 2013
    repi = sqrt(rrup**2 - dep**2) # not correct if rrup != rhypo
    #print 'crust', mag, dists.rjb[0]
    AA13imt = atkinson_adams_2013(mag, dists.rjb[0], crust_ty = crust_ty) # assume Rjb = Repi
    '''
    return Bea97imt, Zea06imt, CB08imt, CY08imt, Bea11imt, BA11imt, AA13imt, Aea14imt, Bea14imt, CY14imt, \
           sites, rup, dists
    '''       
    return Bea97imt, Zea06imt, CB08imt, CY08imt, Bea11imt, BA11imt, AA13imt, Aea14imt, Bea14imt, CB14imt, CY14imt

# calls and calculates candidate intraslab GMPEs - values returned in ln(g)
def inslab_gsims(mag, dep, ztor, dip, rake, rrup, rjb, vs30):
    from openquake.hazardlib.gsim.zhao_2006 import ZhaoEtAl2006SSlab
    from openquake.hazardlib.gsim.atkinson_boore_2003 import AtkinsonBoore2003SSlab
    from openquake.hazardlib.gsim.youngs_1997 import YoungsEtAl1997SSlab
    from atkinson_adams_2013 import atkinson_adams_2013
    from openquake.hazardlib.gsim.base import RuptureContext, SitesContext, DistancesContext
    from numpy import array, sqrt, log, exp
    
    crust_ty = 'inslab'
    
    sites = SitesContext()
    sites.vs30 = array([float(vs30)])
    sites.vs30measured = False
    #sites.z1pt0 = exp(28.5 - (3.82/8.)*log(sites.vs30**8 + 378.7**8)) # in m; from ChiouYoungs2008
    sites.z1pt0 = exp((-7.15 / 4.)*log((sites.vs30**4 + 571.**4) / (1360.**4 + 571.**4))) # in m; from ChiouYoungs2014
    sites.z2pt5 = (519 + 3.595 * sites.z1pt0) / 1000. #in km; from Kaklamanos etal 2011

    rup = RuptureContext()
    rup.mag = mag
    rup.hypo_depth = dep
    rup.dip = dip
    rup.rake = rake
    rup.ztor = rup.hypo_depth

    dists = DistancesContext()
    dists.rrup = array([rrup])
    dists.rjb = array([rjb])
    dists.rx = sqrt(dists.rrup**2 - rup.hypo_depth**2) # this is not correct, but good enough for now

    gmpe = ZhaoEtAl2006SSlab()
    Zea06imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)

    gmpe = AtkinsonBoore2003SSlab()
    AB03imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)

    gmpe = YoungsEtAl1997SSlab()
    Yea97imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)

    # prepare Atkinson & Adams 2013
    #repi = sqrt(rrup**2 - dep**2)
    AA13imt = atkinson_adams_2013(mag, dists.rjb[0], crust_ty = crust_ty)

    return Yea97imt, AB03imt, Zea06imt, AA13imt

# calls and calculates candidate interface GMPEs - values returned in ln(g)
def interface_gsims(mag, dep, ztor, dip, rake, rrup, rjb, vs30):
    from openquake.hazardlib.gsim.zhao_2006 import ZhaoEtAl2006SInter, ZhaoEtAl2006SInterCascadia
    from openquake.hazardlib.gsim.ghofrani_atkinson_2014 import GhofraniAtkinson2014, GhofraniAtkinson2014Cascadia
    from openquake.hazardlib.gsim.atkinson_boore_2003 import AtkinsonBoore2003SInter
    from openquake.hazardlib.gsim.youngs_1997 import YoungsEtAl1997SInter
    from openquake.hazardlib.gsim.abrahamson_2015 import AbrahamsonEtAl2015SInter
    from openquake.hazardlib.gsim.atkinson_macias_2009 import AtkinsonMacias2009
    from atkinson_adams_2013 import atkinson_adams_2013
    from openquake.hazardlib.gsim.base import RuptureContext, SitesContext, DistancesContext
    from numpy import array, sqrt, log, exp
    
    crust_ty = 'interface'
    
    '''
    # test data
    vs30 = [760.]
    mag = 6.5
    dep = 10.
    dip = 45.
    rake = 90.
    rrup = 100.
    '''
    
    sites = SitesContext()
    sites.vs30 = array([float(vs30)])
    sites.vs30measured = False
    sites.z1pt0 = exp((-7.15 / 4.)*log((sites.vs30**4 + 571.**4) / (1360.**4 + 571.**4))) # in m; from ChiouYoungs2014
    sites.z2pt5 = (519 + 3.595 * sites.z1pt0) / 1000. #in km; from Kaklamanos etal 2011
    sites.backarc = [0.]
    
    rup = RuptureContext()
    rup.mag = mag
    rup.hypo_depth = dep
    rup.dip = dip
    rup.rake = rake
    rup.ztor = rup.hypo_depth
    
    dists = DistancesContext()
    dists.rrup = array([rrup])
    dists.rjb = array([rjb])
    dists.rx = sqrt(dists.rrup**2 - rup.hypo_depth**2) # this is not correct, but good enough for now
    
    gmpe = ZhaoEtAl2006SInter()
    # gmpe = ZhaoEtAl2006SSlab() # for testing only
    Zea06imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)
    
    gmpe = ZhaoEtAl2006SInterCascadia()
    # gmpe = ZhaoEtAl2006SSlab() # for testing only
    Zea06CISimt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)
    
    gmpe = AtkinsonBoore2003SInter()
    AB03imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)
    
    gmpe = YoungsEtAl1997SInter()
    Yea97imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)
    
    gmpe = AbrahamsonEtAl2015SInter()
    Aea15imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)
    
    gmpe = AtkinsonMacias2009()
    AM09imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)
    
    gmpe = GhofraniAtkinson2014()
    GA14imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)
    
    gmpe = GhofraniAtkinson2014Cascadia()
    GA14CISimt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)
    
    # prepare Atkinson & Adams 2013
    #print 'interface',  mag, dists.rjb[0]
    AA13imt = atkinson_adams_2013(mag, dists.rrup[0], crust_ty = crust_ty) # note AA13 model uses Rrup
    
    return Yea97imt, AB03imt, Zea06imt, Zea06CISimt, AM09imt, AA13imt, GA14imt, GA14CISimt, Aea15imt

# calls and calculates candidate SCR GMPEs
def scr_gsims(mag, dep, ztor, dip, rake, rrup, rjb, vs30):
    from openquake.hazardlib.gsim.toro_2002 import ToroEtAl2002
    from openquake.hazardlib.gsim.campbell_2003 import Campbell2003
    from openquake.hazardlib.gsim.atkinson_boore_2006 import AtkinsonBoore2006
    from openquake.hazardlib.gsim.chiou_youngs_2008 import ChiouYoungs2008
    from openquake.hazardlib.gsim.somerville_2009 import SomervilleEtAl2009NonCratonic
    from openquake.hazardlib.gsim.somerville_2009 import SomervilleEtAl2009YilgarnCraton
    from openquake.hazardlib.gsim.pezeshk_2011 import PezeshkEtAl2011
    from openquake.hazardlib.gsim.allen_2012 import Allen2012
    from openquake.hazardlib.gsim.boore_2014 import BooreEtAl2014
    from openquake.hazardlib.gsim.yenier_atkinson_2015 import YenierAtkinson2015CEUS
    from openquake.hazardlib.gsim.base import RuptureContext, SitesContext, DistancesContext
    from numpy import array, sqrt, log, exp
    
    crust_ty = 'intraplate'
    
    sites = SitesContext()
    sites.vs30 = array([float(vs30)])
    sites.vs30measured = False
    sites.z1pt0 = exp((-7.15 / 4.)*log((sites.vs30**4 + 571.**4) / (1360.**4 + 571.**4))) # in m; from ChiouYoungs2014
    sites.z2pt5 = (519 + 3.595 * sites.z1pt0) / 1000. #in km; from Kaklamanos etal 2011
    

    rup = RuptureContext()
    rup.mag = mag
    rup.hypo_depth = dep
    rup.dip = dip
    rup.rake = rake
    rup.ztor = rup.hypo_depth
    rup.stress_drop = 300. # bar for YA15 CEUS

    dists = DistancesContext()
    dists.rrup = array([rrup])
    dists.rjb = array([rjb])
    dists.rx = sqrt(dists.rrup**2 - rup.hypo_depth**2) # this is not correct, but good enough for now

    gmpe = ToroEtAl2002()
    Tea02imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)    
    Tea02imt['sa'] = Tea02imt['sa'][:-2] # remove 3-4 s as not originally defined by author
    Tea02imt['per'] = Tea02imt['per'][:-2]

    gmpe = Campbell2003()
    C03imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)

    gmpe = AtkinsonBoore2006()
    AB06imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)
    
    gmpe = ChiouYoungs2008()
    CY08imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)
    
    gmpe = SomervilleEtAl2009NonCratonic()
    Sea09imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)

    gmpe = SomervilleEtAl2009YilgarnCraton()
    Sea09YCimt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)

    gmpe = PezeshkEtAl2011()
    Pea11imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)

    gmpe = Allen2012()
    A12imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)
    
    gmpe = BooreEtAl2014()
    Bea14imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)
    
    gmpe = YenierAtkinson2015CEUS()
    YA15imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)

    return Tea02imt, C03imt, AB06imt, CY08imt, Sea09imt, Sea09YCimt, Pea11imt, A12imt, Bea14imt, YA15imt

def allen2012_gsim(mag, dep, rrup):
    from openquake.hazardlib.gsim.allen_2012 import Allen2012
    from openquake.hazardlib.gsim.base import RuptureContext, SitesContext, DistancesContext
    from numpy import array
    
    crust_ty = 'intraplate'
    
    sites = SitesContext()

    rup = RuptureContext()
    rup.mag = mag
    rup.hypo_depth = dep

    dists = DistancesContext()
    dists.rrup = array([rrup])
    
    gmpe = Allen2012()
    A12imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)
    
    return A12imt
    
def jdfn_gsim(mag, dep, ztor, dip, rake, rrup, rjb, rhypo, vs30):
    from openquake.hazardlib.gsim.gsim_table import GMPETable
    from openquake.hazardlib.gsim.base import RuptureContext, SitesContext, DistancesContext
    from numpy import array, sqrt, log, exp

    crust_ty = 'gmpetable'
    
    sites = SitesContext()
    sites.vs30 = array([float(vs30)])
    sites.vs30measured = False
    sites.z1pt0 = exp((-7.15 / 4.)*log((sites.vs30**4 + 571.**4) / (1360.**4 + 571.**4))) # in m; from ChiouYoungs2014
    sites.z2pt5 = (519 + 3.595 * sites.z1pt0) / 1000. #in km; from Kaklamanos etal 2011

    rup = RuptureContext()
    rup.mag = mag
    rup.hypo_depth = dep
    rup.dip = dip
    rup.rake = rake
    rup.ztor = rup.hypo_depth

    dists = DistancesContext()
    dists.rrup = rrup
    dists.rjb = rjb
    dists.rhypo = rhypo
    dists.rx = sqrt(dists.rrup**2 - rup.hypo_depth**2) # this is not correct, but good enough for now
        
    gmpe = GMPETable(gmpe_table = "/home/trandolph/GSCFRISK_GMPE_tables/WinslabD30_med_clC.hdf5")
    jdfnimt = get_pga_sa_gmpe(gmpe, sites, rup, dists, crust_ty)
    
    return jdfnimt

def aa13_gsims(mag, dep, rrup, rjb, rhypo, vs30):
    from openquake.hazardlib.gsim.gsim_table import GMPETable
    from atkinson_adams_2013 import atkinson_adams_2013
    from openquake.hazardlib.gsim.base import RuptureContext, SitesContext, DistancesContext
    from numpy import array, sqrt, log, exp

    sites = SitesContext()
    sites.vs30 = array([float(vs30)])
    sites.vs30measured = False
    sites.z1pt0 = exp((-7.15 / 4.)*log((sites.vs30**4 + 571.**4) / (1360.**4 + 571.**4))) # in m; from ChiouYoungs2014
    sites.z2pt5 = (519 + 3.595 * sites.z1pt0) / 1000. #in km; from Kaklamanos etal 2011

    rup = RuptureContext()
    rup.mag = mag
    rup.hypo_depth = dep

    dists = DistancesContext()
    dists.rrup = rrup
    dists.rjb = rjb
    dists.rhypo = rhypo
    dists.rx = sqrt(dists.rrup**2 - rup.hypo_depth**2) # this is not correct, but good enough for now
    
    """ 
    gmpe = GMPETable(gmpe_table = "/Users/tallen/Documents/GEM/openquake/oq-engine/2015_gsc_nshm/gm_tables/Wcrust_med_clC.hdf5")
    AA13_crustal = get_pga_sa_gmpe(gmpe, sites, rup, dists, crust_ty)
    
    gmpe = GMPETable(gmpe_table = "/Users/tallen/Documents/GEM/openquake/oq-engine/2015_gsc_nshm/gm_tables/Woffshore_med_clC.hdf5")
    AA13_offshore = get_pga_sa_gmpe(gmpe, sites, rup, dists, crust_ty)
    """
    
    AA13_crustal = atkinson_adams_2013(mag, rjb, crust_ty = 'wcrust') # assume Rjb = Repi
    AA13_offshore = atkinson_adams_2013(mag-0.5, rjb, crust_ty = 'wcrust') # same as offshore
    
    return AA13_crustal, AA13_offshore

# calls and calculates candidate active crustal GMPEs - values returned in ln(g)
def bssa_gsim(mag, dep, ztor, dip, rake, rrup, rjb, vs30):
    from openquake.hazardlib.gsim.boore_2014 import BooreEtAl2014
    from openquake.hazardlib.gsim.base import RuptureContext, SitesContext, DistancesContext
    from atkinson_adams_2013 import atkinson_adams_2013
    from fault_tools import mag2rupwid_WC94
    from numpy import array, sqrt, log, exp
    
    crust_ty = 'wcrust'

    sites = SitesContext()
    sites.vs30 = array([float(vs30)])
    sites.vs30measured = 0
    #sites.z1pt0 = exp(28.5 - (3.82/8.)*log(sites.vs30**8 + 378.7**8)) # in m; from ChiouYoungs2008
    sites.z1pt0 = exp((-7.15 / 4.)*log((sites.vs30**4 + 571.**4) / (1360.**4 + 571.**4))) # in m; from ChiouYoungs2014
    sites.z2pt5 = (519 + 3.595 * sites.z1pt0) / 1000. #in km; from Kaklamanos etal 2011
    

    rup = RuptureContext()
    rup.mag = mag
    rup.hypo_depth = dep
    rup.dip = dip
    rup.rake = rake
    rup.ztor = rup.hypo_depth
    rup.width = mag2rupwid_WC94(mag, 'all') # should be checked

    dists = DistancesContext()
    dists.rrup = array([rrup])
    dists.rjb = array([rjb])    
    dists.rx = sqrt(dists.rrup**2 - rup.hypo_depth**2) # this is not correct, but good enough for now

    gmpe = BooreEtAl2014()
    Bea14imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)
        
    return Bea14imt
