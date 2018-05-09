# converts g to cm/s**2
def g2cgs(gin):
    from scipy.constants import g
    from numpy import array
    return array(gin) * 100. * g


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

    try:
        imt['pga'] = gmpe.get_mean_and_stddevs(sites, rup, dists, PGA(), [StdDev.TOTAL])
        try:
            imt['pgv'] = gmpe.get_mean_and_stddevs(sites, rup, dists, PGV(), [StdDev.TOTAL])
        except:
            imt['pgv'] = nan
    except:
        imt['pga'] = nan
        imt['pgv'] = nan

    imt['sa'] = []
    imt['sig'] = []
    imt['inter'] = []
    imt['intra'] = []
    
    for p in imt['per']:
        
        satmp = array([gmpe.get_mean_and_stddevs(sites, rup, dists, SA(p), [StdDev.TOTAL])[0][0]])
        '''
        try:
            satmp = array([gmpe.get_mean_and_stddevs(sites, rup, dists, SA(p), [StdDev.TOTAL])[0][0]])
        except:
            print gmpe.get_mean_and_stddevs(sites, rup, dists, SA(p), [StdDev.TOTAL])
        '''

        try:
            if isnan(satmp[0]):
                imt['sa'].append(satmp[1])
            else:
                imt['sa'].append(satmp[0])
        except:
            if isnan(satmp[0][0]):
                imt['sa'].append(satmp[0][1])
            else:
                imt['sa'].append(satmp[0][0])
                
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
    from openquake.hazardlib.gsim.zhao_2006 import ZhaoEtAl2006SSlab, ZhaoEtAl2006SSlabCascadia 
    from openquake.hazardlib.gsim.atkinson_boore_2003 import AtkinsonBoore2003SSlab, AtkinsonBoore2003SSlabCascadiaNSHMP2008
    from openquake.hazardlib.gsim.youngs_1997 import YoungsEtAl1997SSlab
    from openquake.hazardlib.gsim.garcia_2005 import GarciaEtAl2005SSlab
    #from openquake.hazardlib.gsim.abrahamson_2015 import AbrahamsonEtAl2015SSlab
    from openquake.hazardlib.gsim.megawati_pan_2010 import MegawatiPan2010
    #from openquake.hazardlib.gsim.zhao_2016 import ZhaoEtAl2016SSlab #, ZhaoEtAl2006SSlabCascadia 
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
    sites.backarc = [0.]
    
    rup = RuptureContext()
    rup.mag = mag
    rup.hypo_depth = dep
    rup.dip = dip
    rup.rake = rake
    rup.ztor = rup.hypo_depth

    dists = DistancesContext()
    dists.rrup = array([rrup])
    dists.rhypo = array([rrup]) # assume only interested in closest piont to the fault
    dists.rjb = array([rjb])
    dists.rx = sqrt(dists.rrup**2 - rup.hypo_depth**2) # this is not correct, but good enough for now
    dists.rvolc = array([100.])
    #print dists.rrup

    gmpe = ZhaoEtAl2006SSlab()
    Zea06imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)
    
    gmpe = ZhaoEtAl2006SSlabCascadia()
    Zea06CISimt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)

    gmpe = AtkinsonBoore2003SSlab()
    AB03imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)
    
    gmpe = AtkinsonBoore2003SSlabCascadiaNSHMP2008()
    AB03CISimt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)

    gmpe = YoungsEtAl1997SSlab()
    Yea97imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)
    
    gmpe = GarciaEtAl2005SSlab()
    Gea05imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)
    
    gmpe = MegawatiPan2010()
    MP10imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)
    
    '''
    gmpe = AbrahamsonEtAl2015SSlab()
    Aea15imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)
    
    gmpe = ZhaoEtAl2016SSlab()
    Zea16imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)
    '''

    # prepare Atkinson & Adams 2013
    #repi = sqrt(rrup**2 - dep**2)
    #AA13imt = atkinson_adams_2013(mag, dists.rjb[0], crust_ty = crust_ty)

    return Yea97imt, AB03imt, AB03CISimt, Gea05imt, Zea06imt, Zea06CISimt, MP10imt #, Aea15imt #, Zea16imt #, AA13imt #, Aea15imt, Zea06CISimt, 

# calls and calculates candidate interface GMPEs - values returned in ln(g)
def interface_gsims(mag, dep, ztor, dip, rake, rrup, rjb, vs30):
    from openquake.hazardlib.gsim.zhao_2006 import ZhaoEtAl2006SInter, ZhaoEtAl2006SInterCascadia
    from openquake.hazardlib.gsim.ghofrani_atkinson_2014 import GhofraniAtkinson2014, GhofraniAtkinson2014Cascadia
    from openquake.hazardlib.gsim.atkinson_boore_2003 import AtkinsonBoore2003SInter
    from openquake.hazardlib.gsim.youngs_1997 import YoungsEtAl1997SInter
    from openquake.hazardlib.gsim.abrahamson_2015 import AbrahamsonEtAl2015SInter
    from openquake.hazardlib.gsim.atkinson_macias_2009 import AtkinsonMacias2009
    from openquake.hazardlib.gsim.megawati_pan_2010 import MegawatiPan2010
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
    dists.rhypo = array([rrup])
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
    
    gmpe = MegawatiPan2010()
    MP10imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)
    
    gmpe = GhofraniAtkinson2014()
    GA14imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)
    
    gmpe = GhofraniAtkinson2014Cascadia()
    GA14CISimt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)
    
    # prepare Atkinson & Adams 2013
    #print 'interface',  mag, dists.rjb[0]
    #AA13imt = atkinson_adams_2013(mag, dists.rrup[0], crust_ty = crust_ty) # note AA13 model uses Rrup
    
    return Yea97imt, AB03imt, Zea06imt, Zea06CISimt, AM09imt, MP10imt, GA14imt, GA14CISimt, Aea15imt

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
    from openquake.hazardlib.gsim.shahjouei_pezeshk_2016 import ShahjoueiPezeshk2016
    from atkinson_adams_2013 import atkinson_adams_2013
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
    
    gmpe = ShahjoueiPezeshk2016()
    SP16imt = get_pga_sa(gmpe, sites, rup, dists, crust_ty)
    
    crust_ty = 'ena'
    
    #AA13imt = atkinson_adams_2013(mag, dists.rjb[0], crust_ty = crust_ty)
    #AA13imt = []

    return Tea02imt, C03imt, AB06imt, CY08imt, Sea09imt, Sea09YCimt, Pea11imt, A12imt, Bea14imt, YA15imt, SP16imt # AA13imt, 

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
    
def gaull1990_gsim(mag, dep, rhypo):
    
    from openquake_local.hazardlib.gsim.gaull_1990 import GaullEtAL1990WesternAustralia, \
                                                    GaullEtAL1990SoutheasternAustralia, \
                                                    GaullEtAL1990PGAfromPGVWesternAustralia, \
                                                    GaullEtAL1990PGAfromPGVSoutheasternAustralia
    '''
    from gaull_1990 import GaullEtAL1990WesternAustralia, \
                                                    GaullEtAL1990SoutheasternAustralia, \
                                                    GaullEtAL1990PGAfromPGVWesternAustralia, \
                                                    GaullEtAL1990PGAfromPGVSoutheasternAustralia
    '''                                               
    from openquake.hazardlib.gsim.base import RuptureContext, SitesContext, DistancesContext
    from openquake.hazardlib.imt import PGA
    from openquake.hazardlib.const import StdDev
    from numpy import array
    
    sites = SitesContext()

    rup = RuptureContext()
    rup.mag = mag # input MW - GaullEtAL1990 corrects to ML
    rup.hypo_depth = dep

    dists = DistancesContext()
    dists.rhypo = array([rhypo])
    
    gmpe = GaullEtAL1990WesternAustralia()
    G90WAimt = {'pga': gmpe.get_mean_and_stddevs(sites, rup, dists, PGA(), [StdDev.TOTAL])}

    gmpe = GaullEtAL1990SoutheasternAustralia()
    G90SEAimt = {'pga': gmpe.get_mean_and_stddevs(sites, rup, dists, PGA(), [StdDev.TOTAL])}
    
    gmpe = GaullEtAL1990PGAfromPGVWesternAustralia()
    G90WA_PGVimt = {'pga': gmpe.get_mean_and_stddevs(sites, rup, dists, PGA(), [StdDev.TOTAL])} 
    
    gmpe = GaullEtAL1990PGAfromPGVSoutheasternAustralia()
    G90SEA_PGVimt = {'pga': gmpe.get_mean_and_stddevs(sites, rup, dists, PGA(), [StdDev.TOTAL])}
    
    return G90WAimt, G90SEAimt, G90WA_PGVimt, G90SEA_PGVimt
    
def hdf5_gsim(mag, dep, ztor, dip, rake, rrup, rjb, rhypo, vs30, hdf5file):
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
        
    gmpe = GMPETable(gmpe_table = hdf5file) # use full path
    hdf5imt = get_pga_sa_gmpe(gmpe, sites, rup, dists, crust_ty)
    
    return hdf5imt

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

# reads sa data files and returns period (T) and acceleration (SA) vectors
def read_sa(safile):
    from numpy import array
    from scipy.constants import g

    lines = open(safile).readlines()
    chan = lines[0].split('.')[-2]
    datestr = lines[1].split('\t')[-1]    
    sta = lines[2].split('\t')[-1]
    stla = float(lines[3].split('\t')[-1])
    stlo = float(lines[4].split('\t')[-1])
    sr = float(lines[5].split('\t')[-1].split()[0])
    eqla = float(lines[6].split('\t')[-1])
    eqlo = float(lines[7].split('\t')[-1])
    dep = float(lines[8].split('\t')[-1].split()[0])
    mag = float(lines[9].split('\t')[-1])
    rhyp = float(lines[10].split('\t')[-1].split()[0])
    azim = float(lines[11].split('\t')[-1].split()[0])
    pga = float(lines[12].split('\t')[-1].split()[0]) / (1000. * g) # convert mm/s**2 to g
    pgv = float(lines[13].split('\t')[-1].split()[0]) / 10. # convert mm/s to cm/s

    SA = []
    T = []
    
    for line in lines[24:]:
        dat = line.strip().split('\t')
        T.append(float(dat[0]))
        SA.append(float(dat[1]))
    
    T = array(T)
    SA = array(SA) / g

    rec = {'chan':chan, 'datestr':datestr, 'sta':sta, 'stla':stla, 'stlo':stlo, \
           'sr':sr, 'eqla':eqla, 'eqlo':eqlo, 'dep':dep, 'mag':mag, 'rhyp':rhyp, \
           'azim':azim, 'pga':pga, 'pgv':pgv, 'per':T, 'sa':SA}
    
    return rec

# script to find extrapolation ratio based on input gmm
def get_extrap_ratio(extrapDat, targetDat, maxPer, extrapPer):
    '''
    assume extrapDat in log space
    extrapDat = GMM from which to determine ratio at target extrpolation period
    targetDat = GMM for which to extrapolate GMM for
    maxPer =  maximum period
    extrapPer = period to extrapolate to
    '''
    
    from numpy import interp, log
    
    # get interpolation ratio from host GMM
    maxGM  = interp(log(maxPer), log(extrapDat['per']), extrapDat['sa'])
    extrapGM = interp(log(extrapPer), log(extrapDat['per']), extrapDat['sa'])
    logRat = maxGM - extrapGM
    
    # now return extrapolated GM
    return targetDat['sa'][-1] - logRat

# scrift to make gmm tables for updating GMMs
def gsim2table(gmmClass, gmmName, mags, distances, depth, vs30, vs30ref, extrapPeriods, interpPeriods, rtype, folder):
    '''
    gmmName = OQ GMM string (e.g. 'GhofraniAtkinson2014Cascadia')
    gmmClass = model class - parsed from calling gsim2table
    mags = numpy array of magnitudes
    dists = numpy array of distances
    vs30 = target vs30 (if required)
    vs30ref = reference vs30 for GMM - if using GMM with vs30 parameterisation, set vs30ref = vs30
    extrapPeriod = list of periods to extrapolate to (in log-log space)
    interpPeriods = True/False - log-log interpolation to NBCC Periods
    rtype = distance metric string (e.g. rrup, ryhpo, rjb)
    folder = folder in which tables are written
    '''
    
    #eval('from openquake.hazardlib.gsim.'+gmmPy+' import '+gmmClass)
    from openquake.hazardlib.gsim.base import RuptureContext, SitesContext, DistancesContext
    from numpy import array, sqrt, log, log10, exp, interp
    from seyhan_stewart_2014 import seyhan_stewart_siteamp
    from atkinson_boore_site_2006 import atkinson_boore_siteamp
    from os import path
    
    ss14vsref = 760. # m/s
    ab06vsref = 760. # m/s
    
    if gmmName == 'AtkinsonMacias2009' or gmmName == 'GhofraniAtkinson2014Cascadia' \
       or gmmName == 'ZhaoEtAl2006SInterCascadia' or gmmName == 'AbrahamsonEtAl2015SInter':
        crust_ty = 'interface'
        # import extrap GMM
        from openquake.hazardlib.gsim.abrahamson_2015 import AbrahamsonEtAl2015SInter
        gmpeExtrap = AbrahamsonEtAl2015SInter()
        
    elif gmmName == 'GarciaEtAl2005SSlab' or gmmName == 'ZhaoEtAl2006SSlabCascadia' \
       or gmmName == 'AtkinsonBoore2003SSlabCascadia' or gmmName == 'AbrahamsonEtAl2015SSlab':
        crust_ty = 'inslab'
        # import extrap GMM
        from openquake.hazardlib.gsim.abrahamson_2015 import AbrahamsonEtAl2015SSlab
        gmpeExtrap = AbrahamsonEtAl2015SSlab()
        
    
    tabtxt = ''
    for i, m in enumerate(mags):
        for d in distances:
            sites = SitesContext()
            sites.vs30 = array([float(vs30)])
            sites.vs30measured = 0
            #sites.z1pt0 = exp(28.5 - (3.82/8.)*log(sites.vs30**8 + 378.7**8)) # in m; from ChiouYoungs2008
            sites.z1pt0 = exp((-7.15 / 4.)*log((sites.vs30**4 + 571.**4) / (1360.**4 + 571.**4))) # in m; from ChiouYoungs2014
            sites.z2pt5 = (519 + 3.595 * sites.z1pt0) / 1000. #in km; from Kaklamanos etal 2011
            sites.backarc = [False] # assume not backarc
            
            rup = RuptureContext()
            rup.mag = m
            rup.hypo_depth = depth
            #rup.dip = dip
            #rup.rake = rake
            #rup.ztor = rup.hypo_depth
            #rup.width = mag2rupwid_WC94(mag, 'all') # should be checked
            
            # Assume all distances are equivalent - distance metric will be stored in tbale header for further interp
            dists = DistancesContext()
            dists.rrup = array([d])
            dists.rjb = array([d])
            dists.rhypo = array([d]) 
            dists.rvolc = array([100.]) # assume backarc distance of 100 km
            
            # get SA for given mag & distance
            gmmDat = get_pga_sa(gmmClass, sites, rup, dists, crust_ty)
            
            #######################################################################################
            
            # do site class correction if needed
            '''
            note: this is a simplification - to do properly where vs30ref > 760 and vs30 < vsref, 
                  should first correct to 760 using AB06, then amplify by SS14
                  
                  I think this is now taken care of, but will need to check!
            '''
            # if vs30 == vs30ref, then using GMM-native amp factors
            if vs30 != vs30ref:
                
                if vs30 < vs30ref:
                    # first SS14 amplification factors to correct to GMM vs30ref - assume no more than 1100 m/s
                    tmpamps = []
                    for t in gmmDat['per']:
                        tmpamps.append(seyhan_stewart_siteamp(vs30ref, t, exp(gmmDat['pga'][0]))[0])
                    
                    vsrefSAcorr = array(tmpamps)
                    vsrefPGAcorr = seyhan_stewart_siteamp(vs30ref, 0.0, exp(gmmDat['pga'][0]))
                    vsrefPGVcorr = seyhan_stewart_siteamp(vs30ref, -1.0, exp(gmmDat['pga'][0]))
                    
                    # correct PGA at GMM refernce vs30 to SS14 760 m/s
                    refPGA_SS14 = log(exp(gmmDat['pga'][0]) / vsrefPGAcorr)
                    
                    # now get SS14 amplification factors from 760 to target vs30 m/s
                    tmpamps = []
                    for t in gmmDat['per']:
                        tmpamps.append(seyhan_stewart_siteamp(vs30, t, exp(refPGA_SS14))[0]) 
                    
                    vstargSAcorr = array(tmpamps)
                    vstargPGAcorr = seyhan_stewart_siteamp(vs30, 0.0, exp(refPGA_SS14))
                    vstargPGVcorr = seyhan_stewart_siteamp(vs30, -1.0, exp(refPGA_SS14))
                    
                # else, use AB06 for harder sites
                else:
                    # first get AB06 amplification factors to correct to GMM vs30ref - assume no more than 1100 m/s
                    tmpamps = []
                    for t in gmmDat['per']:
                        #print 'PGA1', t, m, d, vs30ref, vs30, gmmDat['pga'][0], atkinson_boore_siteamp(vs30ref, t, exp(gmmDat['pga'][0]))
                        try:
                            tmpamps.append(atkinson_boore_siteamp(vs30ref, t, exp(gmmDat['pga'][0]))[0])
                        except:
                            tmpamps.append(atkinson_boore_siteamp(vs30ref, t, exp(gmmDat['pga'][0]))) # have no idea why i have to do this!
                    
                    vsrefSAcorr = array(tmpamps)
                    vsrefPGAcorr = atkinson_boore_siteamp(vs30ref, 0.0, exp(gmmDat['pga'][0]))
                    vsrefPGVcorr = atkinson_boore_siteamp(vs30ref, -1.0, exp(gmmDat['pga'][0]))
                    
                    # correct PGA at given vs30 to GMM refernce vs30
                    refPGA_AB06 = log(exp(gmmDat['pga'][0]) / vsrefPGAcorr)
                    
                    # now get AB06 amplification factors from 760 to target vs30 m/s
                    tmpamps = []
                    for t in gmmDat['per']:
                        #print 'PGA2', t, m, d, vs30ref, vs30, refPGA_AB06, atkinson_boore_siteamp(vs30, t, exp(refPGA_AB06))
                        tmpamps.append(atkinson_boore_siteamp(vs30, t, exp(refPGA_AB06)[0]))
                    
                    vstargSAcorr = array(tmpamps)
                    vstargPGAcorr = atkinson_boore_siteamp(vs30, 0.0, exp(refPGA_AB06))
                    vstargPGVcorr = atkinson_boore_siteamp(vs30, -1.0, exp(refPGA_AB06))
                        
                # now correct gmmDat
                gmmDat['sa'] = list(log(exp(gmmDat['sa']) * vstargSAcorr / vsrefSAcorr))
                gmmDat['pga'][0][0] = log(exp(gmmDat['pga'][0][0]) * vstargPGAcorr / vsrefPGAcorr)
                
                try:
                    doPGV = True
                    gmmDat['pgv'][0][0] = log(exp(gmmDat['pgv'][0][0]) * vstargPGVcorr / vsrefPGVcorr)
                except:
                    doPGV = False
            
            #######################################################################################
            
            # determine if extrapolation should be performed - if so extrapolate based on well-known GMMs (empeExtrap)
            for ep in extrapPeriods:
                if ep > max(gmmDat['per']):
                    #print maxPer, ep
                    maxPer = max(gmmDat['per'])
                    # get extrapolated sa
                    extrapDat = get_pga_sa(gmpeExtrap, sites, rup, dists, crust_ty)
                    extrapSA = get_extrap_ratio(extrapDat, gmmDat, maxPer, ep)
                    
                    gmmDat['per'].append(ep)
                    gmmDat['sa'].append(extrapSA)
                    gmmDat['sig'].append(gmmDat['sig'][-1]) # extrapolate sigma
                    
            
            #######################################################################################
            # interpolate to standard NBCC periods
            if interpPeriods == True:
                targetPeriods = array([0.1, 0.2, 0.3, 0.5, 1., 2., 5., 10.])
                
                # interpolate SA - already in natural log
                gmmDat['sa']  = interp(log(targetPeriods), log(gmmDat['per']), gmmDat['sa'])
                
                # interpolate sigma - already in natural log
                gmmDat['sig'] = interp(log(targetPeriods), log(gmmDat['per']), gmmDat['sig'])
                
                # now reset periods
                gmmDat['per'] = targetPeriods
                
            #######################################################################################
            
            # set text for mag/dist
            # convert ln g to cm/s**2 as required for table builder
            sa = log10(g2cgs(exp(gmmDat['sa'])))
            
            #print gmmDat['sa']
            sastr = ' '.join([str('%0.3f' % x) for x in sa])
            
            #######################################################################################
            
            # get log10 Sa in g
            sag = log10(exp(gmmDat['sa']))
            
            # estimate PGV in case needed
            # get coefs for calculating PGV - PGV is log cm/s, Sa is log g
            if crust_ty == 'interface':
                # use ~2 sec
                c0 = 0.897
                c1 = 2.10
                
                sa20 = interp(log10(2.0), log10(gmmDat['per']), sag)  # check logs are ok here
                proxyPGV = c0 * sa20 + c1
                
                # use Sa2.0 sigma as a proxy for PGV sigma (keep in ln)
                proxyPGVsigma = interp(log10(2.0), log10(gmmDat['per']), gmmDat['sig'])
                
            else:
                # for other crust types based on BSSA at SA 0.5 sec
                c0 = 1.07
                c1 = 1.83
                    
                sa05 = interp(log10(0.5), log10(gmmDat['per']), sag)
                proxyPGV = c0 * sa05 + c1
                
                # use Sa0.5 sigma as a proxy for PGV sigma (keep in ln)
                proxyPGVsigma = interp(log10(0.5), log10(gmmDat['per']), gmmDat['sig'])
            
            try:
                tabtxt += ' '.join((str('%0.2f' % m), str('%0.2f' % d))) + ' ' + sastr + ' ' \
                          + str('%0.3f' % log10(g2cgs(exp(gmmDat['pga'][0])))) + ' ' \
                          + str('%0.3f' % log10(exp(gmmDat['pgv'][0][0]))) + '\n'
                doPGV = True
                
            except:
                tabtxt += ' '.join((str('%0.2f' % m), str('%0.2f' % d))) + ' ' + sastr + ' ' \
                          + str('%0.3f' % log10(g2cgs(exp(gmmDat['pga'][0])))) + ' ' \
                          + str('%0.3f' % proxyPGV) + '\n'
                doPGV = False
                      
    # get header info
    header  = ' '.join((gmmName, 'as implemented in OpenQuake, distance is', rtype+'.','Log10 hazard values in cgs units\n'))
    
    if doPGV == True:
        header += ' '.join(('         ', str(len(mags)), str(len(distances)), str(len(sa)+2), ': nmag, ndist, nperiod')) + '\n'
        header += ' '.join(('         ', ' '.join([str('%0.3f' % x) for x in gmmDat['per']]), 'PGA PGV')) + '\n'
        header += ' '.join(('         ', ' '.join([str('%0.3f' % x) for x in gmmDat['sig']]), \
                            str('%0.3f' % gmmDat['pga'][1][0]), str('%0.3f' % gmmDat['pgv'][1][0]))) + '\n' # natural log
                            
    else:
        header += ' '.join(('         ', str(len(mags)), str(len(distances)), str(len(sa)+1), ': nmag, ndist, nperiod')) + '\n'
        header += ' '.join(('         ', ' '.join([str('%0.3f' % x) for x in gmmDat['per']]), 'PGA PGV')) + '\n'
        header += ' '.join(('         ', ' '.join([str('%0.3f' % x) for x in gmmDat['sig']]), \
                            str('%0.3f' % gmmDat['pga'][1][0]), str('%0.3f' % proxyPGVsigma))) + '\n' # natural log
    
    # write to file
    print '\nWriting table:', '.'.join((gmmName,'vs'+str(int(vs30)),'h'+str(int(depth)),'txt'))
    filename = '.'.join((gmmName,'vs'+str(int(vs30)),'h'+str(int(depth)),'txt'))
    f = open(path.join(folder, filename), 'wb')
    f.write(header+tabtxt)
    f.close()
    return tabtxt, gmmDat