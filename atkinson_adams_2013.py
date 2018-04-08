# converts cm/s^2 to g
def cmpss2g(cmpss):
    return (cmpss / 100) / 9.81

# Interpolates across AA13 look-up tables
def atkinson_adams_2013(mw, repi, **kwargs):

    from numpy import array, where, log10, interp, unique, round, shape, log, ones, exp

    '''
    requires mw & repi

    valid "crust_ty" values are:
        offshore
        ena
        inslab
        interface
        wcrust
    '''

    '''
    # for testing
    mw = 6.5
    rjb = 125.89
    repi = 125.89
    '''

    # set default crustal type
    crust_ty = 'wcrust'

    for key in (['crust_ty']):
       if key in kwargs:
           if key == 'crust_ty':
               crust_ty = kwargs[key]

    if crust_ty == 'offshore':
        gmpetab = '/Users/tallen/Documents/Code/pycode/gmpe/GMPEt_offshore_med.dat'
    elif crust_ty == 'wcrust':
        gmpetab = '/Users/tallen/Documents/Code/pycode/gmpe/GMPEt_Wcrust_med.dat'
    elif crust_ty == 'ena':
        gmpetab = '/Users/tallen/Documents/Code/pycode/gmpe/GMPEt_ENA_med.dat'
    elif crust_ty == 'inslab':
        gmpetab = '/Users/tallen/Documents/Code/pycode/gmpe/GMPEt_Inslab_med.dat'
    elif crust_ty == 'interface':
        gmpetab = '/Users/tallen/Documents/Code/pycode/gmpe/GMPEtInterface_med_combo.dat'
        
    # set arrays
    mwlist = []
    rjblist = []
    repilist = []
    savals = []

    # read appropriate data
    lines = open(gmpetab, 'rb').readlines()

    # get period array
    pertxt  = lines[7].strip().split()[3:]
    freqs = [float(x) for x in pertxt]

    # set mw, rjb, repi & psa values
    for line in lines[8:]:
        dat = line.strip().split()
        mwlist.append(float(dat[0]))
        rjblist.append(float(dat[1]))
        repilist.append(float(dat[2]))
        savals.append([float(x) for x in dat[3:]])
    
    if crust_ty == 'interface':
        dr = 0.15
        savals = 10**array(savals)
    else:
        dr = 0.1 
        

    # get mag indices
    mwlist = array(mwlist)
    rjblist = array(rjblist)
    repilist = array(repilist)
    savals = array(savals)
    index = where((mwlist >= mw-0.25) & (mwlist <= mw+0.25) \
                  & (log10(repilist) >= round(log10(repi)-dr,decimals=2)) \
                  & (log10(repilist) <= round(log10(repi)+dr,decimals=2)))[0]
    #print repilist[index], mwlist[index]
    #print log10(repilist)

    # make temp arrays
    mwstrip = mwlist[index]
    rjbstrip = rjblist[index]
    repistrip = repilist[index]
    sastrip = log10(savals[index])

    # for each Repi, interpolate spectral values across magnitude
    urepi = unique(repistrip)
    interpsa = []
    interpJB = []
    
    for ue in urepi:
        index2 = where((repistrip == ue))[0]
        for i in range(0, len(index2)):
            # make sa vect
            tmpinterpsa = []
            for j in range(len(sastrip[0])):
                tmpsa = []
                for k in range(0, len(index2)):
                    tmpsa.append(sastrip[index2[k]][j])

                # now interpolate across Mw
                intval = interp(mw, mwstrip[index2], tmpsa)
                tmpinterpsa.append(intval)

        interpsa.append(array(tmpinterpsa))
        interpJB.append(interp(mw, mwstrip[index2], rjbstrip[index2]))

    # now we have interpolated across magnitude, interpolate across Repi for each T
    finalsa = []
    for t in range(0,len(interpsa[0])):
        tmpsa = []
        for k in range(0, shape(interpsa)[0]):
            tmpsa.append(interpsa[k][t])

        intval = 10**(interp(log(repi), log(urepi), tmpsa))
        finalsa.append(intval)

    # get equivalent JB distance - this is needed for wcrust model
    equivJBdist = interp(repi, urepi, interpJB)
    
    pga = log(finalsa[-2] / 981.) # in g
    pgv = finalsa[-1] # in cm/s?
    sa  = log(array(finalsa[:-2]) / (100 * 9.81))
    T   = 1. / array(freqs[:-2])

    # set constant sigma
    sig = ones(len(sa)) * log(10**0.27)
    
    #AA13imt = {'pga':pga, 'pgv':pgv, 'sa':sa, 'per':T}
    AA13imt = {'pga':pga, 'pgv':pgv, 'sa':sa[::-1], 'per':T[::-1], 'sig':sig, 'eqivJB':equivJBdist}
    	
    #print exp(AA13imt['sa'])*981., equivJBdist

    return AA13imt













