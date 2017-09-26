# these functions calculate local magnitudes
import numpy as np


def get_component_orient(filename):

    orient = filename.split('.')

    if orient[2][2] == 'Z':
        comp = 0
    elif orient[2][2] == 'E' or orient[2][2] == 'N':
        comp = 1

    return comp

# get Richer (1935; 1958), extended by Eiby & Muir (1968)
def calc_R35(comp, logA, repi):
    if repi <= 7.5:
        logA0 = 1.4
    elif repi > 7.5 and repi <= 12.5:
        logA0 = 1.5
    elif repi > 12.5 and repi <= 17.5:
        logA0 = 1.6
    elif repi > 17.5 and repi <= 22.5:
        logA0 = 1.7
    elif repi > 22.5 and repi <= 27.5:
        logA0 = 1.9
    elif repi > 27.5 and repi <= 32.5:
        logA0 = 2.1
    elif repi > 32.5 and repi <= 37.5:
        logA0 = 2.3
    elif repi > 37.5 and repi <= 42.5:
        logA0 = 2.4
    elif repi > 42.5 and repi <= 47.5:
        logA0 = 2.5
    elif repi > 47.5 and repi <= 52.5:
        logA0 = 2.6
    elif repi > 52.5 and repi <= 57.5:
        logA0 = 2.7
    elif repi > 57.5 and repi <= 72.5:
        logA0 = 2.8
    elif repi > 72.5 and repi <= 87.5:
        logA0 = 2.9
    elif repi > 87.5 and repi <= 105:
        logA0 = 3.0
    elif repi > 105 and repi <= 125:
        logA0 = 3.1
    elif repi > 125 and repi <= 145:
        logA0 = 3.2
    elif repi > 145 and repi <= 165:
        logA0 = 3.3
    elif repi > 165 and repi <= 185:
        logA0 = 3.4
    elif repi > 185 and repi <= 205:
        logA0 = 3.5
    elif repi > 205 and repi <= 215:
        logA0 = 3.6
    elif repi > 215 and repi <= 225:
        logA0 = 3.65
    elif repi > 225 and repi <= 245:
        logA0 = 3.7
    elif repi > 245 and repi <= 265:
        logA0 = 3.8
    elif repi > 265 and repi <= 285:
        logA0 = 3.9
    elif repi > 285 and repi <= 305:
        logA0 = 4.0
    elif repi > 305 and repi <= 325:
        logA0 = 4.1
    elif repi > 325 and repi <= 345:
        logA0 = 4.2
    elif repi > 345 and repi <= 375:
        logA0 = 4.3
    elif repi > 375 and repi <= 395:
        logA0 = 4.4
    elif repi > 395 and repi <= 425:
        logA0 = 4.5
    elif repi > 425 and repi <= 465:
        logA0 = 4.6
    elif repi > 465 and repi <= 505:
        logA0 = 4.7
    elif repi > 505 and repi <= 555:
        logA0 = 4.8
    elif repi > 555 and repi <= 605:
        logA0 = 4.9
    # Extend beyond intended distance (Eiby & Muir, 1968)
    elif repi > 605 and repi <= 625:
        logA0 = 5.0
    elif repi > 625 and repi <= 675:
        logA0 = 5.1
    elif repi > 675 and repi <= 725:
        logA0 = 5.2
    elif repi > 725 and repi <= 775:
        logA0 = 5.3
    elif repi > 775 and repi <= 825:
        logA0 = 5.4
    elif repi > 825 and repi <= 875:
        logA0 = 5.5
    elif repi > 875 and repi <= 925:
        logA0 = 5.55
    elif repi > 925 and repi <= 975:
        logA0 = 5.6
    elif repi > 975 and repi <= 1050:
        logA0 = 5.7
    elif repi > 1050 and repi <= 1150:
        logA0 = 5.8
    elif repi > 1150 and repi <= 1250:
        logA0 = 5.9
    elif repi > 1250 and repi <= 1350:
        logA0 = 6.0
    elif repi > 1350 and repi <= 1450:
        logA0 = 6.05
    elif repi > 1450 and repi <= 1550:
        logA0 = 6.15
    elif repi > 1550 and repi <= 1650:
        logA0 = 6.2

    R35 = logA + logA0

    magstr = 'R35:\t' + str("%0.1f" % R35)
    if comp == 0: # if vertical
        magstr = magstr + '*'
    #print magstr

    return R35

# Greenhalgh and Singh (1986) - Sth Australia
def calc_GS86(comp, logA, repi):
    GS86 = 0.7 + logA + 1.10 * np.log10(repi) + 0.0013 * repi

    magstr = 'GS86:\t' + str("%0.1f" % GS86)
    if comp == 1: # if horizontal
        magstr = magstr + '*'
    print magstr

    return GS86

# Hutton & Boore (1987) - SoCal
def calc_HB87(comp, logA, rhyp):
    logA0 = 1.11 * np.log10(rhyp / 100.) + 0.00189 * (rhyp - 100) + 3.0
    HB87 = logA + logA0

    magstr = 'HB87:\t' + str("%0.1f" % HB87)
    if comp == 0: # if vertical
        magstr = magstr + '*'
    print magstr

    return HB87

# Bakun & Joyner (1984) - SoCal
def calc_BJ84(comp, logA, rhyp):
    logA0 = np.log10(rhyp) + 0.00301 * rhyp + 0.70
    BJ84 = logA + logA0

    magstr = 'HB87:\t' + str("%0.1f" % BJ84)
    if comp == 0: # if vertical
        magstr = magstr + '*'
    print magstr

    return BJ84

# Gaull & Gregson (1991) - WA
def calc_GG91(comp, logA, rhyp):
    GG91 = logA + 1.137 * np.log10(rhyp) + 0.000657 * rhyp + 0.66
    
    magstr = 'GG91:\t' + str("%0.1f" % GG91)
    if comp == 1: # if horizontal
        magstr = magstr + '*'
    print magstr

    return GG91


# Michael-Lieba & Malafant (1992) - SEA Aust
def calc_MLM92(comp, logA, rhyp):
    if comp == 0: # if vertical
        MLM92 = logA + 1.34 * np.log10(rhyp / 100.) + 0.00055 * (rhyp - 100) + 3.13
    elif comp == 1: # if horizontal
        MLM92 = logA + 1.34 * np.log10(rhyp / 100.) + 0.00055 * (rhyp - 100) + 3.0

    magstr = 'MLM92:\t' + str("%0.1f" % MLM92)
    print magstr

    return MLM92

# Wilkie et al. (1996) - SE Aust
def calc_WGW96(comp, logA, rhyp):
    WGW96 = logA + 0.75 + np.log10(rhyp) + 0.0056*rhyp*np.exp(-0.0013*rhyp)
    
    magstr = 'WGW96:\t' + str("%0.1f" % WGW96)
    if comp == 1: # if horizontal
        magstr = magstr + '*'
    print magstr

    return WGW96
    
# Allen 2010 (unpub) - SE Aust
def calc_A10(comp, logA, rhyp):
    c1 = -1.5602
    c2 = -0.2792
    c3 = -2.0468
    r1 = 90.
    r2 = 150.
    if rhyp <= r1:
        logASEA = c1*np.log10(rhyp)
    elif rhyp > r1 and rhyp <= r2:
        logASEA = c1*np.log10(r1) + c2*np.log10(rhyp/r1)
    elif rhyp > r2:
        logASEA = c1*np.log10(r1) + c2*np.log10(r2/r1) + c3*np.log10(rhyp/r2)

    A10 = logA - (logASEA + 0.0619)
    
    magstr = 'A10:\t' + str("%0.1f" % A10)
    if comp == 1: # if horizontal
        magstr = magstr + '*'
    print magstr

    return A10
    
def calc_A16(logA, rhyp, comp):
    from numpy import loadtxt, log10
    coeffile = '//Users//tallen//Dropbox//Magnitudes//2016_working//A16_coeffs.txt'
    dat = loadtxt(coeffile)
    b1 = dat[0]
    b2 = dat[1]
    b3 = dat[2]
    r1 = dat[3]
    r2 = dat[4]
    c0 = dat[5] - 3.0 # norm 
    c1 = dat[6] # h/v correction - constant
    c1l = dat[7] # h/v correction - distance dependent
    c2l = dat[8] # h/v correction - distance dependent
    
    if rhyp <= r1:
        logASEA = b1*log10(rhyp) + c0
    elif rhyp > r1 and rhyp <= r2:
        logASEA = b1*log10(r1) + b2*log10(rhyp/r1) + c0
    elif rhyp > r2:
        logASEA = b1*log10(r1) + b2*log10(r2/r1) + b3*log10(rhyp/r2) + c0

    #A16 = logA - logASEA + c1 # constant H/V
    A16 = logA - logASEA + (c1l * log10(rhyp) + c2l) # linear H/V
    
    magstr = 'A16:\t' + str("%0.1f" % A16)
    if comp == 1: # if horizontal
        magstr = magstr + '*'
    print magstr
    
    return A16

'''    
# calc Yenier 2017
def calc_Y17(logA, rhyp, comp):
    from numpy import loadtxt, log10
    
    r1 = 100.
    r2 = 220.
    b1 = 1.399
    b2 = 0.727
    b3 = 1.806
    c1 = 0.102
    c2 = 4.354
    c3 = -1.579
    
    if rhyp <= r1:
        logY17 = b1*log10(rhyp) + 0.001*rhyp + c1
    elif rhyp > r1 and rhyp <= r2:
        logY17 = b2*log10(rhyp) + 0.001*rhyp + c2
    elif rhyp > r2:
        logY17 = b3*log10(rhyp) + 0.001*rhyp + c3

    Y17 = logA - logY17
    
    magstr = 'A17:\t' + str("%0.1f" % A16)
    if comp == 0: # if vertical
        magstr = magstr + '*'
    print magstr
    
    return Y17
'''    


def main(filename, logA, rhyp, eqdep):

    # get component orientation
    comp = get_component_orient(filename)

    # calculate epicentral distance
    repi = np.sqrt(rhyp**2 - eqdep**2)

    # header text
    print '\nML calculated from record: ' + filename \
        + ' at Rhyp = ' + str("%0.1f" % rhyp) + ' km'

    # get Richer (1935; 1958), extended by Eiby & Muir (1968)
    R35 = calc_R35(comp, logA, repi)

    # calculate Greenhalgh & Singh (1986)
    GS86 = calc_GS86(comp, logA, repi)

    # calculate Hutton & Boore (1987)
    HB87 = calc_HB87(comp, logA, rhyp)
    
    # calculate Bakun & Joyner (1984)
    BJ84 = calc_BJ84(comp, logA, rhyp)

    # calculate GG91
    GG91 = calc_GG91(comp, logA, rhyp)    
    
    # calculate MLM92
    MLM92 = calc_MLM92(comp, logA, rhyp)
    
    # calculate WGW96
    WGW96 = calc_WGW96(comp, logA, rhyp)
    
    # calculate A10
    A10 = calc_A10(comp, logA, rhyp)

    # add disclaimer
    print '\n* Does not use correct component'

    # now write ML data to file
    import write_data
    write_data.write_ML_dat(filename, rhyp, repi, logA, R35, HB87, MLM92, A10)  
    
    # export logA values for further analysis
    
    from datetime import datetime
    logAfile = 'log_A_values.csv'
    lines = open(logAfile).read()

    # append new line                      
    dt = datetime.strptime(filename.split('.')[0], '%Y%m%d%H%M')
    
    newline = ','.join((dt.strftime('%Y-%m-%d %H:%M:%S'), filename.split('.')[1], filename.split('.')[2], \
                        str('%0.1f' % rhyp), str('%0.1f' % repi), str('%0.4f' % logA)))
    
    f = open(logAfile, 'wb')
    f.write(lines + '\n' + newline)
    f.close()
    
    