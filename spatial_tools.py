# -*- coding: utf-8 -*-
"""
Created on Wed Aug 08 13:34:23 2012

@author: u56903
"""

import math

# this function gets earthquake-to-site distance (in km) on a sphere
'''
def distance(eqloc, stnloc):
    # Haversine formula example in Python
    lat1, lon1 = eqloc
    lat2, lon2 = stnloc
    radius = 6371. # km

    dlat = math.radians(lat2-lat1)
    dlon = math.radians(lon2-lon1)
    a = math.sin(dlat/2) * math.sin(dlat/2) + math.cos(math.radians(lat1)) \
        * math.cos(math.radians(lat2)) * math.sin(dlon/2) * math.sin(dlon/2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    d = radius * c

    return d
'''
# this function gets earthquake-to-site azimuth (in degrees) on a sphere
def azimuth(eqloc, stnloc):
    lat1 = math.radians(eqloc[0])
    lon1 = math.radians(eqloc[1])
    lat2 = math.radians(stnloc[0])
    lon2 = math.radians(stnloc[1])

    dlon = lon2-lon1

    y = math.sin(dlon) * math.cos(lat2)
    x = math.cos(lat1) * math.sin(lat2) \
        - math.sin(lat1) * math.cos(lat2) * math.cos(dlon)

    az = math.atan2(y, x)

    return math.degrees(az)

# this function writes new earthquake parameters to "eventlist.dat"
def write_evlist(evdate, eqla, eqlo, eqdep, eqmag):
    from os import getcwd
    datestr = evdate.strftime("%Y%m%d%H%M")
    if getcwd().startswith('/nas'):
        eqlist = '//nas//users//u56903//unix//eventlist.dat'
    else:
        eqlist = '/Users/trev/Documents/Earthquake_Data/eventlist.dat'
        
    newtxt = '\t'
    joinstr = (datestr, str(eqmag), str(eqlo), str(eqla), str(eqdep))
    newtxt = newtxt.join(joinstr)
    neweq = open(eqlist,'a')
    neweq.write('\n'+newtxt)
    neweq.close()

# this function asks user for earthquake parameters
def get_eq_params(evdate):
    
    try:
        # get earthquake latitude
        var = input('\n'+'Enter earthquake latitude (decimal degrees) > ')
        eqla = float(var)
        
        # get earthquake latitude
        var = input('\n'+'Enter earthquake longitude (decimal degrees) > ')
        eqlo = float(var)
        
        # get earthquake latitude
        var = input('\n'+'Enter earthquake depth (km) > ')
        eqdep = float(var)
        
        # get earthquake latitude
        var = input('\n'+'Enter earthquake magnitude (MW) > ')
        eqmag = float(var)
    
    except:
        # get earthquake latitude
        var = input('\n'+'Enter earthquake latitude (decimal degrees) > ')
        eqla = float(var)
        
        # get earthquake latitude
        var = input('\n'+'Enter earthquake longitude (decimal degrees) > ')
        eqlo = float(var)
        
        # get earthquake latitude
        var = input('\n'+'Enter earthquake depth (km) > ')
        eqdep = float(var)
        
        # get earthquake latitude
        var = input('\n'+'Enter earthquake magnitude (MW) > ')
        eqmag = float(var)

    # now, write new event info to event file
    write_evlist(evdate, eqla, eqlo, eqdep, eqmag)

    return eqla, eqlo, eqdep, eqmag

# function to ask for user input and get earthquake location,
# then calculate distance
def get_eq_distance(stlo, stla, evdate):
    import datetime as dt
    from os import getcwd
    from mapping_tools import distance

    # before asking for earthquake parameters, check event file
    if getcwd().startswith('/nas'):
        eqlist = '//nas//users//u56903//unix//eventlist.dat'
    else:
        eqlist = '/Users/trev/Documents/Earthquake_Data/eventlist.dat'

    # get max/min times to return (+/- 10 mins)
    mindate = evdate - dt.timedelta(minutes=60)
    maxdate = evdate + dt.timedelta(minutes=60)

    # read file and get events
    events = open(eqlist).readlines()
    i = 0
    tlist = []
    mlist = []
    lolist = []
    lalist = []
    dlist = []
    for line in events[1:]:
        eq = line.split('\t')
        if len(eq) == 5:
            tmpdate = dt.datetime.strptime(eq[0], "%Y%m%d%H%M")
            if tmpdate >= mindate and tmpdate <= maxdate:
                i += 1
                if i == 1:
                    print('\nSelect earthquake location from list...\n')
                    print('   ' + events[0].strip('\n'))
                print(str(i) + ') ' + line)
                # populate arrays
                tlist.append(dt.datetime.strptime(eq[0], "%Y%m%d%H%M"))
                mlist.append(float(eq[1]))
                lolist.append(float(eq[2]))
                lalist.append(float(eq[3]))
                dlist.append(float(eq[4]))

    # if only one earthquake returned, assume it is earthquake of interest
    # note: if two earthquakes with similar origin, eventlist.dat will need
    # to modified manually
    if i == 1:
        print('\nAutomatically associating earthquake parameters above...')
        print("Modify 'eventlist.dat' manually if this is not what you want!")

        evdate = tlist[0]
        eqmag = mlist[0]
        eqlo = lolist[0]
        eqla = lalist[0]
        eqdep = dlist[0]


    # select from multiple events
    elif i > 0:
        print(str(i+1) + ') None')
        var = input('\n'+'Enter selection > ')
        seleq = int(var)

        # ask for user input
        if seleq == i+1:
            eqla, eqlo, eqdep, eqmag = get_eq_params(evdate)

        elif seleq <= i:
            ind = i - 1
            evdate = tlist[ind]
            eqmag = mlist[ind]
            eqlo = lolist[ind]
            eqla = lalist[ind]
            eqdep = dlist[ind]

    # if no earthquakes in list
    else:
        eqla, eqlo, eqdep, eqmag = get_eq_params(evdate)

    # set variable for input to distance function
    repi = distance(eqla, eqlo, stla, stlo)[0]
    azim = distance(eqla, eqlo, stla, stlo)[1]
    rhyp = math.sqrt(repi**2 + eqdep**2)

#    az = azimuth(eqloc, stnloc)
#    print(az

    return rhyp, azim, eqla, eqlo, eqdep, eqmag, evdate
