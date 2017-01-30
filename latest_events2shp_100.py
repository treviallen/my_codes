#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
Created on Wed Jan 14 12:04:14 2015

Script to parse latest events page and export to shapefile

Usage:
    python latest_events2shp.py outshp
    
    where:
        outshp is output file name (e.g. latest_events.shp)

@author: tallen
"""

from lxml import html
import requests
import shapefile
from sys import argv


try:
    outshp = argv[1]
    
    # parse webpage
    page = requests.get('http://www.earthquakescanada.nrcan.gc.ca/recent/maps-cartes/index-eng.php')
    tree = html.fromstring(page.text)
    
    # get event data
    date = tree.xpath('//td[@headers="date"]/text()')
    time = tree.xpath('//td[@headers="time"]/text()')
    lat = tree.xpath('//td[@headers="lat"]/text()')
    lon = tree.xpath('//td[@headers="lon"]/text()')
    dep = tree.xpath('//td[@headers="depth"]/text()')
    mag = tree.xpath('//td[@headers="mag"]/text()')
    felt = tree.xpath('//td[@headers="felt"]/span[contains(@style, "color")]/text()')
    reg = tree.xpath('//td[@headers="region"]/text()')

    
    # convert relevant fields to float
    lat = [float(x) for x in lat]
    lon = [float(x) for x in lon]
    dep = [float(x)*1000 for x in dep]
    mag = [float(x) for x in mag]
    felt = [int(x.lower() == 'yes') for x in felt]
    
    # set shapefile headers
    w = shapefile.Writer(shapefile.POINT)
    w.field('DATE','C', 12)
    w.field('TIME','C', 10)
    w.field('LAT','F', 13, 6)
    w.field('LON','F', 13, 6)
    w.field('DEP','F', 13, 6)
    w.field('MAG','F', 13, 6)
    w.field('FELT','N', 4)
    w.field('REGION','C', 100)

    
    # now loop thru records and add to shapefile
    # the point has lon, lat, elevation=-depth, measure=mag
    for i in range(0, len(date)):
        w.point(lon[i], lat[i], -dep[i], mag[i])
        w.record(date[i], time[i], lat[i], lon[i], -dep[i], mag[i], felt[i], reg[i])
    
    print 'Writing shapefile...'
    w.save(outshp)
    
    # write projection file
    prjfile = outshp.strip().split('.')[0]+'.prj'
    f = open(prjfile, 'wb')
    f.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]')
    f.close()
    
except:
    print '\n python latest_events2shp.py outshp\n'
