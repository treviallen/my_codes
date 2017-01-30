# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 13:37:32 2013

@author: tallen
"""

def shp2gmt(shpfile, outfile, **kwargs):
    '''
    kwargs decides what to to with the following header:
        header="none"
        header="num", field="field name"
        header="str", field="field name"
    '''

    import shapefile
     
    print 'Reading shapefile...'
    sf = shapefile.Reader(shpfile)
    shapes = sf.shapes()
    records = sf.records()
#    fields =sf.fields
    nrec = len(records)
    
    # now make output text file    
    # first overwirte
    f = open(outfile,'w')
    f.close()
    all_str = ''
    
    tmpfile = open(outfile,'a')
    dosimple = True
    
    # loop through polygons
    for k in range(0,nrec):
    # write polygon to temp file
    # if want to include name
    
    # if want to include quantity
    #all_str = '> -Z' + str(norm_complete[k,0]) + '\n'

        # check to see if shape has multiple parts        
        p = 0
        if len(shapes[k].points) != 0:
            parts = shapes[k].parts
            parts.append(len(shapes[k].points)-1)
            for part in range(0,len(parts)-1):
                pt_str = ''
                all_str = ''
                if dosimple == True:
                    all_str = '>' + '\n'
                    
                while p < parts[part+1]:
                    pt_str = pt_str + str("%0.5f" % shapes[k].points[p][0])+"\t" \
                         +str("%0.5f" % shapes[k].points[p][1])+"\n"
                    p += 1
                    
                all_str = all_str + pt_str
                tmpfile.write(all_str)
    
    print 'Writing to file...'    
    
    tmpfile.close()