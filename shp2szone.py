# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 13:37:32 2013

Note: Now use write_areal_source.write_cmpfile and write_areal_source.write_betainp
to generate *.cmp and beta.pl files from input tables

@author: tallen
"""

def make_szone_inputs(folder, zone, catalogue, **kwargs):
    from os import path
    for key in kwargs:
        if key == 'head':
            header = kwargs[key]
    # check if header is defined
    try:
        header
    except NameError:
        header = None
    
    # if needed, get header from .zon file
    if header == None:
        zonfile = path.join(folder,zone + '.zon')
        zontxt = open(zonfile).readlines()
        header = zontxt[0]

    
    '''
    now make input file for szonegmt
    '''
    # note - assumes catalogue is in dir above zone directories
    # potentially use path.abspath() if ever have problems with this!    
    szonetxt = zone + '\n' + zone + '.zon' + '\n' \
             + zone + '.cmp' + '\n' + path.join('..','..',catalogue) + '\n'
        
    # write szonegmt.inp
    outpath = path.join(folder, 'szonegmt.inp')
    f = open(outpath,'wb')
    f.write(szonetxt)
    f.close()
    

def shp2szone(shpfile, **kwargs):
    '''
    kwargs decides what to to with the following header:
        header="none"
        header="num", field="field name"
        header="str", field="field name"
    '''

    import shapefile
    from os import mkdir, path

    print 'Reading shapefile...'
    sf = shapefile.Reader(shpfile)
    shapes = sf.shapes()
    records = sf.records()
    fields = sf.fields
    fields = fields[1:]
    
    # search for region code and name
    for i in range(0,len(fields)):
#        print fields
        if fields[i][0] == 'CODE':
            icode = i
            
        if fields[i][0] == 'NAME':
            iname = i
#    icode = 1
#    iname = 2
    
    # get model path from shapefile
    shp_path = path.basename(shpfile).split('.')
    shp_path = shp_path[0]
    
    # now loop through records and make input files
    k = 0
    print 'Number of zones: ',len(records)
    for rec in records:
#        print len(rec), rec
        reg_code = rec[icode]
        reg_name = rec[iname]
        print 'Making inputs for source zone: ' + reg_code
        
        # set output filename
        
        path_str = reg_code
        zonefile = reg_code + '.zon'
        outpath = path.join(shp_path,path_str,zonefile)
        
        # make header
        header = reg_code + ' - ' + reg_name + '\n'
        footer = '99        99' + '\n'
        
        # now get coords
        pt_str = ''
        for kk in range(len(shapes[k].points)):
            # output format is: lat lon
            pt_str = pt_str + str("%0.3f" % shapes[k].points[kk][1])+"\t" \
                     +str("%0.3f" % shapes[k].points[kk][0])+"\n"
        
        all_str = header + pt_str + footer
        
        # now check to make sure output path exists
        try:
            f = open(outpath,'wb')
        except:
            if path.exists(shp_path) == False:
                mkdir(shp_path)
            
            mkdir(path.join(shp_path,path_str))
            f = open(outpath,'wb')
        
        # write zone file to file
        f.write(all_str)
        f.close()
        
        '''
        This functionality moved to make_cmp_betapl_inputs.py
        '''
        # now make all zone inputs
#        folder = path.join(shp_path,path_str)
#        make_szone_inputs(folder, reg_code, catalogue, head=header)
                
        k += 1

# using this assumes polygon files already built        
def shp2szone_simple(shpfile, catalogue):
    from os import listdir, path
    
    # get model path from shapefile
    shp_path = path.basename(shpfile).split('.')
    shp_path = shp_path[0]

    dirs = listdir(shp_path)
    for zone in dirs:
        folder = path.join(shp_path,zone)
        try:
            print 'Making inputs for source zone: ' + zone
            # now make all zone inputs
            make_szone_inputs(folder, zone, catalogue) 
                    
        except:
            print zone + " is not a folder"            