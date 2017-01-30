def parse_zone_shp(shpfile):
    import shapefile
    from os import path

    print 'Reading shapefile...'
    sf = shapefile.Reader(shpfile)
    shapes = sf.shapes()
    records = sf.records()
    fields = sf.fields
    fields = fields[1:]

    # get model path from shapefile
    shp_path = path.basename(shpfile).split('.')
    shp_path = shp_path[0]

    # get zone name and code field number
    for i in range(0,len(fields)):
        if fields[i][0] == 'NAME':
            iname = i

        if fields[i][0] == 'CODE':
            icode = i

    # now loop through records and make input files
    k = 0
    zondict = []
    for rec in records:
        tmpdict = {}
        print 'Making inputs for source zone: ' + rec[iname]

        tmpdict['code'] = rec[icode]
        tmpdict['name'] = rec[iname]
        tmpdict['poly'] = shapes[k].points

        zondict.append(tmpdict)
        k+=1

    return zondict


def write_eqrm_szone(hazparams, zondict, outxml, **kwargs):
    import shapefile

    # test data
    bval = 1.
    aval = 1.
    area = 10000.
    rec_min_mag = 3.3
    rec_max_mag = 5.4
    gen_min_mag = 4.5
    mag_sample_bins = 15
    nevents = 5000

    dip = 35
    delta_dip = 0
    azimuth = 180
    delta_azimuth = 180
    depth_top_seismogenic = 7
    depth_bottom_seismogenic = 15.60364655

    # start writing xml
    xmltxt = '<source_model_zone magnitude_type="Mw">\n'
    for zone in zondict:
        xmltxt += "    <zone area = '"+str(area)+"' event_type = 'crustal fault'>\n"
        xmltxt += "        <geometry\n" \
                + "            dip = '"+str(dip)+"'\n" \
                + "            delta_dip = '"+str(delta_dip)+"\n" \
                + "            azimuth = '"+str(azimuth)+"'\n" \
                + "            delta_azimuth = '"+str(delta_azimuth)+"'\n" \
                + "            depth_top_seismogenic = '"+str(depth_top_seismogenic)+"'\n" \
                + "            depth_bottom_seismogenic = '"+str(depth_bottom_seismogenic)+"'>\n" \
                + "            <boundary>\n"

        # now get zone polygon data
        for i in range(0,len(zone['poly'])):
            xmltxt += "                "+str("%0.3f" % zone['poly'][i][1])+"  " \
                                        +str("%0.3f" % zone['poly'][i][0])+"\n"

        xmltxt += "            </boundary>\n" \
                + "        </geometry>\n" \
                + "        <recurrence_model\n" \
                + "            distribution = 'bounded_gutenberg_richter'\n" \
                + "            recurrence_min_mag = '"+str(rec_min_mag)+"'\n" \
                + "            recurrence_max_mag = '"+str(rec_max_mag)+"'\n" \
                + "            recurrence_max_mag = '"+str(rec_max_mag)+"'\n" \
                + "            A_min = '"+str(aval)+"'\n" \
                + "            b = '"+str(bval)+"'>\n" \
                + "            <event_generation\n" \
                + "                generation_min_mag = '"+str(gen_min_mag)+"'\n" \
                + "                number_of_mag_sample_bins = '"+str(mag_sample_bins)+"'\n" \
                + "                number_of_events = '"+str(nevents)+"'/>\n" \
                + "        </recurrence_model>\n" \
                + "    </zone>\n" \

    xmltxt += "</source_model_zone>"

    return xmltxt