def mergeseed(mseedfiles, outfile):

    '''
    mseedfiles = tuple of filenames
    '''
    
    from obspy.core import read
    
    #mseedfiles = ('2012171_105400_0a903_2_1.seed', '2012171_105600_0a903_2_1.seed')
    #outfile = '201210191054.GEES.EHE.seed'
    
    st = read(mseedfiles[0])
    
    if len(mseedfiles) > 1:
        for i in range(1,len(mseedfiles)):
            st += read(mseedfiles[i])
            
    # now write to file
    st.write(outfile, format='MSEED')