# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 13:27:04 2013

Looks recursively in zone directories and finds all *.ps files

Takes a single argument, which is the model path

Requirements: ghostscript

Best used in DOS/Linux command line rather than iPython

Usage:

    python recursive_pdfjoin.py <model path>

e.g.:
    python recursive_pdfjoin.py 20130319_EA_H_Model_v3

@author: tallen
"""


from fnmatch import filter
from os import path, walk, system
from sys import argv

try:
    folder = argv[1]
    outpdf = folder+'.pdf'

    #matches = []
    psfiles = ''
    print '\n'+'Finding postscripts...'
    for root, dirnames, filenames in walk(folder):
        for filename in filter(filenames, '.pdf'):
            #matches.append(path.join(root, filename))
            psfiles = psfiles + ' ' + path.join(root, filename)

    # now call ghostscript to merge
    print '\n'+'Joining postscripts...'

    # initially assume DOS 32 bit
    try:
        system('gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER ' \
               + '-sOutputFile=' + outpdf +' ' + psfiles)

        print '\n'+'Postscripts joined successfully in ' + outpdf

    # now try unix/linux
    except:
        system('gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER ' \
               + '-sOutputFile=' + outpdf +' ' + psfiles)

        print '\n'+'Postscripts joined successfully in ' + outpdf

# Usage statement
except:
    print '\nUsage: python recursive_pdfjoin.py <model path>'
