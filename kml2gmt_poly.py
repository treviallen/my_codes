#*******************************************************************************
#Program: kml2gmt_poly.py
#
#Looks for Google earth kml files of polygons in any target directory and
#converts them to a single gmt polygon file.
#
#Useage: python kml2gmt_poly.py <dir> <outfile>
#*******************************************************************************

import os
import sys

kmldir  = sys.argv[1]
outfile = sys.argv[2]

# look for kml files
dirpath = os.path.join(kmldir)
kml_list=os.listdir(dirpath)

kmltxt = ''
coord_txt = []

for fname in kml_list:

	if fname.find('.kml')>0:
		print 'Reading '+fname
		fpath = os.path.join(kmldir,fname)
		kmlread = open(fpath).readlines()

		# get coordinates
		coord_index = -1
		for line in kmlread:
			if coord_index > 0:
				line = line.replace(',0 ',',0,')
				coord_txt = line.lstrip().rstrip().strip("'").split(',')
				coord_index = -1

			start_index = line.find('<coordinates>')

			if start_index >= 0:
				coord_index = 1

		txt = ''
		k = 0
		for i in range(len(coord_txt)):
			k = k + 1
			if k == 1 and i < len(coord_txt)-1:
				num1 = str(coord_txt[i])
				num2 = str(coord_txt[i+1])
				txt = txt + num1 + ',' + num2 +'\n'
			if k == 3:
				k = 0
	kmltxt = kmltxt + '> ' + fname + '\n' + txt

gmt_file = open(outfile,'w')
gmt_file.write(kmltxt)
gmt_file.close()
