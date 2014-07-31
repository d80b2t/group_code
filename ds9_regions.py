from astropy.io import fits as pf

class Region(object):
	def __init__(self,filepath=None):
		self.filepath = filepath
		self.datafile = pf.open(self.filepath)
#Put circles of constant radius around the good sources (flags=0)
	def circle(self):
		hdu=self.datafile[1].data
		file=open('Combined_ch1ch2_regions.reg','w') #Creates a file as a writable file
		file.write('# Region file format: DS9 version 4.1\n')
		file.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
		file.write('fk5\n')
		for i in range(len(hdu.X_IMAGE_1)): #Iterates over the x,y coordinates from Sextractor
 			if (hdu.FLAGS_1[i]==0) & (hdu.FLAGS_2[i]==0):
				x=hdu.ALPHA_J2000_1[i]
				y=hdu.DELTA_J2000_1[i]
				file.write('circle('+str(x)+','+str(y)+','+ '5"'+')\n') #Region files in DS9 are written in this format.

		file.close()
#Put regions around any of the flagged images (2,4,6...)
	def flags(self):
		hdu=self.datafile[1].data
		file=open('Combined_ch1ch2_regions.reg','a')
		for i in range(len(hdu.X_IMAGE_1)): #Iterates over the x,y coordinates from Sextractor:
			if (hdu.FLAGS_1[i]==2) or (hdu.FLAGS_2[i]==2):
				x=hdu.ALPHA_J2000_1[i]
				y=hdu.DELTA_J2000_1[i]
				file.write('circle('+str(x)+','+str(y)+','+ '5"'+')'+''+'#color=red width=3\n') 
			if (hdu.FLAGS_1[i]==4) or (hdu.FLAGS_2[i]==4):
				x=hdu.ALPHA_J2000_1[i]
				y=hdu.DELTA_J2000_1[i]
				file.write('circle('+str(x)+','+str(y)+','+ '5"'+')'+''+'#color=blue width=3\n')
		file.close()

#To use, call r=Region('File.fits'), r.circle(), r.flags() after importing this document (from ds9_regions import *)



