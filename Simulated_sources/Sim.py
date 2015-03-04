import sys
import os
import numpy as np
from astropy.io import fits as pf
from astropy import wcs
from pylab import *

################## DEFINE THE 2D GAUSSIAN #######################################
def Gaussian2d(x,y,amplitude,xc,yc,sigmax,sigmay):
	x=amplitude*np.exp(-((x-xc)**2/(2.*sigmax**2)+(y-yc)**2/(2.*sigmay**2)))
	return x



###################################### DEFINE THE SIMULATES SOURCES CODE ######################################
def simulated(zeromat,dat,header):
	fact=1.0*10**12 * 1./(180/np.pi)**2 * 1.0/(3600.)**2 * (0.6)**2
	#Create random pixel values (x,y) around which we'll build the gaussian
	yc=np.random.randint(0,np.shape(dat)[0]-1)
	xc=np.random.randint(0,np.shape(dat)[1]-1)
	while np.isnan(dat[yc,xc])==True:
		yc=np.random.randint(0,np.shape(dat)[0]-1)
		xc=np.random.randint(0,np.shape(dat)[1]-1)
	
	#Set random amplitude for the gaussian (in magnitudes)
	mag=np.random.uniform(17,28)
	#Convert magnitude to flux(uJy) then to flux (MJy/sr) divide by 4 and make the source 4 pixels in range
	F=10**((mag-23.9)/-2.5)
	f=F/fact

	
	#Calculate the sigma using FWHM=2*root(2ln(2))sigma with Spitzer FWHM=1.7" and pix scale=0.6"/pix
	sigmax=sigmay=1.7/(2*np.sqrt(2*log(2)))/0.6
	
	#Create the Gaussian Matrix
	r = 10*1.7 #(5 * FWHM of Spitzer in radius)
	G=0 
	for i in range(int(round(r))):
		for j in range(int(round(r))):
			x=xc-i
			y=yc-j
			g=Gaussian2d(x,y,f,xc-2,yc-2,sigmax,sigmay)
			G+=g
			zeromat[y,x]=g
			
	#Convert to magnitude using the summed pixel values from each iteration of the gaussian
	magreal= -2.5*np.log10(G*fact)+23.9
	
	#Find RA and DEC in J2000 coordinates using the header from the original mosaic
	w=wcs.WCS(header[0].header)
	pixcoord=[[xc,yc]]
	coords=np.array(pixcoord)
	world=w.wcs_pix2world(coords,0)  
	ra=world[0][0]
	dec=world[0][1]
	
	############################ Write to ASCII file #########################################################
	############### CAUTION:: Will append to any file with the name below so be sure that the name of the file in the open() statement is correct!! #################################
	infile=open('position.ascii','a')
	infile.write(str(ra)+'\t'+str(dec)+'\t'+str(magreal)+'\t'+str(G)+'\t'+str(xc)+'\t'+str(yc)+'\n')
	infile.close()
	info=[ra,dec,magreal,G]
	print info
	return zeromat


################### IMPLEMENT THE CODE HERE #######################################################################
###### BE SURE TO CHANGE THE NAME OF THE INPUT FILE TO MATCH THE ONE YOU WANT TO RUN #########
dat= pf.open('aor3_ch1_mosaic.fits')[0].data
header=pf.open('aor3_ch1_mosaic.fits')
sim=np.zeros(np.shape(dat))
for i in range(4000):
	simulated(sim,dat,header)
	

new=dat+sim
hdulist=pf.PrimaryHDU(data=new,header=header[0].header)
hdulist.writeto('aor3_ch1_mosaic_simsource.fits')





##################PYTHON TEST PLOTTING###################################
'''
imshow(sim)
colorbar()
show()
'''

######### If we dont need gaussian profiles, we can use this tophat profile as well ######################################


'''
	#For Tophat sources (non-gaussian)
	f=F/fact/9
	#add the random source catalog to the original mosaic and output fits file
		
	sim[x,y]=f
	sim[x+1,y]=f
	sim[x+1,y+1]=f
	sim[x,y+1]=f
	sim[x-1,y]=f
	sim[x-1,y-1]=f
	sim[x,y-1]=f
	sim[x+1,y-1]=f
	sim[x-1,y+1]=f
'''
