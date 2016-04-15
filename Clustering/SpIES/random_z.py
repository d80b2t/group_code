import os
import numpy as np
from astropy.io import fits as pf
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.interpolate as interpolate
from astroML.density_estimation import EmpiricalDistribution

dfile = '/Users/johntimlin/Clustering/SpIES/all_z/SDSS_QSO_inSpIES_3arc_20150619.fits' #Put the path to the data RA,DEC file here
rfile = '/Users/johntimlin/Clustering/SpIES/all_z/allz_rand_spies.txt' #Put the path to the random RA,DEC file here

#Read in the data
data=pf.open(dfile)[1].data
	

#Read in the random catalog
randat=open(rfile,'r') 
randheader=randat.readline()
randoms=randat.readlines()	
randra=np.zeros(len(randoms))
randdec=np.zeros(len(randoms))
#randz=np.zeros(len(randoms))

#From the Random file, pull out the RA and DEC values
for i in range(len(randoms)):
	randra[i]=float(randoms[i].split()[0])
	randdec[i]=float(randoms[i].split()[1])
	#randz[i]=float(randoms[i].split()[2])



gdx=(data.FLAG_2MASS_ch1==0) & (data.FLAG_2MASS_ch2==0) #Subset of SDSS quasars that do not fall within SpIES bright star radius 

factor=len(randra)/float(len(data.ZSPEC[gdx])) #The ratio of randoms to post cuts data

#Get bin edges for histograms
bins=np.linspace(min(data.ZSPEC[gdx]),max(data.ZSPEC[gdx]),60) #Creates steps from the smallest to largest redshift to use as bin edges for the data
binsmid = bins[:-1] + np.diff(bins)/2. #find the center of the bins

bin=np.linspace(min(data.ZSPEC[gdx]),max(data.ZSPEC[gdx]),len(bins))#Creates steps from the smallest to largest redshift to use as bin edges for the randoms
binmid = bin[:-1] + np.diff(bin)/2. #find the center of the bins

#The astroML.density_estimation.EmpiricalDistribution function takes a set of data, in this case redshifts, and Empirically learns the distribution of that data (must be 1-D). the rvs method then draws random variables from the input shape, so .rvs(shape). In our case, 'shape' is 1-D and is the length of the random data

randsz = EmpiricalDistribution(data.ZSPEC[gdx]).rvs(len(data.ZSPEC[gdx])*factor)

#Histogram both the data and random redshifts
dat,xd = np.histogram(data.ZSPEC[gdx], bins=bins)
rand,xr = np.histogram(randsz, bins=bin)

#Plot to make sure that the two distributions overlap
plt.figure(1)
plt.title('Redshift Distributions')
plt.plot(binsmid,dat,linestyle='steps-mid',label='data', linewidth=3)
plt.plot(binmid,rand/factor,linestyle='steps-mid',label='scaled randoms',linewidth=2)
plt.xlabel('z')
plt.ylabel('Number')
plt.legend()

#Make the random data into 3 arrays for RA, DEC, Z
RA=np.array(randra)
DEC=np.array(randdec)
Z=np.array(randsz)


#Write to Fits file
outfile = 'SpIES_random_allz.fits'
tbhdu=pf.BinTableHDU.from_columns([pf.Column(name='RA',format='D',array=RA),
pf.Column(name='DEC',format='D',array=DEC),
pf.Column(name='Z',format='D',array=Z)])
		
prihdr=pf.Header()
prihdr['COMMENT']="400,000 random points in the SpIES dual-band footprint"
prihdu=pf.PrimaryHDU(header=prihdr)

hdulist=pf.HDUList([prihdu,tbhdu])
hdulist.writeto(outfile)





plt.show()