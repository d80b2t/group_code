#LSS Clustering Code::
import numpy as np
from astropy.io import fits as pf
import astropy.units as u
from astropy.cosmology import angular_diameter_distance, comoving_distance, Planck13
from astropy.cosmology import FlatLambdaCDM as FLCDM
from astropy import coordinates as coord
from scipy import integrate, constants
import time

start=time.time()

#Import data 
obsdat=open('3yr_obj_v4_Sample8_nota01ao2s01_radecz.cat','r')
randat=open('rand_LRG_Sm8_spe.dat','r')
data=obsdat.readlines()
randoms=randat.readlines()
end1=time.time()
print end1-start, 'Files read'
#Create arrays of zeros which we use to fill with data
obsra=np.zeros(len(data))
obsdec=np.zeros(len(data))
obsz=np.zeros(len(data))
randra=np.zeros(len(randoms))
randdec=np.zeros(len(randoms))
randz=np.zeros(len(randoms))

for i in range(len(data)):
	obsra[i]=float(data[i].split()[0])
	obsdec[i]=float(data[i].split()[1])
	obsz[i]=float(data[i].split()[2])

for i in range(len(randoms)):
	randra[i]=float(randoms[i].split()[0])
	randdec[i]=float(randoms[i].split()[1])
	randz[i]=float(randoms[i].split()[2])


#Convert redshift to Mpc
#Create the 2SLAQ cosmology (H0=100h^-1, Om=0.3, Ol=0.7, Tcmb=2.725)
slaqcosmo=FLCDM(100,0.3,2.725)
obsX=slaqcosmo.comoving_distance(obsz)
randX=slaqcosmo.comoving_distance(randz)
end2=time.time()
print end2-start, 'Comoving distances calculated'



#Convert RA/DEC to cartesian coords
cutcoord=coord.SkyCoord(ra=obsra*u.degree,dec=obsdec*u.degree, distance=obsX,frame='icrs')
randcoord=coord.SkyCoord(ra=randra*u.degree,dec=randdec*u.degree, distance=randX,frame='icrs')
cx=cutcoord.cartesian.x
cy=cutcoord.cartesian.y
cz=cutcoord.cartesian.z
rx=randcoord.cartesian.x
ry=randcoord.cartesian.y
rz=randcoord.cartesian.z


tbhdu=pf.BinTableHDU.from_columns([pf.Column(name='datx',format='E',array=cx),
pf.Column(name='daty',format='E',array=cy),pf.Column(name='datz',format='E',array=cz),pf.Column(name='randx',format='E',array=rx),
pf.Column(name='randy',format='E',array=ry), pf.Column(name='randz',format='E',array=rz)])

prihdr=pf.Header()
prihdr['COMMENT']="Catalog of 100 square degree of random and observed galaxies in CMASS North"
prihdu=pf.PrimaryHDU(header=prihdr)

hdulist=pf.HDUList([prihdu,tbhdu])
hdulist.writeto('Cartesian_coords.fits')


