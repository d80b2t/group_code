import numpy as np
from astropy.io import fits as pf
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM as FLCDM
#from astropy.cosmology import comoving_distance, Planck13
from astropy import coordinates as coord
from scipy import integrate, constants
import time

start=time.time()

dfile = '/Users/johntimlin/Clustering/SpIES/all_z/SDSS_QSO_inSpIES_3arc_20150619.fits' #Put the path to the data RA,DEC file here
rfile = '/Users/johntimlin/Clustering/SpIES/all_z/SpIES_random_allz.fits' #Put the path to the random RA,DEC file here


#Import RA,DEC,Z for data and randoms
obsdat=pf.open(dfile)[1].data
randat=pf.open(rfile)[1].data

end1=time.time()
print end1-start, 'Files read'

#Create arrays of zeros which we use to fill with data
gdx=(obsdat.FLAG_2MASS_ch1==0) & (obsdat.FLAG_2MASS_ch2==0)
obsra=obsdat.RA[gdx]
obsdec=obsdat.DEC[gdx]
obsz=obsdat.ZSPEC[gdx]
randra=randat.RA
randdec=randat.DEC
randz=randat.Z


#Convert redshift to Mpc
#Create the BOSS cosmology (H0=100h^-1, Om=0.274, Ol=0.726, Tcmb=2.725)
#This uses the astropy built in classes, the commented section below is my own integrator. 
slaqcosmo=FLCDM(100,0.274,2.725)
obsX=slaqcosmo.comoving_distance(obsz)
randX=slaqcosmo.comoving_distance(randz)
'''
#Integrator for comoving distance
def Comoving(x, Om, Og, Ol, Ok):
	F=3000*(x**-2.0)*(Om/x**3.0+Og/x**4.0+Ol+Ok/x**2.0)**-0.5
	return F

obsdist=[integrate.quad(Comoving, 1.0/(1.0+z), 1.0, args=(0.274,0.0,0.726,0.0)) for z in obsz] #quad(f(x),xlower,xupper, args)
randdist=[integrate.quad(Comoving, 1.0/(1.0+z), 1.0, args=(0.274,0.0,0.726,0.0)) for z in randz] #quad(f(x),xlower,xupper, args)
oX=np.array(obsdist)
rX=np.array(randdist)

obsX=oX[:,0]
randX=rX[:,0]

'''

end2=time.time()
print end2-start, 'Comoving distances calculated'



#Convert RA/DEC to cartesian coords using astropy units and skycoord
cutcoord=coord.SkyCoord(ra=obsra*u.degree,dec=obsdec*u.degree, distance=obsX,frame='icrs')
randcoord=coord.SkyCoord(ra=randra*u.degree,dec=randdec*u.degree, distance=randX,frame='icrs')
#The .cartesian.x converts the ra,dec, comoving coord to x (y or z)
cx=cutcoord.cartesian.x
cy=cutcoord.cartesian.y
cz=cutcoord.cartesian.z
rx=randcoord.cartesian.x
ry=randcoord.cartesian.y
rz=randcoord.cartesian.z


#Write the converted coordinates to a file (the structure of the file matters in the upcoming codes...

outfile = 'SpIES_allz_cartesian_coords_FLCDM.fits'

tbhdu=pf.BinTableHDU.from_columns([pf.Column(name='datx',format='D',array=cx),
pf.Column(name='daty',format='D',array=cy),pf.Column(name='datz',format='D',array=cz),pf.Column(name='randx',format='D',array=rx),
pf.Column(name='randy',format='D',array=ry), pf.Column(name='randz',format='D',array=rz)])

prihdr=pf.Header()
prihdr['COMMENT']="Catalog of all redshift quasars in the SpIES field"
prihdu=pf.PrimaryHDU(header=prihdr)

hdulist=pf.HDUList([prihdu,tbhdu])
hdulist.writeto(outfile)


