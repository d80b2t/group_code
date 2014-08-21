#LSS Clustering Code::
import numpy as np
from astropy.io import fits as pf
import astropy.units as u
from astropy.cosmology import angular_diameter_distance, comoving_distance, Planck13
from astropy import coordinates as coord
from scipy import integrate, constants
import time
from scipy.spatial.distance import pdist, cdist
from sklearn.metrics.pairwise import pairwise_distances
#Import data here
#galdat=pf.open('galaxy_DR10v8_CMASS_North.fits')[1].data
#galhead=pf.open('galaxy_DR10v8_CMASS_North.fits')[1].header

start=time.time()

#Import data region cut (100 sq deg)
cutdat=pf.open('galaxy_DR10v8_CMASS_North_100sqdegsq.fits')[1].data
randat=pf.open('random1_DR10v8_CMASS_North_100sqdegsq.fits')[1].data
end1=time.time()
print end1-start, 'Files read'
#Find the RA, DEC, Z
#galra=galdat.RA[:100]
#galdec=galdat.DEC[:100]
#galz=galdat.Z[:100]

cutra=cutdat.RA
cutdec=cutdat.DEC
cutz=cutdat.Z
randra=randat.RA
randdec=randat.DEC
randz=randat.Z


#Convert redshift to Mpc (Using Planck Results stored in astropy.cosmology)
#comoving=Planck13.comoving_distance(galz)
cutX=Planck13.comoving_distance(cutz)
randX=Planck13.comoving_distance(randz)
end2=time.time()
print end2-start, 'Comoving distances calculated'

'''
#From working cosmology code
def Horizon(x, t, Om, Og, Ol, Ok):
	F=(x**-2)*(Om/x**3+Og/x**4+Ol+Ok/x**2)**-0.5
	return F
Om=0.307
Og=0
Ol=0.693
Ok=1-(Om+Og+Ol)
t=linspace(0,1,100000)
DA=[]
chi=[]
for i in galz:
	com=integrate.quad(Horizon,1/(1+i),1.0,args=(t,Om,Og,Ol,Ok))
	X=3e5/67.8 * com[0]
	DA.append(X/(1+i))
	chi.append(X)
#for i in range(len(DA)):
#	print DA[i]*u.Mpc-dist[i]
print chi[0], chi[1]
print comoving[0], comoving[1]
'''

#Convert RA/DEC to cartesian coords
#c=coord.SkyCoord(ra=galra*u.degree,dec=galdec*u.degree, distance=comoving,frame='icrs')
cutcoord=coord.SkyCoord(ra=cutra*u.degree,dec=cutdec*u.degree, distance=cutX,frame='icrs')
randcoord=coord.SkyCoord(ra=randra*u.degree,dec=randdec*u.degree, distance=randX,frame='icrs')
cx=cutcoord.cartesian.x
cy=cutcoord.cartesian.y
cz=cutcoord.cartesian.z
rx=randcoord.cartesian.x
ry=randcoord.cartesian.y
rz=randcoord.cartesian.z

end3=time.time()
print end3-start, 'Cartesian Coordinates found'

cutX=np.concatenate((np.array([cx]).T,np.array([cy]).T,np.array([cz]).T),axis=1)
print 'cutX'
randX=np.concatenate((np.array([rx]).T,np.array([ry]).T,np.array([rz]).T),axis=1)
print 'randX'
#print pdist(cutX)
print pairwise_distances(randX)
end1=time.time()
print end1-start
'''
#Test the cartesian coordinate converter (ra and dec must be specified in degrees because np.cos() must be in radians.)
r=comoving
theta=galdec*u.degree
phi=galra*u.degree

x=r*np.cos(theta)*np.cos(phi)
y=r*np.cos(theta)*np.sin(phi)
z=r*np.sin(theta)

print x[0], c.cartesian.x[0]
print y[0], c.cartesian.y[0]
print z[0], c.cartesian.z[0]
print x[1], c.cartesian.x[1]
print y[1], c.cartesian.y[1]
print z[1], c.cartesian.z[1]

sep=np.sqrt((x[0]-x[1])**2+(y[0]-y[1])**2+(z[0]-z[1])**2)

print sep, c[0].separation_3d(c[1])
'''
'''
#Create random points using the distribution of LRG's in this data set x50

#Find DD pairs
DD=[]
for i in range(len(cutcoord)):
	for j in range(len(cutcoord)):
		if j>i:
			DD.append(cutcoord[i].separation_3d(cutcoord[j]))



#Find RR pairs
RR=[]
for i in range(len(randcoord)):
	for j in range(len(randcoord)):
		if j>i:
			RR.append(randcoord[i].separation_3d(randcoord[j])))
print len(RR)

#Find DR pairs
DR=[]
for i in range(len(c)):
	for j in range(len(c)):
		if i+j<len(c) and i+j!=i:
			s.append(c[i].separation_3d(c[i+j]))
print len(s)
'''
#Construct Landay-Szalay correlation function


