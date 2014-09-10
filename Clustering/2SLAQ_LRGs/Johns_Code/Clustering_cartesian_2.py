#LSS Clustering Code::
import numpy as np
from astropy.io import fits as pf
import astropy.units as u
from astropy.cosmology import angular_diameter_distance, comoving_distance, Planck13
from astropy import coordinates as coord
from scipy import integrate, constants
import time
from scipy.spatial.distance import pdist, cdist

#Import data here
obs=pf.open('Cartesian_coords.fits')[1].data
bdx=(obs.datx!=0)

#Combine the x,y,z data from the imported file in the form [[x1,y1,z1],[x2,y2,z2]...]
obsX=np.concatenate((np.array([obs.datx[bdx]]).T,np.array([obs.daty[bdx]]).T,np.array([obs.datz[bdx]]).T),axis=1)
randX=np.concatenate((np.array([obs.randx]).T,np.array([obs.randy]).T,np.array([obs.randz]).T),axis=1)

#rad=np.linspace(-0.8,-0.4,5)
#rad=np.linspace(-0.4,0.5,10)
rad=np.linspace(0.5,1.4,10)
#rad=np.linspace(1.4,2.3,10)
radii=10**rad
dd=[]
rr=[]
dr=[]
start=time.time()
#Separate the randoms into chunks of length 2000 to avoid memory errors...
steps=np.linspace(0,len(randX)-1,len(randX)/2000.)
print len(obsX)
print len(randX)
for r in range(len(radii)-1):
	start1=time.time()
	rad1=radii[r]
	rad2=radii[r+1]
	#Data Data pairs
	adx=(pdist(obsX)>=rad1) & (pdist(obsX)<=rad2)
	dd.append(len(pdist(obsX)[adx]))

	#Random Random pairs
	mat=[]
	RR=0
	DR=0
	tr=0
	#This loop is designed to split the data file into smaller chunks so that the pdist() function can be run without giving a memory error, which gives the distances between objects in the same chunk of the file
	for i in range(len(steps)-1): 
		#Set the values in the step array to be integers so that I can use them as indicies in separating the data into chunks	
		a=int(steps[i]) 
		b=int(steps[i+1])
		#This makes sure to not obs off the last value in the file which is why it is open ended.
		if i==len(steps)-2:
			val=randX[a:] #Finds the set of numbers starting (and including) with the integer 'a' through the last point
		else:	
			val=randX[a:b] #Finds the set of numbers starting (and including) with the integer 'a' and ending (not including) integer 'b'
		mat.append(val) #Append these to a list for use in the cdist loops
		gdx=(pdist(val)>=rad1) & (pdist(val)<=rad2) #Set the pdist values less than or equal to some r value (to be determined later...) for DD, RR, DR pairs
		RR+=len(pdist(val)[gdx]) #additive number of points in the random random paircounts (can add radius info here with gdx)
		tr+=len(pdist(val))

	#This loop finds the distances between the objects in one chunk of data and another chunk of data, so two arrays need to be input here. Mat is the list of arrays from above loop
	for i in range(len(mat)):
		for j in range(len(mat)):
			if i<j: #makes sure I do not double count the distances between points
				ddx=(cdist(mat[i],mat[j])>=rad1) & (cdist(mat[i],mat[j])<=rad2)	#Calculates the distances between the elements in mat[i] with those in mat[j]
				RR+=np.size(cdist(mat[i],mat[j])[ddx]) #Output is 2D array, size gives total number of elements. (Can add r with ddx)
				tr+=np.size(cdist(mat[i],mat[j]))
		#Data Random pairing
		cdx=(cdist(obsX,mat[i])>=rad1) & (cdist(obsX,mat[i])<=rad2)
		DR+=np.size(cdist(obsX,mat[i])[cdx]) #DR pairing for these bins

	rr.append(RR)
	dr.append(DR)
	end=time.time()
	print 'Loop', r+1,'of', len(radii)-1,'Done in',end-start1
	print tr

print 'Done in', end-start,'s'
print 'Data-Data pairs', dd
print 'Rand-Rand pairs', rr
print 'Data-Rand pairs', dr
print 'Separations', radii
sep=[]
for i in range(len(radii)-1):
	a=radii[i]
	b=radii[i+1]
	sep.append((a+b)/2.)
print sep
linsep=[]
for i in range(len(rad)-1):
	c=rad[i]
	d=rad[i+1]
	linsep.append((c+d)/2.)
print 'linear separations', linsep

file=open('Pair_counts_2slaq_small.txt','a')
#file.write('S DD RR DR\n')
for i in range(len(dd)):
	S=sep[i]
        dat=dd[i]
        rand=rr[i]
        drand=dr[i]
        file.write(str(S)+'\t'+str(dat)+'\t'+str(rand)+'\t'+str(drand)+'\n')

file.close()










