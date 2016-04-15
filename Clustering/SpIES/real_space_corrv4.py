import os
from time import time
import astropy.io.fits as pf
import numpy as np
from numpy.linalg import norm
from astroML.correlation import two_point, bootstrap_two_point
from multiprocessing import Pool
from functools import partial


#Open the cartesian coordinates used to compute the redshift space correlation function

start1=time()
#Path to the data on my machine
file='/home/jtimlin/CLUSTERING_IN_SPIES/all_z/SpIES_allz_cartesian_coords_FLCDM.fits'
data=pf.open(file)[1].data

# Separate the data and random points and stack so that I have an array of [x,y,z] coordinates for each quasar

detected = data.datx > 0

datx = data.datx[detected]
daty = data.daty[detected]
datz = data.datz[detected]
datvect = np.concatenate((np.array([datx]).T,np.array([daty]).T,np.array([datz]).T),axis = 1)

randx = data.randx
randy = data.randy
randz = data.randz
randvect = np.concatenate((np.array([randx]).T,np.array([randy]).T,np.array([randz]).T),axis = 1)


# Find the separation angle between quasars
# Calculate pi = |s1-s2|cos(theta/2). For small angles (<2 degrees) cos(theta/2)=1
# Calculate r_p = (s1+s2) * sin(theta/2). For small angles (<2 degrees) sin(theta/2) = theta/2

#Define a function to calculate rp and pi for the data
def data_loop(N,rpradii,piradii):
	DD = np.zeros((len(rpradii)-1,len(piradii)-1))
	for i in range(N):
		if i % 1000==0:
			print i/float(len(datvect))  #Step counter
		vect = datvect[i]  #The x,y,z cartesian vector of a quasar in the data wrt the origin [x1,y1,z1]
		rest = datvect[i+1:]  #The remaining array of [x,y,z] vectors in the data set (i.e.[[x2,y2,z2],[x3,y3,z3],...
		dprod = [np.dot(vect,j) for j in rest]  #The dot product of the first vector and the remaining vectors
		angle = np.arccos(np.divide(dprod, (norm(vect) * norm(rest, axis = 1)))) # The angle between the first vector and the remaining vectors
		gdx = angle*180./np.pi <2 #Look for correlations on small scales
		p = np.abs(norm(vect) - norm(rest[gdx],axis=1))  #Calculate pi akin to Peebles 1980 text book for small angles
		r = np.abs(norm(vect) + norm(rest[gdx], axis=1)) * angle[gdx]/2.0 #Calculate rp akin to Peebles 1980 text book for small angles
		#Build a DD matrix with increasing rp across the matrix and increasing pi down the matrix
		for k in range(len(piradii)-1):
			for j in range(len(rpradii)-1):
				gdat = (p > piradii[k]) & (p<=piradii[k+1]) & (r>rpradii[j]) & (r<=rpradii[j+1]) #Cut data in pi + rp space
				DD[k,j] += float(len(p[gdat]))

	return DD
	
#Define function to calculate rp and pi for randoms (See function above for comments)
def rand_loop(rpradii,piradii,count):
	RR = np.zeros((len(rpradii)-1,len(piradii)-1))
	for i in count:
		s = time()
		vect = randvect[i]
		rest = randvect[i+1:]
		dprod = np.array([np.dot(vect,j) for j in rest])
		angle = np.arccos(np.divide(dprod,(norm(vect) * norm(rest, axis = 1))))
		gdx = angle*180./np.pi < 2
		#Calculate all pi and rp values for small angles
		p = np.abs(norm(vect) - norm(rest[gdx],axis=1))
		r = np.abs(norm(vect) + norm(rest[gdx], axis=1)) * angle[gdx]/2.0
		#Build an RR matrix with increasing rp across the matrix and increasing pi down the matrix
		for k in range(len(piradii)-1):
			for j in range(len(rpradii)-1):
				grand = (p > piradii[k]) & (p<=piradii[k+1]) & (r>rpradii[j]) & (r<=rpradii[j+1]) #Cut data in pi + rp space
				RR[k,j] += float(len(p[grand]))
		e = time()
		print "RR loop",i,"done in",e-s,"s"
	return RR
	

	
#Define function to calculate rp and pi for the data - random pairs (See function above for comments)
def dat_rand_loop(rpradii,piradii,count):
	DR = np.zeros((len(rpradii)-1,len(piradii)-1))
	for i in count:
		s = time()
		vect = datvect[i]
		rest = randvect
		dprod = np.array([np.dot(vect,j) for j in rest])
		angle = np.arccos(np.divide(dprod,(norm(vect) * norm(rest, axis = 1))))
		gdx = angle*180./np.pi < 2
		p = np.abs(norm(vect) - norm(rest[gdx],axis=1))
		r = np.abs(norm(vect) + norm(rest[gdx], axis=1)) * angle[gdx]/2.0
		#Build an RR matrix with increasing rp across the matrix and increasing pi down the matrix
		for k in range(len(piradii)-1):
			for j in range(len(rpradii)-1):
				gdrand = (p > piradii[k]) & (p<=piradii[k+1]) & (r>rpradii[j]) & (r<=rpradii[j+1]) #Cut data in pi + rp space
				DR[k,j] += float(len(p[gdrand]))
		e = time()
		print "DR loop",i,"done in",e-s,"s"
	return DR

##############################################################
#Make arrays of Rp and Pi

RP=np.linspace(0,1.75,20)
rpradiiinit=10**RP
rpradii=np.insert(rpradiiinit,0,-1)
PI=np.linspace(0,1.75,20)
piradiiinit=10**PI
piradii = np.insert(piradiiinit,0,-1)

#Run the data conversion loop
data = data_loop(len(datvect),rpradii,piradii)

DDpairs = data #Returned DD matrix

end1= time()

print "Data conversion done in",end1-start1,"s"

start2=time()

########################################################################
#Begin multiprocess of randoms

nproc = 10  #Number of processors on which to run
vals = np.arange(len(randvect))  #an array of integers the length of the data to split for multiprocessing
ev=len(vals)/nproc #Define an equal split of data...NOTE!!!! THIS MUST BE INT DIVISION FOR THIS TO WORK!
sploc=[]
for i in range(nproc-1):
	sploc.append((i+1)*ev)
ct = np.split(vals,sploc)

#ct = np.split(vals,nproc) #split the vals array into nproc to run a subset of the data on each processor

#Call the Pool command to multiprocess
if __name__ == '__main__':
	p = Pool(nproc)
	func = partial(rand_loop,rpradii,piradii)
	randoms = p.map(func, ct) #Map the random_loop routine to each processor with a single array in ct to run over
	RRdummy = np.zeros((len(rpradii)-1,len(piradii)-1))
	for i in range(len(randoms)):
		RRdummy = np.add(RRdummy,randoms[i])
	RRpairs = RRdummy

end2=time()
print "Random conversion done in",end2-start2,"s"





##################################################
#Begin multiprocess of data-randoms
start3=time()
nproc2 = 10  #Number of processors on which to run
vals2 = np.arange(len(datvect))  #an array of integers the length of the data to split for multiprocessing
c=len(vals2)/nproc2 #Define an equal split of data...NOTE!!!! THIS MUST BE INT DIVISION FOR THIS TO WORK!
array=[]
for i in range(nproc2-1):
	array.append((i+1)*c)
ct2 = np.split(vals2,array)


#ct2 = np.split(vals2,nproc2) #split the vals array into nproc to run a subset of the data on each processor

#Call the Pool command to multiprocess
if __name__ == '__main__':
	p = Pool(nproc2)
	funcdr = partial(dat_rand_loop,rpradii,piradii)
	datrands = p.map(funcdr, ct2) #Map the data-random_loop routine to each processor with a single array in ct to run over
	DRdummy = np.zeros((len(rpradii)-1,len(piradii)-1))
	for i in range(len(datrands)):
		DRdummy = np.add(DRdummy,datrands[i])
	DRpairs = DRdummy
	
end3=time()

print "Data Random conversion done in",end3-start3,"s"
start4=time()
#Calculate XI matrix 

factor = len(randvect) * 1.0/len(datvect)
xi = (pow(factor,2)*DDpairs - (2 * factor * DRpairs) + RRpairs) / RRpairs

for i in range(len(piradii)-1):
	for j in range(len(rpradii)-1):
	
		DD = DDpairs[i,j]
		DR = DRpairs[i,j]
		RR = RRpairs[i,j]
		XI = xi[i,j]
		#Write/Append to file depending on whether file exists or not
		if os.path.exists("SpIES_RealSpace_Clustering_v8.txt"):
			infile=open("SpIES_RealSpace_Clustering_v8.txt" ,"a")
			PIavg1 = (piradii[i] + piradii[i+1])/2.0
			RPavg1 = (rpradii[j] + rpradii[j+1])/2.0
			infile.write(str(PIavg1) + ' '+ str(RPavg1)+ ' '+ str(DD)+ ' ' + str(DR) + ' ' + str(RR) + ' '+str(XI)+'\n')
			infile.close()

		else:
			infile=open("SpIES_RealSpace_Clustering_v8.txt" ,"w")
			infile.write("PIavg RPavg DD DR RR XI \n")
			PIavg1 = (piradii[i] + piradii[i+1])/2.0
			RPavg1 = (rpradii[j] + rpradii[j+1])/2.0
			infile.write(str(PIavg1) + ' '+ str(RPavg1)+ ' '+ str(DD)+ ' ' + str(DR) + ' ' + str(RR) + ' '+str(XI)+'\n')
			infile.close()
    		

end4=time()

print 'Full 2D 2PCF calculated in',end4-start4,"s"


print "Total Time=",end4-start1,"s"
