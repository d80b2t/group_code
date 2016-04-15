from astropy.io import fits as pf
from astroML.correlation import two_point, bootstrap_two_point
import numpy as np
from matplotlib import *
from pylab import *
from time import time

#import the cartesian file here
dfile = '/Users/johntimlin/Clustering/SpIES/all_z/SpIES_allz_cartesian_coords_FLCDM.fits' #Put the path to the data RA,DEC file here

start=time()
#Import data here
obs=pf.open(dfile)[1].data
bdx=(obs.datx!=0) #This file cuts out the extra rows in the data column which were put in to match the number of randoms (fits tables have to have the same number of rows to be saved properly) 

############################################################################################
#COMPUTE THE CORRELATION FUNCTION
#Combine the x,y,z data from the imported file in the form [[x1,y1,z1],[x2,y2,z2]...]
obsX=np.concatenate((np.array([obs.datx[bdx]]).T,np.array([obs.daty[bdx]]).T,np.array([obs.datz[bdx]]).T),axis=1)
randX=np.concatenate((np.array([obs.randx]).T,np.array([obs.randy]).T,np.array([obs.randz]).T),axis=1)

print 'Data Read in'

#Define the radii of the annuli over which to run the analysis (same as Ross 2009 radii)
rad=np.arange(0.1,4.1,0.1)
radii=10**rad
#Use the astroML two_point function to compute the Landy-Szalay estimator of the two point correlation function. Since we have random data, import that in the data_R section and turn off random_state (None)
#This function will create a tree of points and compute pairwise distances between points. It then counts how many are within the radii
#I had to adjust the astroML code to return the DD, DR and RR pairs along with XI for use in the errors later

TPCF=two_point(obsX,radii,method='landy-szalay',data_R=randX,random_state=None)

xifull=np.asarray(TPCF[0])
DDfull=np.asarray(TPCF[1])
DRfull=np.asarray(TPCF[2])
RRfull=np.asarray(TPCF[3])
print 'Full 2PCF calculated'
end = time()

print end-start,'s'

#############################################################
#calculate Poisson Errors using Equation 13 in Ross 2009
#CALCULATE POISSON ERRORS USING EQUATION 13 IN ROSS 2009
Poisson_err = (1+xifull)/np.sqrt(1.0*DDfull)

print 'Poisson Errors calculated'

########################################################
#CALCULATE THE CENTERS OF THE BINS
sep=[]
for i in range(len(radii)-1):
	a=radii[i]
	b=radii[i+1]
	sep.append((a+b)/2.)
print 'Separations Calculated'


############################################################
###JACKKNIFE CODE
##FIRST SEPARATE THE X,Y,X DATA AND RANDOMS
xdat=np.array([obs.datx[bdx]])
ydat=np.array([obs.daty[bdx]])
zdat=np.array([obs.datz[bdx]])
xrand=np.array([obs.randx])
yrand=np.array([obs.randy])
zrand=np.array([obs.randz])

xd = xdat[0]
yd = ydat[0]
zd = zdat[0]
xr = xrand[0]
yr = yrand[0]
zr = zrand[0]

start2=time()

#SPLIT INTO THE 'SPLITVAL' AMOUNT OF SUBSETS (note, these return an unequal amount of data)
xsplitval=5
ysplitval=3
#Now I split the data at the max and min x and y data to include all of the points in the jackknife routine
xspvals= np.linspace(min(np.hstack((xr,xd))),max(np.hstack((xr,xd))),xsplitval)
yspvals= np.linspace(min(np.hstack((yr,yd))),max(np.hstack((yr,yd))),ysplitval)


#I am iterating over the 
XIjk=[]
RRjk=[]
for i in range(xsplitval-1):
	for j in range(ysplitval-1):
		#Make the proper cuts to pull out a small subset of the data
		cut1 = (xd>=xspvals[i])&(xd<=xspvals[i+1])&(yd>=yspvals[j])&(yd<=yspvals[j+1])
		cut2 = (xr>=xspvals[i])&(xr<=xspvals[i+1])&(yr>=yspvals[j])&(yr<=yspvals[j+1])
		#Invert these cuts so I am including all but the above cuts in my analysis
		dcut=np.invert(cut1)
		rcut=np.invert(cut2)
		#Separate the data from the part I leave out and combine to get [[x1,y1,z1,],[x2,y2,z2]...]
		xdatsplit=xd[dcut]
		ydatsplit=yd[dcut]
		zdatsplit=zd[dcut]
		xrandsplit=xr[rcut]
		yrandsplit=yr[rcut]
		zrandsplit=zr[rcut]
		oX=np.concatenate((np.array([xdatsplit]).T,np.array([ydatsplit]).T,np.array([zdatsplit]).T),axis=1)
		rX=np.concatenate((np.array([xrandsplit]).T,np.array([yrandsplit]).T,np.array([zrandsplit]).T),axis=1)
		#Runs the correlation function here and output the Xi and RR data
		TPCF=two_point(oX,radii,method='landy-szalay',data_R=rX, random_state=None)
		xi=np.asarray(TPCF[0])
		dd=np.asarray(TPCF[1])
		dr=np.asarray(TPCF[2])
		rr=np.asarray(TPCF[3])
		RRjk.append(rr)
		XIjk.append(xi)

##Compute the variances (Here I only do the main diagonal, I should update to include the full covariance matrix)
c=0
C=[]
for i in range(len(XIjk)):
	print c
	sig = (np.sqrt(RRjk[i]/(1.0*RRfull))*(XIjk[i]-xifull))**2
	print sig
	c += sig
	if i+1 == len(XIjk):
		C.append(c)
	
stdev=C[0]**0.5

end2=time()
print 'Done Jackknife'
print 'C=', C
print 'stdev=', stdev
print end2-start2,'s'
#Write the initial run to a file (this is the non-Jackknife run)
file=open('SpIES_fullcorr_test_v2.txt','w')
file.write('s DD DR RR Xifull sigmajk\n')
for i in range(len(sep)):
	s= sep[i]
	dd= DDfull[i]
	dr= DRfull[i]
	rr=RRfull[i]
	xif=xifull[i]
	sigma=stdev[i]
	file.write(str(s)+' '+str(dd)+' '+str(dr)+' '+str(rr)+' '+str(xif)+' '+str(sigma)+'\n')
file.close()
#Writhe the Jackknife data to a table
f=open('SpIES_jackknife_results_v2.txt','w')
f.write('s RRjk Xijk \n')
for i in range((xsplitval-1)*(ysplitval-1)):
	for j in range(len(RRjk[i])):
		rrjk=RRjk[i][j]
		xijk=XIjk[i][j]
		se=sep[j]
		f.write(str(se)+' '+str(rrjk)+' '+str(xijk)+'\n')
f.close()


'''
#Test plot with a comparison to the ROSS 2009 data

#OPEN THE ROSS 2009 DATA TO COMPARE
data=open('k_output_UNI22.dat','r')

r=[]
Xi=[]

for i in data.readlines():
	val=i.split()
	r.append(float(val[0]))
	Xi.append(float(val[2]))
	
R=np.array(r)

print 'Plotting'


figure(1)
plot(10**R,Xi, linewidth=2, label='Ross2009')
scatter(sep,xifull,color='g',label='SpIES JK')
errorbar(sep,xifull,yerr=stdev,linestyle="None",linewidth=2,color='g')
xscale('log')
yscale('log')
xlabel(r's(h$^{-1}$Mpc)')
ylabel(r'$\xi$(s)')
legend(scatterpoints=1)
savefig('SpIES_Jackknife_test_v2.png')
show()

'''

