import numpy as np
from astropy.io import fits as pf
from matplotlib import *
#use('Agg')
from pylab import *

cut=pf.open('Cartesian_coords.fits')[1].data
bdx=(cut.datx!=0)
#Combine the x,y,z data from the imported file in the form [[x1,y1,z1],[x2,y2,z2]...]
cutX=np.concatenate((np.array([cut.datx[bdx]]).T,np.array([cut.daty[bdx]]).T,np.array([cut.datz[bdx]]).T),axis=1)
randX=np.concatenate((np.array([cut.randx]).T,np.array([cut.randy]).T,np.array([cut.randz]).T),axis=1)
print len(cutX),len(randX), len(cutX)*len(randX)
print len(cutX)*(len(cutX)-1)/2
#Boss data
BOSS=open('Anderson_2013_CMASSDR10_corrfunction_x0x2_postrecon.dat','r')
Head=BOSS.readline()
data=BOSS.readlines()
S=[]
xiboss=[]
r2xiboss=[]
for i in data:
	S.append(float(i.split()[0]))
	xiboss.append(float(i.split()[1]))
	r2xiboss.append(float(i.split()[0])**2*float(i.split()[1]))
#Import 2slaq check
#Other file 'k_output_2SLAQ_Edin_cor_FtF.dat'
Slaq=open('k_output_jack_perl_full_newcor_v2_ed.dat','r')
slaqdat=Slaq.readlines()
R=[]
xislaq=[]
r2xislaq=[]
DR=[]
DD=[]
RR=[]
for i in range(len(slaqdat)-1):
	#print slaqdat[i].split()[1], slaqdat[i].split()[6]
	R.append(float(slaqdat[i].split()[1]))
	xislaq.append(float(slaqdat[i].split()[6]))
	r2xislaq.append(float(slaqdat[i].split()[1])**2*float(slaqdat[i].split()[6]))
	DD.append(float(slaqdat[i].split()[4]))
	DR.append(float(slaqdat[i].split()[5]))
	RR.append(float(slaqdat[i].split()[8]))


#Import the File I created
Test=open('Pair_counts_2slaq_small.txt','r')
head=Test.readline()
data=Test.readlines()
sep=[]
dd=[]
rr=[]
dr=[]

for i in range(len(data)):
	sep.append(float(data[i].split()[0]))	
	dd.append(float(data[i].split()[1]))
	rr.append(float(data[i].split()[2]))
	dr.append(float(data[i].split()[3]))
#Create the correlation function estimators
xi=[]
r2xi=[]
phxi=[]
r2phxi=[]
hxi=[]
r2hxi=[]
s=[]
for i in range(len(dd)):
	if dr[i]!=0:
		s.append(sep[i])
		n_d=len(cutX)
		n_r=len(randX)
		datdat=dd[i]/(n_d*(n_d-1)/2.)
		randrand=rr[i]/(n_r*(n_r-1)/2.)
		datrand=dr[i]/(n_d*n_r)
		LS=(randrand+datdat-2*datrand)/randrand #Landy-Szalay
		PH=datdat/randrand - 1			#Peebles-Hauser
		H= (datdat*randrand)/datrand**2 - 1	#Hamilton
		xi.append(LS)
		r2xi.append(sep[i]**2*LS)
		phxi.append(PH)
		r2phxi.append(PH*sep[i]**2)
		hxi.append(H) 
		r2hxi.append(H*sep[i]**2)




figure(6)
loglog(s,xi,label='LS')
loglog(s,phxi,label='PH')
loglog(s,hxi,label='H')
scatter(R,xislaq,color='g',label='2slaq')
xlabel(r's(h$^{-1}$Mpc)')
ylabel(r'$\xi$(s)')
legend()

figure(7)
plot(s,r2xi,label='LS')
plot(s,r2phxi,label='PH')
plot(s,r2hxi,label='H')
scatter(R,r2xislaq,color='g',label='2slaq')
legend()
'''
figure(3)
loglog(sep,dr)
loglog(R,DR,color='g')

figure(4)
loglog(sep,dd)
loglog(R,DD,color='g')

figure(5)
loglog(sep,rr)
loglog(R,RR,color='g')
'''
show()




