import numpy as np
import scipy as sp
import pyfits as pf
from matplotlib import *
from pylab import *

#Import the data 

dr10data=pf.open("SpIES_ch1ch2_dr10Q_match.fits")[1].data
Spdata=pf.open("SpIESch1ch2.fits")[1].data
#alldata=pf.open("Sp_all_ch1ch2_match.fits")[1].data
LRGdata=pf.open("SpIES_DR9gal_match.fits")[1].data


#Define the Aperture corrections
apcor1=0.765 #Aperture correction for channel 1 for the 1.9" aperture
apcor2=0.740 #Aperture correction for channel 2 for the 1.9" aperture
fluxcorr1=1.0*10**12 * 1./(180/np.pi)**2 * 1.0/(3600.)**2 * (0.6)**2/apcor1 #Converting MJy/sr to uJy including the aperture correction in channel 1
fluxcorr2=1.0*10**12 * 1./(180/np.pi)**2 * 1.0/(3600.)**2 * (0.6)**2/apcor2 #Converting MJy/sr to uJy including the aperture correction in channel 2

print 'Uploaded'

#Create AB magnitudes and colors from the SpIES data using M=-2.5log10(flux)+23.9
SpIEScormag1=-2.5*np.log10(Spdata.FLUX_APER_1[:,1]*fluxcorr1)+23.9
SpIEScormag2=-2.5*np.log10(Spdata.FLUX_APER_2[:,1]*fluxcorr2)+23.9
SpIEScolor12=SpIEScormag1-SpIEScormag2

#Create AB magnitudes and colors from the matched DR10 data using M=-2.5log10(flux)+23.9
dr10cormag1=-2.5*np.log10(dr10data.FLUX_APER_1[:40000,1]*fluxcorr1)+23.9
dr10cormag2=-2.5*np.log10(dr10data.FLUX_APER_2[:40000,1]*fluxcorr2)+23.9
dr10color12=dr10cormag1-dr10cormag2
'''
#Create AB magnitudes and colors from the matched Sp all data using M=-2.5log10(flux)+23.9
allcormag1=-2.5*np.log10(alldata.FLUX_APER_1[:,1]*fluxcorr1)+23.9
allcormag2=-2.5*np.log10(alldata.FLUX_APER_2[:,1]*fluxcorr2)+23.9
allcolor12=allcormag1-allcormag2
'''
#Create AB magnitudes and colors from the matched LRG data using M=-2.5log10(flux)+23.9
LRGcormag1=-2.5*np.log10(LRGdata.FLUX_APER_1[:,1]*fluxcorr1)+23.9
LRGcormag2=-2.5*np.log10(LRGdata.FLUX_APER_2[:,1]*fluxcorr2)+23.9
LRGcolor12=LRGcormag1-LRGcormag2

print 'mags converted'

#make the color cuts here and divide by the area of the survey
area=1.63*44
colorcut=linspace(-1,1,50)
numQ=[]
numL=[]
numO=[]
for i in range(len(colorcut)):
	lrg=LRGcolor12>colorcut[i] #Creates indicies where the conditional is true
	qso=dr10color12>=colorcut[i]#Creates indicies where the conditional is true
	spi=SpIEScolor12>=colorcut[i]#Creates indicies where the conditional is true
	numQ.append(len(Spdata.FLUX_APER_2[:,1][qso])/area) #Store the number per area as a function of color
	numL.append(len(LRGdata.FLUX_APER_2[:,1][lrg])/area)#Store the number per area as a function of color
	numO.append(len(Spdata.FLUX_APER_2[:,1][spi])/area)#Store the number per area as a function of color

	
print 'colors cut'
'''
figure(1)
plot(colorcut,numQ)
xlabel('Ch1-Ch2 [AB]')
ylabel('#Quasars/deg^2')
savefig('Num_den_qso.png')

figure(2)
plot(colorcut,numO)
xlabel('Ch1-Ch2 [AB]')
ylabel('#Objects/deg^2')
savefig('Num_den_total.png')

figure(3)
plot(colorcut,numL)
xlabel('Ch1-Ch2 [AB]')
ylabel('#LRG/deg^2')
savefig('Num_den_lrg.png')

figure(4)
scatter(SpIEScormag1,SpIEScolor12,s=1,alpha=0.3,label='All')
scatter(dr10cormag1,dr10color12,s=1,color='r',label='QSO')
scatter(LRGcormag1,LRGcolor12,s=1,color='g',label='LRG')
#scatter(allcormag1,allcolor12,s=1,color='b')
xlabel('Ch1 [AB]')
ylabel('Ch1-Ch2 [AB]')
legend()
savefig('Color_Ch1mag.png')
'''
figure(7)
scatter(LRGcormag1,LRGcolor12,s=1,color='g',label='LRG')
scatter(dr10cormag1,dr10color12,s=1,color='r',label='QSO')
xlabel('Ch1 [AB]')
ylabel('Ch1-Ch2 [AB]')

figure(5)
scatter(LRGcormag1,LRGcolor12,s=1,color='g',label='LRG')
xlabel('Ch1 [AB]')
ylabel('Ch1-Ch2 [AB]')

figure(6)
scatter(dr10cormag1,dr10color12,s=1,color='r',label='QSO')
xlabel('Ch1 [AB]')
ylabel('Ch1-Ch2 [AB]')
show()


