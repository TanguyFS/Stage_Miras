import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.constants import c
from scipy import interpolate
import math
from scipy import integrate
from copy import copy
from scipy.ndimage import uniform_filter1d
from scipy.optimize import least_squares

# constant
Vmax = 200
Vmin = -200
Vsteps = 6 #Spectral resolution
indice_max = 47 
indice_min=32
#Telluric rays
whugg=[3600,4000]
wo2b = [6866,6950]
wo2a= [7590,7724]
wH2oa=[8950,9860]
wH2ob=[8100,8300]

c = c/1000 # Speed of light in km/s
facteurs_mult = 1


#Function that I will minimize 
def integrale(param_libre0,V,NV,NsV,param,Indice,NN,NsN):
	
	#normalisation by the sum of the unknown lande factor
	normal = np.sum(param_libre0[:]) 

	param=np.array(param)
	param_libre=np.array(param_libre0)
	
	#put the free parameters at the right place in the array of lande factor
	for k in range (Inconnue) : 
		param[int(Indice[k])] = param_libre[k]
	
	#Equation of WLA
	Nv2=NV*param 
	mult=NsV*Nv2
	Somme=np.nansum(mult,axis=1)
	sommek=np.nansum(NsV,axis=1)
	raie = Somme/(sommek)
	
	#put the mean at 0
	offset=np.mean(raie)	
	raie_offset = raie-offset

	#Frequency phase to smooth the signal
	fft=np.fft.fft(raie_offset)
	fft=np.fft.fftshift(fft)
	fft[0:25] = (-0-0j)
	fft[40:] = (-0-0j)
	fft=np.fft.ifftshift(fft)
	ifft=np.fft.ifft(fft)

	#Integration between the index corresponding of max Stokes I
	isole = ifft[indice_min:indice_max]
	isole=np.abs(isole) #Absolute value because Stokes V is symetrical

	IV = np.nansum(isole/normal)
	
	
	
	#Repeat all the operations for Stokes N
	NN2=NN*param 
	

	multN=NsN*NN2
	SommeN=np.nansum(multN,axis=1)
	sommekN=np.nansum(NsN,axis=1)
	raieN = SommeN/(sommekN)
	offset=np.mean(raieN)	
	raieN_offset = raieN-offset


	fftN=np.fft.fft(raieN_offset)
	fftN=np.fft.fftshift(fftN)	
	fftN[0:30] = (-0-0j)
	fftN[35:] = (-0-0j)	
	fftN=np.fft.ifftshift(fftN)
	ifftN=np.fft.ifft(fftN)

	isoleN = ifftN[indice_min:indice_max]
	isoleN=np.abs(isoleN)
	IN = np.nansum(isoleN/normal)
	


#Return the integral of Stokes N and the inverse of 
#the integral of Stokes V because I want to minimize Stokes N 
#and maximize Stokes V
	return(IN+1/IV) 


#Velocity grid
V=np.arange(Vmin,Vmax,Vsteps)

#get data from ascii
data = np.genfromtxt('UMon-Feb2021-NeoNarval-V/UMon_2021-02-24-stokesV.ascii')

#observed data
wavelength = data[:,0]
intensity = data[:,1] - 1
sigma=data[:,2]


#Array that will contain the filtered data (I and sigma)
i=[]
s=[]


#filtration of data
for k in (range (len(wavelength))):
 	
 	if wavelength[k]> wo2a[0] and wavelength[k] < wo2a[1]:
 		 i.append(np.nan) 
 		 s.append(np.nan)
 	elif wavelength[k]> wo2b[0] and wavelength[k] < wo2b[1]:
 		 i.append(np.nan) 
 		 s.append(np.nan)

 	elif wavelength[k]> wH2oa[0] and wavelength[k] < wH2oa[1]:
 		 i.append(np.nan)
 		 s.append(np.nan)

 	elif wavelength[k]> wH2ob[0] and wavelength[k] < wH2ob[1]:
 		 i.append(np.nan)
 		 s.append(np.nan)

 	elif wavelength[k]> whugg[0] and wavelength[k] < whugg[1]:
 		 i.append(np.nan)
 		 s.append(np.nan)

 	else :
 		 i.append(intensity[k]) 
 		 s.append(sigma[k])
i=np.array(i)
s=np.array(s)


#get mask
pickle = pickle.load(open('t5000g+3.0.sel_atmo+mol.pickle', 'rb'))

mask = pickle['wave_vac'] #wavelength of all the rays
poids = pickle['depth'] 
weightI=pickle['weight_I']
lande=pickle['lande'] #lande factor
Ismolec = pickle['is_molecule?'] #if the ray is molecular or not 
Weight_Z=pickle['weight_Zeeman']



# #Array that will contain all data from the rays
Nintens=np.zeros((len(V),len(mask)))
Ns=np.zeros((len(V),len(mask)))

for k in (range (len(mask))) :
 	
 	#transform the volicity grid in wavelength
 	wmin= mask[k]*Vmin/c + mask[k]
 	wmax= mask[k]*Vmax/c + mask[k]
 	iii=np.argwhere(np.logical_and(wavelength > wmin-0.2, wavelength < wmax+0.2))
 	w=wavelength[iii]
 	intens = i[iii]
 	sigm=s[iii]	
 	Vk = c*(w-mask[k])/mask[k]
 	try : 
 		 #interpolation
 		 f=interpolate.interp1d(np.squeeze(Vk), np.squeeze(intens))	 
 		 g=interpolate.interp1d(np.squeeze(Vk), np.squeeze(sigm))	
		
 		 y = f(V)
 		 z=g(V)
 		 Nintens[:,k]=y
		
 		 Ns[:,k]=1/z
 	except : 
 		 #column of nan if interpolation doesn't work
 		 a = np.empty(len(V))
 		 a[:]=np.nan
 		 Nintens[:,k]=a
 		 Ns[:,k]=a

#Equation of WLA
multintens=Ns*Nintens
SommeB=np.nansum(multintens,axis=1)
sommeW=np.nansum(Ns,axis=1)

WlaI = SommeB/sommeW


#################### STOKES V





#get data from ascii
dataV = np.genfromtxt('UMon-Feb2021-NeoNarval-V/UMon_2021-02-14-stokesV.ascii')
dataN = np.genfromtxt('UMon-Feb2021-NeoNarval-V/UMon_2021-02-14-stokesN.ascii')


#observed data
wavelength2= dataV[:,0]
intensity2 = dataV[:,1]
sigma2=dataV[:,2]

wavelengthN= dataN[:,0]
intensityN = dataN[:,1]
sigmaN=dataN[:,2]

#Array that will contain all data from the rays
NN=np.zeros((len(V),len(mask)))
NsN=np.zeros((len(V),len(mask)))

NV=np.zeros((len(V),len(mask)))
NsV=np.zeros((len(V),len(mask)))

lande0=np.zeros((len(V),len(mask)))


#Array that will contain the filtered data (V/N and sigma)
iN=[]
iV=[]

sV=[]
sN=[]



#filtration of data
for k in (range (len(wavelength2))):
	
	if wavelengthN[k]<300:
		iV.append(np.nan) 
		sV.append(np.nan)
		iN.append(np.nan) 
		sN.append(np.nan)
	elif wavelengthN[k]>12000:
		iV.append(np.nan) 
		sV.append(np.nan)
		iN.append(np.nan) 
		sN.append(np.nan)

	else :
		iV.append(intensity2[k]) 
		sV.append(sigma2[k])
		iN.append(intensityN[k]) 
		sN.append(sigmaN[k])
		
iV=np.array(iV)
sV=np.array(sV)
iN=np.array(iN)
sN=np.array(sN)


#construct the vector of known and unknown lande factor
Inconnue = len(np.argwhere(lande==99)) #value of all the unknown lande factor 
param = [] #vector of all lande factor
param_libre0 = np.zeros(Inconnue)  #vector of unknown lande factor
Indice = np.zeros(Inconnue) #index at which I found the unknown lande factor in the mask
j=0
for k in range (len(mask)) :
	if lande[k] ==99 :  #if the lande factor is unknown
		param_libre0[j]=1 #initialisation 
		Indice[j]=k #index
		param.append(0) #placeholder 
		j+=1
		
	else : 
		param.append(lande[k])
		


for k in (range (len(mask))) :
	
	#transform the volicity grid in wavelength
	wmin= mask[k]*Vmin/c + mask[k]
	wmax= mask[k]*Vmax/c + mask[k]
	iii=np.argwhere(np.logical_and(wavelengthN > wmin-0.2, wavelengthN < wmax+0.2))
	w=wavelength2[iii]
	intens = iV[iii]
	sigm=sV[iii]	
	
	wN=wavelengthN[iii]
	intensN= iN[iii]
	sigmN=sN[iii]
	
	Vk = c*(w-mask[k])/mask[k]
	try : 
		#interpolation
		f=interpolate.interp1d(np.squeeze(Vk), np.squeeze(intens))	 
		g=interpolate.interp1d(np.squeeze(Vk), np.squeeze(sigm))	
		
		fN=interpolate.interp1d(np.squeeze(Vk), np.squeeze(intensN))	 
		gN=interpolate.interp1d(np.squeeze(Vk), np.squeeze(sigmN))	
		
		y = f(V)
		z=g(V)
		yN = fN(V)
		zN=gN(V)
		
		NV[:,k]=y

		NN[:,k]=yN
	
		NsV[:,k]=1/z #Take the opposite by defintion of the matrix
		NsN[:,k]=1/zN
	
	
	except : 
		 
		#column of nan if interpolation doesn't work
		a = np.empty(len(V)) #Array of nan
		a[:]=np.nan
		NV[:,k]=a
		NsV[:,k]=a
		NN[:,k]=a
		NsN[:,k]=a
	
		
	
bounds = np.ones(Inconnue) #lande factor are limited from -1 to +1

#Using least square to maximize the signal between the index choosen by adjusting the uknown lande factor
res = least_squares(integrale, param_libre0, args=(V,NV,NsV,param,Indice,NN,NsN),bounds = (-bounds,bounds),verbose=2)

param_LS = res.x #results of the unknown Lande factor
normal = 1 #normalization

for k in range (1,Inconnue) : 
	param[int(Indice[k])] = param_LS[k] #Put the new Lande factor in the right place
	
lande_fin = param
param=np.array(param)


NV_fin = NV*lande_fin #multiplication of the rays that had an uknown lande factor, by the value

#Equation of WLA for Stokes V
mult_LS=NsV*NV_fin
Somme_LS=np.nansum(mult_LS,axis=1)
sommek_LS=np.nansum(NsV,axis=1)
WLA_LS = Somme_LS/(sommek_LS*normal)


#Repeat operations for Stokes N
NN_fin = NN*lande_fin

multN_LS=NsN*NN_fin
SommeN_LS=np.nansum(multN_LS,axis=1)
sommekN_LS=np.nansum(NsN,axis=1)
WLAN_LS = SommeN_LS/(sommekN_LS*normal)



		
##############ALL PLOTS 


#plot for Stokes I
fig=plt.figure()

ax1 = plt.subplot(1,1,1)
# ax1.plot(V,Sommeintens+1.1,label='Stokes I SLA',color='forestgreen')
ax1.plot(V,WlaI+1.05,label='Stokes I WLA', color='darkred')
plt.xlabel('v(km/s)')
plt.ylabel('I/Ic') 	
plt.legend()
plt.grid(True)



#Plot to compare Stokes V and N
ax2 = plt.subplot(2,1,1)
ax2.plot(V,WLAN_LS,'b',label='Stokes N WLA with Least Square. Umon TBL 2021-02-27')
ax2.plot(V,WLA_LS,'r',label='Stokes V WLA with Least Square. Umon TBL 2021-02-27')
plt.xlabel('v(km/s)')
plt.ylabel('V/Ic') 	
plt.vlines(V[indice_min], np.min(WLA_LS),np.max(WLA_LS),linestyles='--',color='green')
plt.vlines(V[indice_max], np.min(WLA_LS),np.max(WLA_LS),linestyles='--',color='green')
plt.legend()
plt.title('Comparaison Stokes N et V')
plt.grid(True)
plt.show()


##plot of histogram for the lande factor that I found
# ax2 = plt.subplot(2,1,2)
# ax2.hist(param_LS[1:], bins=100)
	
#plt.hist(param_LS, bins= 500)


	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	