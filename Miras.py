import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.constants import c
from scipy import interpolate
import math
from scipy import integrate
from copy import copy

# constant
Vmax = 200
Vmin = -200
Vsteps = 6
indice_max = 12
whugg=[3600,4000]
wo2b = [6866,6950]
wo2a= [7590,7724]
wH2oa=[8950,9860]
wH2ob=[8100,8300]
c = c/1000

def integrale(V,raie):
	Ntemp = np.abs(raie)
	isole = Ntemp[25:34]
	I1 = integrate.simpson(isole,dx=Vsteps)
	return(I1)



# #get data from ascii
data = np.genfromtxt('Miras-stars-NeoNarval-V/UHer_2024-04-12-stokesI.ascii.sum')

# #observed data
wavelength = data[:,0]
intensity = data[:,1] - 1
sigma=data[:,2]

#Velocity grid
V=np.arange(Vmin,Vmax,Vsteps)

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
pickle = pickle.load(open('t3000g+3.0.sel_atmo+mol.pickle', 'rb'))

mask = pickle['wave_vac']
poids = pickle['depth']
weightI=pickle['weight_I']
lande=pickle['lande']
Ismolec = pickle['is_molecule?']
Weight_Z=pickle['weight_Zeeman']



#Array that will contain all data from the rays
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

Wla = SommeB/sommeW
# ################### Q





#get data from ascii
data2 = np.genfromtxt('Miras-stars-NeoNarval-V/UHer_2024-04-12-stokesV.ascii.sum')
dataN = np.genfromtxt('Miras-stars-NeoNarval-V/UHer_2024-04-12-stokesN.ascii.sum')


#observed data
wavelength2= data2[:,0]
intensity2 = data2[:,1]
sigma2=data2[:,2]

wavelengthN= dataN[:,0]
intensityN = dataN[:,1]
sigmaN=dataN[:,2]

#Array that will contain all data from the rays
NN=np.zeros((len(V),len(mask)))
NsN=np.zeros((len(V),len(mask)))

NQ=np.zeros((len(V),len(mask)))
NsQ=np.zeros((len(V),len(mask)))

#Array that will contain the filtered data (V/N and sigma)
iN=[]
i2=[]

s2=[]
sN=[]



#filtration of data
for k in (range (len(wavelength2))):
	
	if wavelength[k]<300:
		i2.append(np.nan) 
		s2.append(np.nan)
		iN.append(np.nan) 
		sN.append(np.nan)
	elif wavelength[k]>12000:
		i2.append(np.nan) 
		s2.append(np.nan)
		iN.append(np.nan) 
		sN.append(np.nan)

	else :
		i2.append(intensity2[k]) 
		s2.append(sigma2[k])
		iN.append(intensityN[k]) 
		sN.append(sigmaN[k])
		
i2=np.array(i2)
s2=np.array(s2)
iN=np.array(iN)
sN=np.array(sN)




for k in (range (len(mask))) :
	
	#transform the volicity grid in wavelength
	wmin= mask[k]*Vmin/c + mask[k]
	wmax= mask[k]*Vmax/c + mask[k]
	iii=np.argwhere(np.logical_and(wavelength > wmin-0.2, wavelength < wmax+0.2))
	w=wavelength2[iii]
	intens = i2[iii]
	sigm=s2[iii]	
	
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
		
		NQ[:,k]=y

		NN[:,k]=yN
	
		NsQ[:,k]=1/z #Take the opposite by defintion of the matrix
		NsN[:,k]=1/zN
	
	
	except : 
		 
		#column of nan if interpolation doesn't work
		a = np.empty(len(V)) #Array of nan
		a[:]=np.nan
		NQ[:,k]=a
		NsQ[:,k]=a
		NN[:,k]=a
		NsN[:,k]=a
		
	
	#Test if we know the lande of the ray (99 is unknown)
	if lande[k]==99 :

		
		#Take the data from the last ray 
		NQpositif = copy(NQ[:,k]) #positive lande 
		NQnegatif = -copy(NQ[:,k]) #negative lande
		
		NNpositif = copy(NN[:,k])
		NNnegatif = -copy(NN[:,k])

		#Creation of a matrix with k column that contains all data at this index of the loop
		NQtemp = NQ[:,0:k+1]
		NQtemp[:,k]=NQpositif #put the positive value as the last column
		NsQtemp = NsQ[:,0:k+1] #The sigma doesn't change
		
		NNtemp = NN[:,0:k+1]
		NNtemp[:,k]=NNpositif #put the positive value as the last column
		NsNtemp = NsN[:,0:k+1]
		
		#Equation of WLA 
		mult=NsQtemp*NQtemp
		Somme=np.nansum(mult,axis=1)
		sommek=np.nansum(NsQtemp,axis=1)
		Wla2positif = Somme/sommek
		Ipositif = integrale(V, Wla2positif) #integral between the max/min of the pseudo profile created until now
		
		#Same things for Stokes N
		multN=NsNtemp*NNtemp
		SommeN=np.nansum(multN,axis=1)
		sommeNk=np.nansum(NsNtemp,axis=1)
		WlaNpositif = SommeN/sommeNk
		INpositif = integrale(V, WlaNpositif) 
		
		
		NQnegatiftemp = copy(NQ[:,0:k+1])
		NQnegatiftemp[:,k]=NQnegatif #put the negative value as the last column
		NsQtemp = NsQ[:,0:k+1]
		
		NNnegatiftemp = copy(NN[:,0:k+1])
		NNnegatiftemp[:,k]=NNnegatif #put the positive value as the last column
		NsNtemp = NsN[:,0:k+1]
		
		#Equation of WLA 
		mult=NsQtemp*NQnegatiftemp
		Somme=np.nansum(mult,axis=1)
		sommek=np.nansum(NsQtemp,axis=1)
		Wla2negatif = Somme/sommek
		Inegatif = integrale(V, Wla2negatif) #integral between the max/min of the pseudo profile created until now
		
		#Same things for Stokes N
		multN=NsNtemp*NNnegatiftemp
		SommeN=np.nansum(multN,axis=1)
		sommeNk=np.nansum(NsNtemp,axis=1)
		WlaNnegatif = SommeN/sommeNk
		INnegatif = integrale(V, WlaNnegatif)
	
		
		#Test to take the maximum value of the integral between lande +1 or -1
		if Ipositif < Inegatif : 
			NQ[:,k]=NQnegatif
			lande[k] = -1
		else : 
			NQ[:,k]=NQpositif
			lande[k]=1
			
		if INpositif < INnegatif : 
			NN[:,k]=NNnegatif
		else : 
			NN[:,k]=NNpositif
	else : 
		
		NQ[:,k] = NQ[:,k]
		NN[:,k] = NN[:,k]




#Equation of WLA for Stokes V after all the +1/-1
mult=NsQ*NQ
Somme=np.nansum(mult,axis=1)
sommek=np.nansum(NsQ,axis=1)

Wla_with_lande_plus_or_minus = Somme/(sommek)



#Equation of WLA for Stokes N after all the +1/-1
multN=NsN*NN
SommeN=np.nansum(multN,axis=1)
sommekN=np.nansum(NsN,axis=1)

Wla_with_N = SommeN/(sommekN)




##############ALL PLOTS 


# fig=plt.figure()

# ax1 = plt.subplot(1,1,1)
# # ax1.plot(V,Sommeintens+1.1,label='Stokes I SLA',color='forestgreen')
# ax1.plot(V,Wla+1.05,label='Stokes I WLA', color='darkred')
# plt.xlabel('v(km/s)')
# plt.ylabel('I/Ic') 	
# plt.legend()
# plt.grid(True)




ax2 = plt.subplot(1,1,1)
ax2.plot(V,Wla_with_N-0.0012,'b',label='Stokes N WLA. Date : 2014-04-11')
ax2.plot(V,Wla_with_lande_plus_or_minus+0.00025,'r',label='Stokes V WLA with Lande +1 or -1. Date : 2014-04-11')
plt.xlabel('v(km/s)')
plt.ylabel('V/Ic') 	
plt.legend()
plt.title('Comparaison Stokes N et V')
plt.grid(True)
plt.show()

 	
	
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	