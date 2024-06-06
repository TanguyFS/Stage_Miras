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
# intensity = data[:,1] - 1
# sigma=data[:,2]

# i=[]

# s=[]

# for k in (range (len(wavelength))):
# 	
# 	if wavelength[k]> wo2a[0] and wavelength[k] < wo2a[1]:
# 		i.append(np.nan) 
# 		s.append(np.nan)
# 	elif wavelength[k]> wo2b[0] and wavelength[k] < wo2b[1]:
# 		i.append(np.nan) 
# 		s.append(np.nan)

# 	elif wavelength[k]> wH2oa[0] and wavelength[k] < wH2oa[1]:
# 		i.append(np.nan)
# 		s.append(np.nan)

# 	elif wavelength[k]> wH2ob[0] and wavelength[k] < wH2ob[1]:
# 		i.append(np.nan)
# 		s.append(np.nan)

# 	elif wavelength[k]> whugg[0] and wavelength[k] < whugg[1]:
# 		i.append(np.nan)
# 		s.append(np.nan)

# 	else :
# 		i.append(intensity[k]) 
# 		s.append(sigma[k])



V=np.arange(Vmin,Vmax,Vsteps)
# i=np.array(i)
# s=np.array(s)


#get mask
pickle = pickle.load(open('t3000g+3.0.sel_atmo+mol.pickle', 'rb'))

mask = pickle['wave_vac']
poids = pickle['depth']
weightI=pickle['weight_I']
lande=pickle['lande']
Ismolec = pickle['is_molecule?']
Weight_Z=pickle['weight_Zeeman']

# Nintens=np.zeros((len(V),len(mask)))
# Nintens2=np.zeros((len(V),len(mask)))
# Ns=np.zeros((len(V),len(mask)))
# Ns1=np.zeros((len(V),len(mask)))

# for k in (range (len(mask))) :
#  	wmin= mask[k]*Vmin/c + mask[k]
#  	wmax= mask[k]*Vmax/c + mask[k]
#  	iii=np.argwhere(np.logical_and(wavelength > wmin-0.2, wavelength < wmax+0.2))
#  	w=wavelength[iii]
#  	intens = i[iii]
#  	sigm=s[iii]	
#  	Vk = c*(w-mask[k])/mask[k]
#  	try : 
# 		f=interpolate.interp1d(np.squeeze(Vk), np.squeeze(intens))	 
# 		g=interpolate.interp1d(np.squeeze(Vk), np.squeeze(sigm))	
# 		
# 		y = f(V)
# 		z=g(V)
# 		Nintens[:,k]=y
# 		
# 		Ns[:,k]=1/z
#  	except : 
# 		a = np.empty(len(V))
# 		a[:]=np.nan
# 		Nintens[:,k]=a
# 		Ns[:,k]=a

# Sommeintens = np.nanmean(Nintens,axis=1)

# multkkk=Ns*Nintens
# Sommekk=np.nansum(multkkk,axis=1)
# sommekkk=np.nansum(Ns,axis=1)

# Wla = Sommekk/sommekkk
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

NN=np.zeros((len(V),len(mask)))
NsN=np.zeros((len(V),len(mask)))

NQ=np.zeros((len(V),len(mask)))
NQ2=np.zeros((len(V),len(mask)))
NQ3=np.zeros((len(V),len(mask)))
Ns2=np.zeros((len(V),len(mask)))
Ns3=np.zeros((len(V),len(mask)))


iN=[]
i2=[]

s2=[]
sN=[]

for k in (range (len(wavelength2))):
	
	if wavelength[k]<6000:
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
		f=interpolate.interp1d(np.squeeze(Vk), np.squeeze(intens))	 
		g=interpolate.interp1d(np.squeeze(Vk), np.squeeze(sigm))	
		
		fN=interpolate.interp1d(np.squeeze(Vk), np.squeeze(intensN))	 
		gN=interpolate.interp1d(np.squeeze(Vk), np.squeeze(sigmN))	
		
		y = f(V)
		z=g(V)
		yN = fN(V)
		zN=gN(V)
		
		NQ[:,k]=y
		NQ2[:,k]=y
		NQ3[:,k]=y
		NN[:,k]=yN
	
		Ns2[:,k]=1/z
		Ns3[:,k]=1/z
		NsN[:,k]=1/zN
	
	
	except : 
		a = np.empty(len(V))
		a[:]=np.nan
		NQ[:,k]=a
		NQ2[:,k]=a
		Ns2[:,k]=a
		Ns3[:,k]=a
		NN[:,k]=a
		NsN[:,k]=a
		
	
		
	if lande[k]==99 :
# 		NQ2[:,k]=a
		
		NQk = NQ[:,k]
		Ns2k = Ns2[:,k]
		NsNk=NsN[:,k]
		
		
		NQpositif = copy(NQ[:,k])
		NQnegatif = -copy(NQ[:,k])
		
		NNpositif = copy(NN[:,k])
		NNnegatif = -copy(NN[:,k])
		
		Ns2positif = Ns2k
		Ns2negatif = Ns2k
		
		NsNpositif = NsNk
		NsNnegatif = NsNk
		
	
		NQ3[:,k]=NQpositif
		
		NQtemp = NQ[:,0:k+1]
		NQtemp[:,k]=NQpositif
		Ns2temp = Ns2[:,0:k+1]
		Ns2temp[:,k]=Ns2k
		
		NNtemp = NN[:,0:k+1]
		NNtemp[:,k]=NNpositif
		NsNtemp = NsN[:,0:k+1]
		NsNtemp[:,k]=NsNk
		
		mult=Ns2temp*NQtemp
		Somme=np.nansum(mult,axis=1)
		sommek=np.nansum(Ns2temp,axis=1)
		Wla2positif = Somme/sommek
		Ipositif = integrale(V, Wla2positif)
		
		multN=NsNtemp*NNtemp
		SommeN=np.nansum(multN,axis=1)
		sommeNk=np.nansum(NsNtemp,axis=1)
		WlaNpositif = SommeN/sommeNk
		INpositif = integrale(V, WlaNpositif)
		
		
		NQnegatiftemp = copy(NQ[:,0:k+1])
		NQnegatiftemp[:,k]=NQnegatif
		Ns2temp = Ns2[:,0:k+1]
		Ns2temp[:,k]=Ns2k	
		
		NNnegatiftemp = copy(NN[:,0:k+1])
		NNnegatiftemp[:,k]=NNnegatif
		NsNtemp = NsN[:,0:k+1]
		NsNtemp[:,k]=NsNk	
		
		mult=Ns2temp*NQnegatiftemp
		Somme=np.nansum(mult,axis=1)
		sommek=np.nansum(Ns2temp,axis=1)
		Wla2negatif = Somme/sommek
		Inegatif = integrale(V, Wla2negatif)
		
		multN=NsNtemp*NNnegatiftemp
		SommeN=np.nansum(multN,axis=1)
		sommeNk=np.nansum(NsNtemp,axis=1)
		WlaNnegatif = SommeN/sommeNk
		INnegatif = integrale(V, WlaNnegatif)
	
		
		if Ipositif < Inegatif : 
			NQ[:,k]=NQnegatif
			Ns2[:,k]=Ns2negatif
			lande[k] = -1
		else : 
			NQ[:,k]=NQpositif
			Ns2[:,k]=Ns2positif
			lande[k]=1
			
		if INpositif < INnegatif : 
			NN[:,k]=NNnegatif
			NsN[:,k]=NsNnegatif
		else : 
			NN[:,k]=NNpositif
			NsN[:,k]=NsNpositif
	else : 
		
		NQ[:,k] = NQ[:,k]*1
		NN[:,k] = NN[:,k]*1

Sommelande = np.sum(lande)

Somme2 = np.nanmean(NQ2,axis=1)
mult3=Ns3*NQ2
Somme3=np.nansum(mult3,axis=1)
sommek3=np.nansum(Ns3,axis=1)

Wla_without_lande = Somme3/sommek3


mult=Ns3*NQ
Somme=np.nansum(mult,axis=1)
sommek=np.nansum(Ns3,axis=1)

Wla_with_lande_plus_or_minus = Somme/(sommek)

mult=Ns3*NQ3
Somme=np.nansum(mult,axis=1)
sommek=np.nansum(Ns3,axis=1)

Wla_with_lande_plus1 = Somme/(sommek)


multN=NsN*NN
SommeN=np.nansum(multN,axis=1)
sommekN=np.nansum(NsN,axis=1)

Wla_with_N = SommeN/(sommekN)

####### Barre d'erreurs





# We=Ns2
# Se=1/Ns2
# erreur = np.sqrt(np.nansum(We*Se*Se,axis=1)/(np.nansum(We,axis=1)*np.shape(NQ)[1]))







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
# ax2.plot(V,Wla_with_lande_plus1,'b--',label='Stokes V WLA with Lande +1')
plt.xlabel('v(km/s)')
plt.ylabel('V/Ic') 	
plt.legend()
plt.title('Comparaison Stokes N et V')
plt.grid(True)
plt.show()

# ax2=plt.subplot(1,1,1)
# ax2.plot(V,Wla2,label='Stokes V WLA', color='darkred')
# ax2.errorbar(V,Wla2,yerr=erreur,barsabove=True,ecolor="blue")

# plt.xlabel('v(km/s)')
# plt.ylabel('V/Ic') 
# plt.legend()
# plt.grid(True)
 	
	
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	