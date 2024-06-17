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
Vsteps = 3.7
indice_max = 73
indice_min=52
whugg=[3600,4000]
wo2b = [6866,6950]
wo2a= [7590,7724]
wH2oa=[8950,9860]
wH2ob=[8100,8300]
c = c/1000

def integrale(param_libre0,V,NV,NsV,param):
	
	param_libre = param_libre0@matrice_lande
	param=np.array(param)
	param_libre=np.array(param_libre)
# 	lande0=np.transpose(param+param_libre)
	
	lande0=np.nansum(np.dstack((param,param_libre)),2)	
	
	Nv2=NV*lande0
	

	mult=NsV*Nv2
	Somme=np.nansum(mult,axis=1)
	sommek=np.nansum(NsV,axis=1)
	raie = Somme/(sommek)
	
	Ntemp = np.abs(raie)
	Ntemp = uniform_filter1d(Ntemp,size=3)
	isole = Ntemp[indice_min:indice_max]
	I1 = np.nansum(isole)
	return(1/I1)


#Velocity grid
V=np.arange(Vmin,Vmax,Vsteps)

# # #get data from ascii
# data = np.genfromtxt('Miras-stars-NeoNarval-V/UHer_2024-04-12-stokesI.ascii.sum')

# # #observed data
# wavelength = data[:,0]
# intensity = data[:,1] - 1
# sigma=data[:,2]


# #Array that will contain the filtered data (I and sigma)
# i=[]
# s=[]


# #filtration of data
# for k in (range (len(wavelength))):
 	
#  	if wavelength[k]> wo2a[0] and wavelength[k] < wo2a[1]:
# 		 i.append(np.nan) 
# 		 s.append(np.nan)
#  	elif wavelength[k]> wo2b[0] and wavelength[k] < wo2b[1]:
# 		 i.append(np.nan) 
# 		 s.append(np.nan)

#  	elif wavelength[k]> wH2oa[0] and wavelength[k] < wH2oa[1]:
# 		 i.append(np.nan)
# 		 s.append(np.nan)

#  	elif wavelength[k]> wH2ob[0] and wavelength[k] < wH2ob[1]:
# 		 i.append(np.nan)
# 		 s.append(np.nan)

#  	elif wavelength[k]> whugg[0] and wavelength[k] < whugg[1]:
# 		 i.append(np.nan)
# 		 s.append(np.nan)

#  	else :
# 		 i.append(intensity[k]) 
# 		 s.append(sigma[k])
# i=np.array(i)
# s=np.array(s)


#get mask
pickle = pickle.load(open('t5000g+3.0.sel_atmo+mol.pickle', 'rb'))

mask = pickle['wave_vac']
poids = pickle['depth']
weightI=pickle['weight_I']
lande=pickle['lande']
Ismolec = pickle['is_molecule?']
Weight_Z=pickle['weight_Zeeman']



# #Array that will contain all data from the rays
# Nintens=np.zeros((len(V),len(mask)))
# Ns=np.zeros((len(V),len(mask)))

# for k in (range (len(mask))) :
# 	
# 	#transform the volicity grid in wavelength
#  	wmin= mask[k]*Vmin/c + mask[k]
#  	wmax= mask[k]*Vmax/c + mask[k]
#  	iii=np.argwhere(np.logical_and(wavelength > wmin-0.2, wavelength < wmax+0.2))
#  	w=wavelength[iii]
#  	intens = i[iii]
#  	sigm=s[iii]	
#  	Vk = c*(w-mask[k])/mask[k]
#  	try : 
# 		 #interpolation
# 		 f=interpolate.interp1d(np.squeeze(Vk), np.squeeze(intens))	 
# 		 g=interpolate.interp1d(np.squeeze(Vk), np.squeeze(sigm))	
# 		
# 		 y = f(V)
# 		 z=g(V)
# 		 Nintens[:,k]=y
# 		
# 		 Ns[:,k]=1/z
#  	except : 
# 		 #column of nan if interpolation doesn't work
# 		 a = np.empty(len(V))
# 		 a[:]=np.nan
# 		 Nintens[:,k]=a
# 		 Ns[:,k]=a

# #Equation of WLA
# multintens=Ns*Nintens
# SommeB=np.nansum(multintens,axis=1)
# sommeW=np.nansum(Ns,axis=1)

# Wla = SommeB/sommeW
# ################### V





#get data from ascii
data2 = np.genfromtxt('UMon-Feb2021-NeoNarval-V/UMon_2021-02-14-stokesV.ascii')
dataN = np.genfromtxt('UMon-Feb2021-NeoNarval-V/UMon_2021-02-14-stokesN.ascii')


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


Inconnue = len(np.argwhere(lande==99))
matrice_lande = np.zeros((len(np.argwhere(lande==99)),len(lande)))
param = []
param_libre0 = []
j=0
for k in range (len(mask)) :
	if lande[k] ==99 : 
		param_libre0.append(-1)
		param.append(np.nan)
		matrice_lande[j,k]=1
		j+=1
		
	else : 
		param.append(lande[k])
		matrice_lande[:,k]=np.nan*np.ones(Inconnue)
		


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
		
	

res = least_squares(integrale, param_libre0, args=(V,NV,NsV,param),bounds = (-np.ones(Inconnue),np.ones(Inconnue)),verbose=2)

param_LS = res.x

param_fin = param_LS@matrice_lande
param=np.array(param)
param_fin = np.array(param_fin)
lande_fin = np.nansum(np.dstack((param,param_fin)),2)

NV_fin = NV*lande_fin

mult_LS=NsV*NV_fin
Somme_LS=np.nansum(mult_LS,axis=1)
sommek_LS=np.nansum(NsV,axis=1)
WLA_LS = Somme_LS/(sommek_LS)


NN_fin = NN*lande_fin

multN_LS=NsN*NN_fin
SommeN_LS=np.nansum(multN_LS,axis=1)
sommekN_LS=np.nansum(NsN,axis=1)
WLAN_LS = SommeN_LS/(sommekN_LS)


#Equation of WLA for Stokes V after all the +1/-1
multV=NsV*NV
SommeV=np.nansum(multV,axis=1)
sommekV=np.nansum(NsV,axis=1)

Wla_with_lande_plus_or_minus = SommeV/(sommekV)



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
ax2.plot(V,WLAN_LS,'b',label='Stokes N WLA with Least Square. Date : 2022-05-26')
ax2.plot(V,WLA_LS,'r',label='Stokes V WLA with Least Square Date : 2022-05-26')
plt.xlabel('v(km/s)')
plt.ylabel('V/Ic') 	
plt.vlines(V[indice_min], np.min(Wla_with_N),np.max(Wla_with_N),linestyles='--',color='green')
plt.vlines(V[indice_max], np.min(Wla_with_N),np.max(Wla_with_N),linestyles='--',color='green')
plt.legend()
plt.title('Comparaison Stokes N et V')
plt.grid(True)
plt.show()

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	