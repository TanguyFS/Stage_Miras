import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.constants import c
from scipy import interpolate
import math
from scipy import integrate
# constant
Vmax = 200
Vmin = -200
Vsteps = 3.7
indice_max = 12
whugg=[3600,4000]
wo2b = [6866,6950]
wo2a= [7590,7724]
wH2oa=[8950,9860]
wH2ob=[8100,8300]
c = c/1000

def integrale(V,raie):
	Ntemp = np.abs(raie)
	isole = Ntemp[52:73]
	I1 = integrate.simpson(isole,dx=Vsteps)
	return(I1)



#get data from ascii
data = np.genfromtxt('UMon_2014-04-11-stokesI.ascii')

#observed data
wavelength = data[:,0]
intensity = data[:,1] - 1
sigma=data[:,2]

i=[]

s=[]

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



V=np.arange(Vmin,Vmax,Vsteps)
i=np.array(i)
s=np.array(s)


#get mask
pickle = pickle.load(open('t5000g+3.0.sel_atmo+mol.pickle', 'rb'))

mask = pickle['wave_vac']
poids = pickle['depth']
weightI=pickle['weight_I']
lande=pickle['lande']
Ismolec = pickle['is_molecule?']
Weight_Z=pickle['weight_Zeeman']

Nintens=np.zeros((len(V),len(mask)))
Nintens2=np.zeros((len(V),len(mask)))
Ns=np.zeros((len(V),len(mask)))
Ns1=np.zeros((len(V),len(mask)))

for k in (range (len(mask))) :
	wmin= mask[k]*Vmin/c + mask[k]
	wmax= mask[k]*Vmax/c + mask[k]
	iii=np.argwhere(np.logical_and(wavelength > wmin-0.2, wavelength < wmax+0.2))
	w=wavelength[iii]
	intens = i[iii]
	sigm=s[iii]	
	Vk = c*(w-mask[k])/mask[k]
	try : 
		f=interpolate.interp1d(np.squeeze(Vk), np.squeeze(intens))	 
		g=interpolate.interp1d(np.squeeze(Vk), np.squeeze(sigm))	
		
		y = f(V)
		z=g(V)
		Nintens[:,k]=y
		
		Ns[:,k]=1/z
	except : 
		a = np.empty(len(V))
		a[:]=np.nan
		Nintens[:,k]=a
		Ns[:,k]=a

Sommeintens = np.nanmean(Nintens,axis=1)

multkkk=Ns*Nintens
Sommekk=np.nansum(multkkk,axis=1)
sommekkk=np.nansum(Ns,axis=1)

Wla = Sommekk/sommekkk
################### Q




#get data from ascii
data2 = np.genfromtxt('UMon_2014-04-11-stokesV.ascii')

#observed data
wavelength2= data2[:,0]
intensity2 = data2[:,1]
s2=data2[:,2]


NQ=np.zeros((len(V),len(mask)))
NQ2=np.zeros((len(V),len(mask)))
Ns2=np.zeros((len(V),len(mask)))
Ns3=np.zeros((len(V),len(mask)))





for k in (range (len(mask))) :
	
	wmin= mask[k]*Vmin/c + mask[k]
	wmax= mask[k]*Vmax/c + mask[k]
	iii=np.argwhere(np.logical_and(wavelength > wmin-0.2, wavelength < wmax+0.2))
	w=wavelength2[iii]
	intens = intensity2[iii]
	sigm=s2[iii]	
	Vk = c*(w-mask[k])/mask[k]
	try : 
		f=interpolate.interp1d(np.squeeze(Vk), np.squeeze(intens))	 
		g=interpolate.interp1d(np.squeeze(Vk), np.squeeze(sigm))	
		
		y = f(V)
		z=g(V)
		NQ[:,k]=y
		NQ2[:,k]=y

		Ns2[:,k]=1/z
		Ns3[:,k]=1/z


	except : 
		a = np.empty(len(V))
		a[:]=np.nan
		NQ[:,k]=a
		NQ2[:,k]=a
		Ns2[:,k]=a
		Ns3[:,k]=a

		
	if lande[k]==99 :
		
		NQk = NQ[:,k]
		Ns2k = Ns2[:,k]
		
		NQpositif = NQk
		NQnegatif = -NQk
		Ns2positif = Ns2k
		Ns2negatif = -Ns2k
		NQtemp = NQ[:,0:k+1]
		NQtemp[:,k]=NQpositif	
		Ns2temp = Ns2[:,0:k+1]
		Ns2temp[:,k]=NQpositif	
		mult=Ns2temp*NQtemp
		Somme=np.nansum(mult,axis=1)
		sommek=np.nansum(Ns2temp,axis=1)
		Wla2positif = Somme/sommek
		Ipositif = integrale(V, Wla2positif)
		
		
		NQtemp = NQ[:,0:k+1]
		NQtemp[:,k]=NQnegatif	
		Ns2temp = Ns2[:,0:k+1]
		Ns2temp[:,k]=NQnegatif	
		mult=Ns2temp*NQtemp
		Somme=np.nansum(multkkk,axis=1)
		sommek=np.nansum(Ns2temp,axis=1)
		Wla2negatif = Somme/sommek
		Inegatif = integrale(V, Wla2negatif)
		
		if Ipositif < Inegatif : 
			NQ[:,k]=NQnegatif
			Ns2[:,k]=Ns2negatif
		else : 
			NQ[:,k]=NQpositif
			Ns2[:,k]=Ns2positif
	else : 
		
		NQ[:,k] = NQ[:,k]*lande[k]
		
Sommelande = np.sum(lande)

Somme2 = np.nanmean(NQ2,axis=1)
mult3=Ns3*NQ2
Somme3=np.nansum(mult3,axis=1)
sommek3=np.nansum(Ns3,axis=1)

Wla_without_lande = Somme3/sommek3


mult=Ns3*NQ
Somme=np.nansum(mult,axis=1)
sommek=np.nansum(Ns2,axis=1)

Wla_with_lande = Somme/(sommek)



####### Barre d'erreurs





We=Ns2
Se=1/Ns2
erreur = np.sqrt(np.nansum(We*Se*Se,axis=1)/(np.nansum(We,axis=1)*np.shape(NQ)[1]))







# fig=plt.figure()

# ax1 = plt.subplot(2,1,1)
# ax1.plot(V,Sommeintens+1.1,label='Stokes I SLA',color='forestgreen')
# ax1.plot(V,Wla+1.05,label='Stokes I WLA', color='darkred')
# plt.xlabel('v(km/s)')
# plt.ylabel('I/Ic') 	
# plt.legend()
# plt.grid(True)



ax2 = plt.subplot(1,1,1)
ax2.plot(V,Wla_without_lande,label='Stokes V WLA without Lande',color='forestgreen')
ax2.plot(V,Wla_with_lande,label='Stokes V WLA with Lande', color='darkred')
plt.xlabel('v(km/s)')
plt.ylabel('I/Ic') 	
plt.legend()
plt.grid(True)

plt.show()


# ax2=plt.subplot(1,1,1)
# ax2.plot(V,Wla2,label='Stokes V WLA', color='darkred')
# ax2.errorbar(V,Wla2,yerr=erreur,barsabove=True,ecolor="blue")

# plt.xlabel('v(km/s)')
# plt.ylabel('V/Ic') 
# plt.legend()
# plt.grid(True)
 	
	
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	