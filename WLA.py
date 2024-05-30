import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.constants import c
from scipy import interpolate
import math
# constant
Vmax = 200
Vmin = -200
Vsteps = 3.7
whugg=[3600,4000]
wo2b = [6866,6950]
wo2a= [7590,7724]
wH2oa=[8950,9860]
wH2ob=[8100,8300]
c = c/1000

#get data from ascii
data = np.genfromtxt('UMon_2014-04-11-stokesI.ascii')

#observed data
wavelength = data[:,0]
intensity = data[:,1]
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
pickle = pickle.load(open('t5000g+3.0.sel_atmo+mol+arbitrarygeff.pickle', 'rb'))

mask = pickle['wave_vac']
poids = pickle['depth']


Nintens=np.zeros((len(V),len(mask)))
Ns=np.zeros((len(V),len(mask)))


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

Somme = np.nanmean(Nintens,axis=1)

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
Ns2=np.zeros((len(V),len(mask)))





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
		
		Ns2[:,k]=1/z
	except : 
		a = np.empty(len(V))
		a[:]=np.nan
		NQ[:,k]=a
		Ns2[:,k]=a

Somme2 = np.nanmean(NQ,axis=1)

multkkk=Ns2*NQ
Sommekk=np.nansum(multkkk,axis=1)
sommekkk=np.nansum(Ns2,axis=1)

Wla2 = Sommekk/sommekkk



####### Barre d'erreurs





We=Ns2
Se=1/Ns2
erreur = np.sqrt(np.nansum(We*Se*Se,axis=1)/(np.nansum(We,axis=1)*np.shape(NQ)[1]))




























fig=plt.figure()

ax1 = plt.subplot(2,1,1)
ax1.plot(V,Somme,label='Stokes I SLA',color='forestgreen')
ax1.plot(V,Wla,label='Stokes I WLA', color='darkred')
plt.xlabel('v(km/s)')
plt.ylabel('I/Ic') 	
plt.legend()
plt.grid(True)



ax2 = plt.subplot(2,1,2)
ax2.plot(V,Somme2,label='Stokes V SLA',color='forestgreen')
ax2.plot(V,Wla2,label='Stokes V WLA', color='darkred')
plt.xlabel('v(km/s)')
plt.ylabel('I/Ic') 	
plt.legend()
plt.grid(True)



# ax2=plt.subplot(1,1,1)
# ax2.plot(V,Wla2,label='Stokes V WLA', color='darkred')
# ax2.errorbar(V,Wla2,yerr=erreur,barsabove=True,ecolor="blue")

# plt.xlabel('v(km/s)')
# plt.ylabel('V/Ic') 
# plt.legend()
# plt.grid(True)
 	
	
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	