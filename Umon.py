import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.constants import c
from scipy import interpolate

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
data = np.genfromtxt('UMON_2006-02-09-stokesI.ascii')

#observed data
wavelength = data[:,0]
intensity = data[:,1]

i=[]
Nintens=[]

for k in (range (len(wavelength))):
	
	if wavelength[k]> wo2a[0] and wavelength[k] < wo2a[1]:
		i.append(np.nan) 
	elif wavelength[k]> wo2b[0] and wavelength[k] < wo2b[1]:
		i.append(np.nan) 
	elif wavelength[k]> wH2oa[0] and wavelength[k] < wH2oa[1]:
		i.append(np.nan)
	elif wavelength[k]> wH2ob[0] and wavelength[k] < wH2ob[1]:
		i.append(np.nan)
	elif wavelength[k]> whugg[0] and wavelength[k] < whugg[1]:
			i.append(np.nan)
	else :
		i.append(intensity[k]) 



V=np.arange(Vmin,Vmax,Vsteps)
i=np.array(i)



#get mask
pickle = pickle.load(open('t5000g+3.0.sel.pickle', 'rb'))

mask = pickle['wave_vac']
poids = pickle['depth']

for k in (range (len(mask))) :
	wmin= mask[k]*Vmin/c + mask[k]
	wmax= mask[k]*Vmax/c + mask[k]
	iii=np.argwhere(np.logical_and(wavelength > wmin-0.2, wavelength < wmax+0.2))
	w=wavelength[iii]
	intens = i[iii]	
	Vk = c*(w-mask[k])/mask[k]
	try : 
		f=interpolate.interp1d(np.squeeze(Vk), np.squeeze(intens))	 
		
		y = f(V)
		Ntemp = np.concatenate((np.vstack(V),np.vstack(y)),axis=1) 
		
		Nintens.append(Ntemp)
	except : 
		a = np.empty((len(V),2))
		a[:]=np.nan
		Nintens.append(a)

Somme = np.nanmean(Nintens,axis=0)




################### POL




#get data from ascii
data2 = np.genfromtxt('UMON_2006-02-09-stokesQ.ascii')

#observed data
wavelength2= data2[:,0]
intensity2 = data2[:,1]

NQ=[]




for k in (range (len(mask))) :
	wmin= mask[k]*Vmin/c + mask[k]
	wmax= mask[k]*Vmax/c + mask[k]
	iii=np.argwhere(np.logical_and(wavelength2 > wmin-0.2, wavelength2 < wmax+0.2))
	w=wavelength2[iii]
	intens = intensity2[iii]	
	Vk = c*(w-mask[k])/mask[k]
	try : 
		f=interpolate.interp1d(np.squeeze(Vk), np.squeeze(intens))	 
		
		y = f(V)
		Ntemp = np.concatenate((np.vstack(V),np.vstack(y)),axis=1) 
		
		NQ.append(Ntemp)
	except : 

		NQ.append(a)


Somme2 = np.nanmean(NQ,axis=0)


fig=plt.figure()

ax1 = plt.subplot(2,1,1)
ax1.plot(Somme[:,0],Somme[:,1],label='Stokes I', color='forestgreen')
plt.xlabel('v(km/s)')
plt.ylabel('I/Ic') 	
plt.legend()
plt.grid(True)


ax2=plt.subplot(2,1,2)
ax2.plot(Somme2[:,0],Somme2[:,1],label='Stokes Q', color='darkred')
plt.xlabel('v(km/s)')
plt.ylabel('Q/Ic') 
plt.legend()
plt.grid(True)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	