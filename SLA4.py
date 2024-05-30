import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.constants import c
from scipy import interpolate

# constant
Vmax = 100
Vmin = -100
Vsteps = 1
w02 = [6550,6575]
wtell = [6867,6913]
wint1=[6759.5,6762]
wint2=[5338,5339]
c = c/1000

#get data from ascii
data = np.genfromtxt('HD-201601_2012-07-16T07:19:11.908-INTENS.ascii')

#observed data
wavelength = data[:,0]
intensity = data[:,1]

i=[]
Nintens=[]

for k in (range (len(wavelength))):
	
	if wavelength[k]> w02[0] and wavelength[k] < w02[1]:
		i.append(np.nan) 
	elif wavelength[k]> wtell[0] and wavelength[k] < wtell[1]:
		i.append(np.nan) 
	elif wavelength[k]> wint1[0] and wavelength[k] < wint1[1]:
		i.append(np.nan)
	elif wavelength[k]> wint2[0] and wavelength[k] < wint2[1]:
		i.append(np.nan)
	else :
		i.append(intensity[k]) 



V=np.arange(Vmin,Vmax,Vsteps)
i=np.array(i)



#get mask
pickle = pickle.load(open('t8000g+4.5.sel_arbitrary-geff.pickle', 'rb'))

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
		a = np.empty(((int((Vmax-Vmin)/Vsteps)),2))
		a[:]=np.nan
		Nintens.append(a)

Somme = np.nanmean(Nintens,axis=0)
Somme3=np.nanmedian(Nintens,axis=0)




################### POL




#get data from ascii
data2= np.genfromtxt('HD-201601_2012-07-16T07:19:11.908-POL.ascii')

#observed data
wavelength2= data2[:,0]
intensity2 = data2[:,1]

Npol=[]




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
		
		Npol.append(Ntemp)
	except : 

		Npol.append(a)


Somme2 = np.nanmean(Npol,axis=0)


fig=plt.figure()

ax1 = plt.subplot(2,1,1)
ax1.plot(Somme[:,0],Somme[:,1],label='Stokes I', color='forestgreen')
plt.xlabel('v(km/s)')
plt.ylabel('I/Ic') 	
plt.grid(True)


ax2=plt.subplot(2,1,2)
ax2.plot(Somme3[:,0],Somme3[:,1],label='Stokes I median', color='darkred')
plt.xlabel('v(km/s)')
plt.ylabel('V/Ic') 
plt.legend()
plt.grid(True)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	