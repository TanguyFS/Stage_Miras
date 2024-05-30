import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.constants import c


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
data = np.genfromtxt('UHer_2024-03-13-stokesI.ascii')

#observed data
wavelength = data[:,0]
intensity = data[:,1]
s=data[:,2]
w = []
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

pickle = pickle.load(open('t5000g+3.0.sel_atmo+mol+arbitrarygeff.pickle', 'rb'))

mask = pickle['wave_vac']
poids = pickle['depth']
weightI=pickle['weight_I']
lande=pickle['lande']
Ismolec = pickle['is_molecule?']

V=np.arange(Vmin,Vmax,Vsteps)
i=np.array(i)

# plt.subplot(2,1,2)
# plt.plot(wavelength, s, label = "Erreur")
# plt.xlabel('λ(Å)')
# plt.ylabel('Stokes I') 
# plt.legend()
# plt.grid(True)

# plt.subplot(2,1,1)
# plt.plot(wavelength, i, label = "Signal")
# plt.xlabel('λ(Å)')
# plt.ylabel('Stokes I') 
# plt.legend()
# plt.grid(True)