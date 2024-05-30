import numpy as np
import matplotlib.pyplot as plt
import pickle
from math import *

# constant
c = 3*10**5
Vmax = 100
Vmin = -100
Vsteps = 1
w02 = [6550,6575]
wtell = [6867,6913]
wint1=[6759.5,6762]
wint2=[5338,5339]


#get data from ascii
data = np.genfromtxt('HD-201601_2012-07-16T07:19:11.908-INTENS.ascii')

#observed data
wavelength = data[:,0]
intensity = data[:,1]

w = []
i=[]

for k in (range (len(wavelength))):
	
	if wavelength[k]> w02[0] and wavelength[k] < w02[1]:
		i.append(nan) 
	elif wavelength[k]> wtell[0] and wavelength[k] < wtell[1]:
		i.append(nan) 
	elif wavelength[k]> wint1[0] and wavelength[k] < wint1[1]:
		i.append(nan)
	elif wavelength[k]> wint2[0] and wavelength[k] < wint2[1]:
		i.append(nan)
	else :
		i.append(intensity[k]) 



V=np.arange(Vmin,Vmax,Vsteps)




#get mask
pickle = pickle.load(open('t8000g+4.5.sel_arbitrary-geff.pickle', 'rb'))

mask = pickle['wave_vac']
poids = pickle['depth']

plt.plot(wavelength, i)
		

