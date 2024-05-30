import numpy as np
import matplotlib.pyplot as plt
import pickle

# constant
c = 3*10**5
delta_nu = 0.6


#get data from ascii
data = np.genfromtxt('HD-201601_2012-07-16T07:19:11.908-INTENS.ascii')

#observed data
wavelength = data[:,0]
intensity = data[:,1]

#get mask
pickle = pickle.load(open('t8000g+4.5.sel_arbitrary-geff.pickle', 'rb'))

mask = pickle['wave_vac']
poids = pickle['depth']


#keep  useful mask
#lambd0=[]
#for k in range (len(mask)):
#	if (mask[k] > wavelength[0] and mask[k] < wavelength[-1]) :
#		lambd0.append(mask[k])


lambd0 = np.array(mask)


delta_v=[]

Num_raie =[]

intens_raie=[]

k=0
i=0
for j in range (len(lambd0)):
	k=0
	while (wavelength[k]<lambd0[j]+1 and k< len(wavelength)-1):
		if np.logical_and(wavelength[k] > lambd0[j] - 0.3,wavelength[k]< lambd0[j] + 0.3):
			Num_raie.append(i)
			delta_v.append(c*(wavelength[k]-lambd0[j])/lambd0[j])
			intens_raie.append(intensity[k]) 
			weight = poids[j]
		k+=1

	i+=1





Nbis= np.concatenate((np.vstack(Num_raie),np.vstack(delta_v)),axis=1)
N = np.concatenate((Nbis,np.vstack(intens_raie)),axis=1) 


for k in range(len(N)//70):
	A = N[2*k*35:35*(2*k+1),1:3]
	B = N[35*(2*k+1):35*(2*k+2),1:3]
	
	C = (A+B)/2
	

	
C=1-C[np.argsort(C[:,0])]
plt.plot(C[:,0],C[:,1])
plt.show()
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

