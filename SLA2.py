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

Nlambda=[]
intens=[]


for k in range(len(wavelength)):
	if intensity[k]<0.90 :
		Nlambda.append(wavelength[k])
		intens.append(intensity[k])
	else:
		Nlambda.append(wavelength[k])
		intens.append(1)
N = np.concatenate((np.vstack(Nlambda),np.vstack(intens)),axis=1) 



center=[]
deep=[]



for k in range(len(N)-10):

	if N[k,1]!=1 :
		min = 1
		for j in range (k-10,k+10) :
			if N[j,1]<min :
				min = N[j,1]
				m=j
		center.append(N[m,0])
		deep.append(N[m,1])
				
		


V = []
Depth=[]

for j in range (len(center)) :
	v_min = c
	for k in range(len(mask)) :
		delta_v = c*(center[j]-mask[k])/mask[k]

		if np.abs(delta_v) < np.abs(v_min) : 
			v_min = delta_v
			m=k
	V.append(v_min)
	Depth.append(poids[m])
		
C = np.concatenate((np.vstack(V),np.vstack(Depth)),axis=1) 

C=C[np.argsort(C[:,0])]


a=[]
b=[]
d=[]

for k in range(len(C)//2) :
	a.append((C[2*k,0] + C[2*k+1,0])/2)
	b.append((C[2*k,1] + C[2*k+1,1])/2)
	d.append((deep[2*k] + deep[2*k+1])/2)
	




a.append(C[-1,0]) 
b.append(C[-1,1])
d.append(deep[-1])






C = np.concatenate((np.vstack(a),np.vstack(b)),axis=1) 




C = np.unique(C, axis=0)



print(C)
C[:,1] = 1-C[:,1]
plt.scatter(C[:,0],C[:,1],linewidth=0.01)
plt.show()


















