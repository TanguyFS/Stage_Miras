import numpy as np
import matplotlib.pyplot as plt
import pickle


# constant
c = 3*10**5
seuil = 0.2

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
	if intensity[k]<1-seuil or intensity [k] > 1+seuil :
		Nlambda.append(wavelength[k])
		intens.append(intensity[k])
	else:
		Nlambda.append(wavelength[k])
		intens.append(1)
N = np.concatenate((np.vstack(Nlambda),np.vstack(intens)),axis=1) 


Nob=[]
j=0
k=0

while k < len(N) :

	j=0
	if N[k,1]!=1 :
		while N[k,1]!=1 :
			k+=1
			j+=1
		
		Nob.append(N[k-j:k,:])
		
	else :
		k+=1

Center=[]
for k in (range (len(Nob))) :
		  Nk = Nob[k]
		  Nk = Nk[np.argsort(Nk[:,1])]
		  Center.append(Nk[0,0])

Nv_temp=[]


for k in (range(len(Nob))) : 
	mini = c
	for j in (range(len(mask))) :
		diff  = np.abs(Center[k] - mask[j])
		
		if diff < mini :
			  mini =diff
			  m=j
	
	
	
	Nk = Nob[k]
	Nk[:,0] = c*(Nk[:,0] - mask[m])/mask[m]
	Nv_temp.append(Nk)


Pas = 100//1		
V = np.linspace(-50,50,int(Pas))
Nv=[]
			
for k in (range(len(Nv_temp))):	
	Nk = Nv_temp[k]
	intens_temp=[]
	Vtemp=[]
	
	Nsort = Nk[np.argsort(Nk[:,0])]
	Nk = Nv_temp[k]

	maxi = Nsort[-1,0]
	mini = Nsort[0,0]
	y=0
	j=0
	for j in (range(len(V))) : 
		if V[j]<mini:
			Vtemp.append(V[j])
			intens_temp.append(1)
			y+=1
			
		elif V[j]>maxi : 
			Vtemp.append(V[j])
			intens_temp.append(1)
			y+=1

		else : 
			Vtemp.append(Nk[j-y,0])
			intens_temp.append(Nk[j-y,1])
		
	
		
	
		
	Ntemp = np.concatenate((np.vstack(Vtemp),np.vstack(intens_temp)),axis=1) 
	Nv.append(Ntemp)

	
Somme = np.mean(Nv,axis=0)

plt.plot(Somme[:,0],Somme[:,1])
	
	
	
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	