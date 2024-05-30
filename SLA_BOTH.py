
import numpy as np
import matplotlib.pyplot as plt
import pickle

#prendre quelque nanomètre en plus

# constant
c = 3*10**5
seuil = 0.02
seuil2 = 0.07

Reponse = input('Que voulez vous voir ? (Pol, Intens, Both)\n')


if Reponse =='Intens' :

	#get data from ascii
	data = np.genfromtxt('HD-201601_2012-07-16T07:19:11.908-INTENS.ascii')
	
	#observed data
	wavelength = data[:,0]
	intensity = data[:,1]
	
	#get mask
	pickle = pickle.load(open('t8000g+4.5.sel_arbitrary-geff.pickle', 'rb'))
	
	#list with the data of these masks
	mask = pickle['wave_vac']
	poids = pickle['depth']
	
	
	#temporary vector
	Nlambda=[]
	intens=[]
	
	#Smoothing up the data to make the ray selection easier
	for k in range(len(wavelength)):
		if intensity[k]<1-seuil or intensity [k] > 1+seuil :
			Nlambda.append(wavelength[k])
			intens.append(intensity[k])
		else:
			Nlambda.append(wavelength[k])
			intens.append(1)
	N = np.concatenate((np.vstack(Nlambda),np.vstack(intens)),axis=1) 
	
	
	#make an array of all the rays left
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
	
	#find the center of those rays
	Center=[]
	for k in (range (len(Nob))) :
			  Nk = Nob[k]
			  Nk = Nk[np.argsort(Nk[:,1])]
			  Center.append(Nk[0,0])
	
	
	#find the closer ray in the mask and converte wavelength in velocity
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
	
	
	
	#creating our grid of velocity
	Pas =1
	V = np.linspace(-100,100,200)
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
	
	plt.plot(Somme[:,0],Somme[:,1],color ='forestgreen',label='Stokes I')
		
	plt.xlabel('v(km/s)')
	plt.ylabel('I/Ic') 	
		
		
		
if Reponse =='Pol'	: 
	
	
	#get data from ascii
	data2= np.genfromtxt('HD-201601_2012-07-16T07:19:11.908-POL.ascii')
	
	#observed data
	wavelength2= data2[:,0]
	intensity2 = data2[:,1]*3+1
	
	#get mask
	pickle = pickle.load(open('t8000g+4.5.sel_arbitrary-geff.pickle', 'rb'))
	
	mask = pickle['wave_vac']
	poids = pickle['depth']
	
	Nlambda=[]
	intens=[]
	
	
	
	for k in range(len(wavelength2)):
		if intensity2[k]<1-seuil2 or intensity2 [k] > 1+seuil2 :
			Nlambda.append(wavelength2[k])
			intens.append(intensity2[k])
		else:
			Nlambda.append(wavelength2[k])
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
	
	
	Pas = 200//1		
	V = np.linspace(-100,100,int(Pas))
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
	
		
	Somme2 = np.mean(Nv,axis=0)
	
	plt.plot(Somme2[:,0],Somme2[:,1],color='darkred',label='Stokes V')
		
	plt.xlabel('v(km/s)')
	plt.ylabel('V/Ic') 
	
		
if Reponse == 'Both' :	
		
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
	
		
		
		
	
	seuil2 = 0.08

	#get data from ascii
	data2= np.genfromtxt('HD-201601_2012-07-16T07:19:11.908-POL.ascii')
	
	#observed data
	wavelength2= data2[:,0]
	intensity2 = data2[:,1]*3+1
	
	Nlambda=[]
	intens=[]
	
	
	
	for k in range(len(wavelength2)):
		if intensity2[k]<1-seuil2 or intensity2 [k] > 1+seuil2 :
			Nlambda.append(wavelength2[k])
			intens.append(intensity2[k])
		else:
			Nlambda.append(wavelength2[k])
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
	
		
	Somme2 = np.mean(Nv,axis=0)
	
	plt.plot(Somme2[:,0],0.1+(Somme2[:,1]*10+1-10),color='darkred', label = 'Stokes V')
	plt.plot(Somme[:,0],Somme[:,1],color='forestgreen',label='Stokes I')
			
		
	plt.xlabel('v(km/s)')
	plt.ylabel('I                                  V(scaled 30x') 
	

plt.grid(True) 
plt.legend()		
		
		
		
		
		
		
		
		
		
		
		
		