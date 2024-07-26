import numpy as np
import matplotlib.pyplot as plt
import pickle
import seaborn as sns
import plotly.graph_objects as go
import plotly.subplots as sp
import plotly.io as pio

pio.renderers.default = 'pdf'

pick = pickle.load(open('t5000g+3.0.sel_atmo+mol.pickle', 'rb'))
lande=pick['lande'] #lande factor


index= np.argwhere(lande==99)

data1 = pickle.load(open('data/UMon_2021-02-14-1.pickle', 'rb'))
lande1=data1['Lande'] #lande factor

data2 = pickle.load(open('data/UMon_2021-02-24-1.pickle', 'rb'))
lande2=data2['Lande'] #lande factor

data3 = pickle.load(open('data/UMon_2021-02-25-stokesV.pickle', 'rb'))
lande3=data3['Lande'] #lande factor

data4 = pickle.load(open('data/UMon_2021-02-26-stokesV.pickle', 'rb'))
lande4=data4['Lande'] #lande factor

data5 = pickle.load(open('data/UMon_2021-02-27-stokesV.pickle', 'rb'))
lande5=data5['Lande'] #lande factor

data6 = pickle.load(open('data/UMon_2021-02-stokesV-SUM.pickle', 'rb'))
lande6=data6['Lande'] #lande factor




fig = sp.make_subplots(191,10)

for k in range(len(index)) : 
	landek = []
	landek.append(lande1[index[k][0]])
	landek.append(lande2[index[k][0]])
	landek.append(lande3[index[k][0]])
	landek.append(lande4[index[k][0]])
	landek.append(lande5[index[k][0]])
	landek.append(lande6[index[k][0]])

	row = k // 10 + 1
	col = k%10 +1
	fig.add_trace(go.Histogram(x=landek), row=row, col=col)
	

fig.update_layout(height=2000, width=1500, title_text="Histograms")
fig.show()	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	