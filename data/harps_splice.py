import numpy as np
from astropy.io import fits
import sys
from scipy.interpolate import interp1d
import os

f = sys.argv[1]

def splice_spectrum(fi):
    ff = fits.open(fi)
    # Determine Stokes parameter and target name
    star = ff[0].header['OBJECT'].replace(" ", "")
    # Load data: wave, intensity, polarisation, error, null
    w,s,p = ff[1].data['WAVE'], ff[1].data['INTENS_NORM'], ff[1].data['STOKES/Ic']
    e,n = ff[1].data['ERR_INTENS_NORM'],ff[1].data['NULL/Ic']

    w1 = w  
    s1 = s
    p1 = p
    e1 = e   
    n1 = n

    nwl = len(w1) # number of wavelength points  
    # array of delta wave. if negative: signified switch to another order
    dw = w1[1:nwl] -w1[0:nwl-1]
    ii= np.where(dw < 0)[0] # finding where splicing need to be done
    ni = len(ii) # number of splicing need to be done

    w2=w1[0:ii[0]] # isolate first order
    s2=s1[0:ii[0]]
    e2=e1[0:ii[0]]
    p2=p1[0:ii[0]] - np.nanmedian(p1[0:ii[0]])
    n2=n1[0:ii[0]] - np.nanmedian(n1[0:ii[0]])

    for i in range(ni):
        try:
            i1=ii[i]+1
            if i < ni -1: i2 = ii[i+1]
            else: i2 = nwl-1
            w1x=w1[i1:i2]
            s1x=s1[i1:i2]
            e1x=e1[i1:i2]
            p1x=p1[i1:i2] - np.nanmedian(p1[i1:i2])
            n1x=n1[i1:i2] - np.nanmedian(n1[i1:i2])

            wmax=np.nanmax(w2)
            wmin=np.nanmin(w1x)
            i1o=np.where(w2 > wmin) # where blue order overlaps with red
            i2n=np.where(w1x > wmax)
            f=(wmax-w2[i1o])/(wmax-wmin)
            fs1i=interp1d(w1x,s1x)
            fe1i=interp1d(w1x,e1x)
            fp1i=interp1d(w1x,p1x)
            fn1i=interp1d(w1x,n1x)

            s1i = fs1i(w2[i1o])
            e1i = fe1i(w2[i1o])
            p1i = fp1i(w2[i1o])
            n1i = fn1i(w2[i1o])

            s2[i1o]= s2[i1o]*f+s1i*(1.-f)
            p2[i1o]= p2[i1o]*f+p1i*(1.-f)
            n2[i1o]= n2[i1o]*f+n1i*(1.-f)
            e2[i1o]= np.sqrt((e2[i1o]*f)**2.+(e1i*(1.-f))**2.)

            w2=np.append(w2,w1x[i2n])
            s2=np.append(s2,s1x[i2n])
            e2=np.append(e2,e1x[i2n])
            p2=np.append(p2,p1x[i2n])
            n2=np.append(n2,n1x[i2n])
        except:
            continue

    kk = np.squeeze(np.argwhere(~np.isnan(w2)))
    w2,s2,e2,p2,n2 = w2[kk],s2[kk],e2[kk],p2[kk],n2[kk]
    kk = np.squeeze(np.argwhere(~np.isnan(s2)))
    w2,s2,e2,p2,n2 = w2[kk],s2[kk],e2[kk],p2[kk],n2[kk]
    kk = np.squeeze(np.argwhere(~np.isnan(e2)))
    w2,s2,e2,p2,n2 = w2[kk],s2[kk],e2[kk],p2[kk],n2[kk]
    kk = np.squeeze(np.argwhere(~np.isnan(p2)))
    w2,s2,e2,p2,n2 = w2[kk],s2[kk],e2[kk],p2[kk],n2[kk]
    kk = np.squeeze(np.argwhere(~np.isnan(n2)))
    w2,s2,e2,p2,n2 = w2[kk],s2[kk],e2[kk],p2[kk],n2[kk]

    data = {'WAVE' : w2}
    data['INTENS'] = s2
    data['ERR'] = e2
    data['STOKES'] = p2
    data['NULL'] = n2
    data['star'] = star
    data['date-obs'] = ff[0].header['DATE-OBS']
    return data

data =  splice_spectrum(f)

fileoutP = data['star'] + '_' + data['date-obs'] + '-POL' + '.ascii'
fileoutI = data['star'] + '_' + data['date-obs'] + '-INTENS.ascii'
fileoutN = data['star'] + '_' + data['date-obs'] + '-NULL.ascii'
dirname = './'
np.savetxt(dirname+fileoutI,
  np.stack((data['WAVE'],data['INTENS'],data['ERR'])).T,
  fmt='%.8f %.8f %.8f')

np.savetxt(dirname+fileoutP,  
  np.stack((data['WAVE'],data['STOKES'],data['ERR'])).T,   
  fmt='%.8f %.8f %.8f')

np.savetxt(dirname+fileoutN,  
  np.stack((data['WAVE'],data['NULL'],data['ERR'])).T,   
  fmt='%.8f %.8f %.8f')

