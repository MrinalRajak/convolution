

# convolution

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps,quad

def f(t):
    return A/np.sqrt(2*np.pi*sigma1**2)*np.exp(-(t-mu1)**2/(2*sigma1**2))
def g(t):
    return B/np.sqrt(2*np.pi*sigma2**2)*np.exp(-(t-mu2)**2/(2*sigma2**2))

A,mu1,sigma1=1.5,-2.0,1.3   #parameters for 1
B,mu2,sigma2=1.5,2.0,1.3    #parameters for 2
Fs=15.0                     #sampling frequency

T=30.0                      #Time range we are interested
t=np.arange(-T,T,1.0/Fs)
# calculation of overall convolution result
# using simpsons integration
I=[]
for t1 in t:
    prod=lambda tau:f(tau)*g(t1-tau)
    I.append(simps(prod(t),t))

fig,(ax0,ax1)=plt.subplots(nrows=2)
ax0.plot(t,f(t),lw=2,ls='--',c='orange',label='Function f(t)')
ax0.plot(t,g(t),lw=2,ls='-.',c='blue',label='Function g(t)')
ax1.plot(t,I,c='red',label='Integrated data')
cc=np.convolve(f(t),g(t),mode='Same')-(1./Fs)
ax1.plot(t,cc,'o',label='convolused data')
ax0.legend()
ax1.legend()
ax0.set_title('convolution')
plt.show()
         
