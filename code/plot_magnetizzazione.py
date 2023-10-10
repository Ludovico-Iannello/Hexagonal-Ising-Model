import pylab as plt
import numpy as np
from scipy.optimize import curve_fit
plt.rc('font',size=16)

nu = 1
b = 1/8
g = 7/4
alfa = 0
beta_c = 1/2*np.log(2+np.sqrt(3))


magn=np.zeros(20000)
N=np.zeros(20000)
magn_2=np.zeros(20000)
M,E=np.loadtxt('data_L20_histo_bassa.dat', unpack=True)

for i in range (0,20000):
    N[i]=i+1
    magn[i]=M[i]


plt.figure(1)
N=np.linspace(1,len(M),len(M))
plt.minorticks_on()
plt.title('Istogramma Magnetizzazione')
plt.subplot(3,1,(1,2))
plt.xlabel('M [a.u.]')
plt.ylabel('Distribuzione di M')
plt.hist(M,2000,density=True,histtype='step')
plt.subplot(3,1,3)
plt.errorbar(N,M,marker=',',linestyle='',label='L=20')
plt.xlabel(r'tempo montecarlo [a.u.]')
plt.legend()
plt.ylabel('M [a.u.]')

plt.figure(2)
magn_2,E=np.loadtxt('data_L20_histo.dat', unpack=True)
N=np.linspace(1,len(magn_2),len(magn_2))
plt.minorticks_on()
plt.title('Istogramma Magnetizzazione')
plt.subplot(3,1,(1,2))
plt.ylabel('Distribuzione di M')
plt.xlabel('M [a.u.]')
plt.hist(magn_2,2000,density=True,histtype='step')
plt.subplot(3,1,3)
plt.errorbar(N,magn_2,marker=',',linestyle='',label='L=20')
plt.xlabel(r'tempo montecarlo [a.u.]')
plt.legend()
plt.ylabel('M [a.u.]')
plt.show()
