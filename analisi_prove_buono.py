import pylab as plt
import numpy as np
from scipy.optimize import curve_fit

'''
beta=np.zeros(50)
for ibeta in range (0,50):
    beta[ibeta]=0.4+((0.9-0.4)/50)*(ibeta+1)
'''
nu = 1
b = 1/8
g = 7/4
alfa = 0
beta_c = 1/2*np.log(2+np.sqrt(3))

for i in range (1,6):

    beta,M,dM,E,dE,X,dX,C,dC,B,dB=np.loadtxt('ising_L'+str(10*i)+'.dat', unpack=True)
    L=i*10

#magnetizzzione
    plt.figure(1)
    plt.minorticks_on()
    plt.title('Modulo magnetizzazione intorno alla transizione per diversi volumi')
    plt.errorbar(beta,M,dM,marker='.',linestyle='--',lw=0.3,capsize=4,label='L='+str(i*10))
    plt.legend()
    plt.xlabel(r'$\beta$ [a.u.]')
    plt.ylabel('|M| [a.u.]')
#calore specifico
    plt.figure(2)
    plt.minorticks_on()
    plt.title('Calore specifico intorno alla transizione per diversi volumi')
    plt.errorbar(beta,C,dC,marker='.',linestyle='--',lw=0.3,capsize=4,label='L='+str(i*10))
    plt.legend()
    plt.xlabel(r'$\beta$ [a.u.]')
    plt.ylabel(r'$C_v$ [a.u.]')
#suscettività
    plt.figure(3)
    plt.minorticks_on()
    plt.title('Suscettività magnetica intorno alla transizione per diversi volumi')
    plt.errorbar(beta,X,dX,marker='.',linestyle='--',lw=0.3,capsize=4,label='L='+str(i*10))
    plt.legend()
    plt.xlabel(r'$\beta$ [a.u.]')
    plt.ylabel(r'$\chi$ [a.u.]')
##cumulanti
    plt.figure(4)
    plt.minorticks_on()
    plt.title('Cumulante di binder intorno alla transizione per diversi volumi')
    plt.errorbar(beta,B,dB,marker='.',linestyle='--',lw=0.3,capsize=4,label='L='+str(i*10))
    plt.legend()
    plt.xlabel(r'$\beta$ [a.u.]')
    plt.ylabel('Cumulante di Binder [a.u.]')

#FSS suscettività sovrapposte
    plt.figure(5)
    plt.minorticks_on()
    plt.title('Suscettività magnetica F.S.S.')
    plt.errorbar((beta-beta_c)*L**(1/nu),X/(L**(g/nu)),dX/(L**(g/nu)),marker='.',
                 linestyle='--',lw=0.3,capsize=4,label='L='+str(i*10))
    plt.legend()
    plt.xlabel(r'$(\beta-\beta_c)L^{1/\nu}$ [a.u.]')
    plt.ylabel(r'$\chi/L^{\gamma/\nu}$ [a.u.]')

#FSS magnetizzazioni sovrapposte
    if i>1:
        plt.figure(6)
        plt.minorticks_on()
        plt.title('Modulo magnetizzazione F.S.S')
        plt.errorbar((beta-beta_c)*L**(1/nu),M/(L**(-b/nu)),dM/(L**(-b/nu)),marker='.',linestyle='--',
                     lw=0.3,capsize=4,label='L='+str(i*10))
        plt.legend()
        plt.xlabel(r'$(\beta-\beta_c)L^{1/\nu}$ [a.u.]')
        plt.ylabel(r'$|M|/L^{-\beta/\nu}$ [a.u.]')

##fit al variare di L
M=np.zeros((50,i))
dM=np.zeros((50,i))
E=np.zeros((50,i))
dE=np.zeros((50,i))
X=np.zeros((50,i))
dX=np.zeros((50,i))
C=np.zeros((50,i))
dC=np.zeros((50,i))
B=np.zeros((50,i))
dB=np.zeros((50,i))

## array dei valori
L=np.zeros(i)
for i in range (0,len(L)):
    beta,M[:,i],dM[:,i],E[:,i],dE[:,i],X[:,i],dX[:,i],C[:,i],dC[:,i],B[:,i],dB[:,i]=np.loadtxt(
        'ising_L'+str(10*(i+1))+'.dat', unpack=True)
    L[i]=(i+1)*10

## fit beta pseudo critico
beta_pse=np.zeros(5)
for i in range (0,5):
    beta_pse[i]=beta[np.argmax(X[:,i])]

def fit(L,beta_c,xbar,nu):
    return beta_c+xbar*L**(-1/nu)

dbeta_pse=beta_pse*0 + (0.8-0.5)/50/np.sqrt(12)/2
pars, covm = curve_fit(fit, L, beta_pse, sigma=dbeta_pse,absolute_sigma=True)
err = np.sqrt(covm.diagonal())
chisq = np.sum(((beta_pse - fit(L, *pars))/dbeta_pse)**2)
ndof = len(L) - len(pars)

plt.figure(7)
xx=np.linspace(min(L),max(L),1000)
plt.minorticks_on()
plt.title(r'fit $\beta_{pc}$ al variare di L')
plt.errorbar(L,beta_pse,dbeta_pse,marker='.',capsize=4,linestyle='')
plt.plot(xx,fit(xx,*pars))
plt.xlabel(r'L [a.u.]')
plt.ylabel(r'$\beta_{pc}$ [a.u.]')

print('chi quadro = %.3f (%d dof)' % (chisq, ndof))
print('parametri fit: beta_c, xbar, nu')
print(pars)
print('errori fit')
print(err)

## fit bilog per determinare g/nu picco suscettività in funzione di L

chi_max=np.zeros(5)
dchi_max=np.zeros(5)
for i in range (0,5):
    chi_max[i]=max(X[:,i])
    dchi_max[i]=dX[np.argmax(X[:,i]),i]

def fit_2(L,g_su_nu,c):
    return c*L**(g_su_nu)

pars, covm = curve_fit(fit_2, L, chi_max, sigma=dchi_max,absolute_sigma=True)
err = np.sqrt(covm.diagonal())
chisq = np.sum(((chi_max - fit_2(L, *pars))/dchi_max)**2)
ndof = len(L) - len(pars)

plt.figure(8)
xx=np.linspace(min(L),max(L),10000)
plt.minorticks_on()
plt.title(r'fit $\chi_m$ al variare di L')
plt.errorbar(L,chi_max,dchi_max,marker='.',capsize=4,linestyle='')
plt.plot(xx,fit_2(xx,*pars))
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'L [a.u.]')
plt.ylabel(r'$\chi_{m}$ [a.u.]')

print('chi quadro = %.3f (%d dof)' % (chisq, ndof))
print('parametri fit: g/nu,cost')
print(pars)
print('errori fit')
print(err)


plt.show()
