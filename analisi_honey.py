import pylab as plt
import numpy as np
from scipy.optimize import curve_fit

plt.rc('font',size=16)

nu = 1
b = 1/8
g = 7/4
alfa = 0
beta_c = 1/2*np.log(2+np.sqrt(3))

for i in range (1,7):

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

    plt.figure(12)
    plt.minorticks_on()
    plt.title('Calore specifico F.S.S.')
    plt.errorbar((beta-beta_c)*L**(1/nu),C+(2.32058930-max(C)),dC,marker='.',
                 linestyle='--',lw=0.3,capsize=4,label='L='+str(i*10))
    plt.legend()
    plt.xlabel(r'$(\beta-\beta_c)L^{1/\nu}$ [a.u.]')
    plt.ylabel(r'$C/L^{\alpha/\nu}$ + correzioni  [a.u.]')
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

'''
##fit al variare di L
M=np.zeros((50,10))
dM=np.zeros((50,10))
E=np.zeros((50,10))
dE=np.zeros((50,10))
X=np.zeros((50,10))
dX=np.zeros((50,10))
C=np.zeros((50,10))
dC=np.zeros((50,10))
B=np.zeros((50,10))
dB=np.zeros((50,10))

## array dei valori
L=np.zeros(10)
for i in range (0,len(L)+1):
    if i<9:
        beta,M[:,i],dM[:,i],E[:,i],dE[:,i],X[:,i],dX[:,i],C[:,i],dC[:,i],B[:,i],dB[:,i]=np.loadtxt(
        'ising_L'+str(5*(i+2))+'.dat', unpack=True)
        L[i]=(i+2)*5
    if i==10:
        beta,M[:,i-1],dM[:,i-1],E[:,i-1],dE[:,i-1],X[:,i-1],dX[:,i-1],C[:,i-1],dC[:,i-1],B[:,i-1],dB[:,i-1]=np.loadtxt('ising_L'+str(5*(i+2))+'.dat', unpack=True)
        L[i-1]=(i+2)*5

## fit beta pseudo critico
beta_pse=np.zeros(10)
for i in range (0,10):
    beta_pse[i]=beta[np.argmax(X[:,i])]

def fit(L,beta_c,xbar,nu):
    return beta_c+xbar*L**(-1/nu)

dbeta_pse=beta_pse*0 + (0.8-0.5)/50/np.sqrt(12)
pars, covm = curve_fit(fit, L, beta_pse, sigma=dbeta_pse,absolute_sigma=True)
err = np.sqrt(covm.diagonal())
chisq = np.sum(((beta_pse - fit(L, *pars))/dbeta_pse)**2)
ndof = len(L) - len(pars)

plt.figure(7)
xx=np.linspace(min(L),max(L),1000)
plt.minorticks_on()

plt.subplot(3,1,(1,2))
plt.title(r'fit $\beta_{pc}$ al variare di L')
plt.ylabel(r'$\beta_{pc}$ [a.u.]')
plt.errorbar(L,beta_pse,dbeta_pse,marker='.',capsize=4,linestyle='')
plt.plot(xx,fit(xx,*pars))
plt.subplot(3,1,3)
plt.axhline(0)
plt.errorbar(L,((beta_pse - fit(L, *pars))/dbeta_pse),1,marker='.',capsize=4,linestyle='')
plt.xlabel(r'L [a.u.]')
plt.ylabel(r'res.norm')

print('chi quadro = %.3f (%d dof)' % (chisq, ndof))
print('parametri fit: beta_c, xbar, nu')
print(pars)
print('errori fit')
print(err)

## fit bilog per determinare g/nu picco suscettività in funzione di L

chi_max=np.zeros(10)
dchi_max=np.zeros(10)
for i in range (0,10):
    chi_max[i]=max(X[:,i])
    dchi_max[i]=dX[np.argmax(X[:,i]),i]

def fit_log(L,g_su_nu,c):
    return c*L**(g_su_nu)

pars, covm = curve_fit(fit_log, L, chi_max, sigma=dchi_max,absolute_sigma=True)
err = np.sqrt(covm.diagonal())
chisq = np.sum(((chi_max - fit_log(L, *pars))/dchi_max)**2)
ndof = len(L) - len(pars)

plt.figure(8)
xx=np.linspace(min(L),max(L),10000)
plt.minorticks_on()

plt.subplot(3,1,(1,2))
plt.title(r'fit $\chi_m$ al variare di L')
plt.errorbar(L,chi_max,dchi_max,marker='.',capsize=4,linestyle='')
plt.plot(xx,fit_log(xx,*pars))
plt.xscale('log')
plt.yscale('log')
plt.ylabel(r'$\chi_{m}$ [a.u.]')
plt.subplot(3,1,3)
plt.errorbar(L,((chi_max - fit_log(L, *pars))/dchi_max),1,marker='.',capsize=4,linestyle='')
plt.axhline(0)
plt.xlabel(r'L [a.u.]')
plt.xscale('log')
plt.ylabel(r'res.norm')

print('chi quadro = %.3f (%d dof)' % (chisq, ndof))
print('parametri fit: g/nu,cost')
print(pars)
print('errori fit')
print(err)


## fit bilog per determinare g  suscettività in funzione di t


def fit_log(t,g,c):
    return c*t**(-g)
t=np.array([])
chi=np.array([])
dchi=np.array([])
for i in range (0,50):
    if (1/beta[i] - 1/beta_c)>0.05:
        t=np.append(t,1/beta[i] - 1/beta_c)
        chi=np.append(chi,X[i,9])
        dchi=np.append(dchi,dX[i,9])

init=(1.75,2)
pars, covm = curve_fit(fit_log, t, chi,init, sigma=dchi,absolute_sigma=True)
err = np.sqrt(covm.diagonal())
chisq = np.sum(((chi - fit_log(t, *pars))/dchi)**2)
ndof = len(t) - len(pars)

plt.figure(9)
xx=np.linspace(min(t),max(t),10000)
plt.minorticks_on()
plt.title(r'stima di $\gamma$ dal fit ($\chi$ al variare di t)')
plt.errorbar(t,chi,dchi,marker='.',capsize=4,linestyle='')
plt.plot(xx,fit_log(xx,*pars))
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r't [a.u.]')
plt.ylabel(r'$\chi$ [a.u.]')

print('chi quadro = %.3f (%d dof)' % (chisq, ndof))
print('parametri fit: g,cost')
print(pars)
print('errori fit')
print(err)
'''

## fit bilog per determinare b magentizzazione in funzione di t
'''
def fit_log(t,b,c):
    return c*t**(-g)


t=np.array([])
M1=np.array([])
dM1=np.array([])
for i in range (0,50):
    if (1/beta[i] - 1/beta_c)>0.1:
        t=np.append(t,1/beta[i] - 1/beta_c)
        M1=np.append(M1,M[i,9])
        dM1=np.append(dM1,dM[i,9])
init=(1/8,2)
pars, covm = curve_fit(fit_log, t, M1, init,sigma=dM1,absolute_sigma=True)
err = np.sqrt(covm.diagonal())
chisq = np.sum(((M1 - fit_log(t, *pars))/dM1)**2)
ndof = len(t) - len(pars)

plt.figure(10)
xx=np.linspace(min(t),max(t),10000)
plt.minorticks_on()
plt.title(r'stima di $\beta$ dal fit ($\|M|$ al variare di t)')
plt.errorbar(t,M1,dM1,marker='.',capsize=4,linestyle='')
plt.plot(xx,fit_log(xx,*pars))
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r't [a.u.]')
plt.ylabel(r'$|M|$ [a.u.]')

print('chi quadro = %.3f (%d dof)' % (chisq, ndof))
print('parametri fit: g/nu,cost')
print(pars)
print('errori fit')
print(err)

## fit bilog per determinare alfa calore specifico in funzione di t
def fit_log(t,a,c):
    return c*t**(-a)

t=np.array([])
C1=np.array([])
dC1=np.array([])
for i in range (0,50):
    if (1/beta[i] - 1/beta_c)>0.1:
        t=np.append(t,1/beta[i] - 1/beta_c)
        C1=np.append(C1,C[i,9])
        dC1=np.append(dC1,dC[i,9])

pars, covm = curve_fit(fit_log, t, C1, sigma=dC1,absolute_sigma=True)
err = np.sqrt(covm.diagonal())
chisq = np.sum(((C1 - fit_log(t, *pars))/dC1)**2)
ndof = len(t) - len(pars)

plt.figure(11)
xx=np.linspace(min(t),max(t),10000)
plt.minorticks_on()
plt.title(r'stima di $\beta$ dal fit (Calore specifico al variare di t)')
plt.errorbar(t,C1,dC1,marker='.',capsize=4,linestyle='')
plt.plot(xx,fit_log(xx,*pars))
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r't [a.u.]')
plt.ylabel(r'$C_v$ [a.u.]')

print('chi quadro = %.3f (%d dof)' % (chisq, ndof))
print('parametri fit: g/nu,cost')
print(pars)
print('errori fit')
print(err)
'''

plt.show()
