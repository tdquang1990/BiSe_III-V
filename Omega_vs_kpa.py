import numpy as np
import matplotlib.pyplot as plt
from Gmatrix import Gmatrix
from structures import *

mode = 'TM'
n=1000
kmin =1e5
kmax = 1e7
omegamin = 5   
omegamax = 300
Nqw =1
kpa = np.linspace(kmin,kmax,n)
omega = np.linspace(omegamin,omegamax,n)
R = np.zeros((n,n),dtype=complex)

for i in range(0,n):
    print(i)
    (d,eps,_,sigma,N) = structuremqw(omega[i],Nqw)
    for j in range(0,n):
        (_,R[i,j],_,_,_) = Gmatrix(d,eps,sigma,N,omega[i],kpa[j],mode)
        
plt.pcolor(kpa/1e7,0.03*omega,abs(np.imag(R)), cmap='hot',vmin=0,vmax=0.5)
plt.xlabel('$k_x (10^{5} cm^{-1})$',fontsize=14,fontname="Times New Roman")
plt.ylabel('$\omega (THz)$',fontsize=14,fontname="Times New Roman")
cbar=plt.colorbar()
cbar.ax.set_title('Im(r)')
