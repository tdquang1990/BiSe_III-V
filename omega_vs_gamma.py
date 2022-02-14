import numpy as np
import matplotlib.pyplot as plt
from Gmatrix import Gmatrix
from structures import *

mode = 'TM'
kpa = 3e6
Nqw =1
Ns = 1.5e12 # Doping per QW
n=1000
gmin =0.1
gmax = 2
omegamin = 5  
omegamax = 300   
gamma = np.linspace(gmin,gmax,n)
omega = np.linspace(omegamin,omegamax,n)
R = np.zeros((n,n),dtype=complex)

for i in range(0,n):
    print(i)
    for j in range(0,n):
        (d,eps,_,sigma,N) = structure_materials(omega[i],gamma[j],Ns,Nqw)
        (_,R[i,j],_,_,_) = Gmatrix(d,eps,sigma,N,omega[i],kpa,mode)
        
plt.pcolor(gamma,0.03*omega,abs(np.imag(R)), cmap='hot',vmin=0,vmax=0.3)
plt.xlabel('$\gamma (THz)$',fontsize=14,fontname="Times New Roman")
plt.ylabel('$\omega (THz)$',fontsize=14,fontname="Times New Roman")
cbar=plt.colorbar()
cbar.ax.set_title('Im(r)')
