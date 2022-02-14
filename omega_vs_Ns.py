import numpy as np
import matplotlib.pyplot as plt
from Gmatrix import Gmatrix
from structures import *

mode = 'TM'
kpa = 5e6
Nqw =1
gamma = 1
n=1000
Nsmin =1e11
Nsmax = 5e12
omegamin = 5  
omegamax = 300   
Ns= np.linspace(Nsmin,Nsmax,n)
omega = np.linspace(omegamin,omegamax,n)
R = np.zeros((n,n),dtype=complex)

for i in range(0,n):
    print(i)
    for j in range(0,n):
        (d,eps,_,sigma,N) = structure_materials(omega[i],gamma,Ns[j],Nqw)
        (_,R[i,j],_,_,_) = Gmatrix(d,eps,sigma,N,omega[i],kpa,mode)
        
plt.pcolor(Ns,0.03*omega,abs(np.imag(R)), cmap='hot',vmin=0,vmax=0.5)
plt.xlabel('$Ns$',fontsize=14,fontname="Times New Roman")
plt.ylabel('$\omega (THz)$',fontsize=14,fontname="Times New Roman")
cbar=plt.colorbar()
cbar.ax.set_title('Im(r)')
