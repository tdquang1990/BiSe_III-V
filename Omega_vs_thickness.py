import numpy as np
import matplotlib.pyplot as plt
from Gmatrix import Gmatrix
from structures import *

mode = 'TM'
kpa = 5e6
Nqw =1
n=1000
dmin =0
dmax = 200
omegamin = 5  
omegamax = 300    
dt = np.linspace(dmin,dmax,n)
omega = np.linspace(omegamin,omegamax,n)
R = np.zeros((n,n),dtype=complex)

for i in range(0,n):
    print(i)
    (eps,_,sigma,N) = structuremqw_thickness(omega[i],Nqw)
    for j in range(0,n):
        dqw = np.array([100, 25]*Nqw)
        dti = np.array([0,dt[j],10])
        d = np.concatenate((dti,dqw,np.array([100,0])),axis=0)*1e-9
        (_,R[i,j],_,_,_) = Gmatrix(d,eps,sigma,N,omega[i],kpa,mode)
        
plt.pcolor(dt,0.03*omega,abs(np.imag(R)), cmap='hot',vmin=0,vmax=0.5)
plt.xlabel('$d_{sp} (nm)$',fontsize=14,fontname="Times New Roman")
plt.ylabel('$\omega (THz)$',fontsize=14,fontname="Times New Roman")
cbar=plt.colorbar()
cbar.ax.set_title('Im(r)')
