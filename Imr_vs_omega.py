import numpy as np
import matplotlib.pyplot as plt
from Gmatrix import Gmatrix
from structures import *

mode = 'TM'
n=1000
kpa =4e6
omegamin = 1  
omegamax = 300    
Nqw =1
omega = np.linspace(omegamin,omegamax,n)
R = np.zeros(n,dtype=complex)
for i in range(0,n):
    print(i)
    (d,eps,_,sigma,N) = structuremqw(omega[i],Nqw)
    (_,R[i],_,_,_) = Gmatrix(d,eps,sigma,N,omega[i],kpa,mode)
    
        
plt.semilogy(0.03*omega,abs(np.imag(R)),'b')
plt.xlabel('$\omega (THz)$',fontsize=14,fontname="Times New Roman")
plt.ylabel('$Im(r)$ (log.)',fontsize=14,fontname="Times New Roman")
plt.xlim((0, 9))


