import numpy as np
import matplotlib.pyplot as plt
from EMquantities import *

materials = 'BiSe'
n=10000
omegamin = 5  
omegamax = 250    
omega = np.linspace(omegamin,omegamax,n)
epsxx = np.zeros(n,dtype=complex)
epszz = np.zeros(n,dtype=complex)

for i in range(0,n):
    (eps,_) = EMquantities(materials,omega[i])
    epsxx [i] = eps[0]
    epszz[i] = eps[2]
    
plt.subplot(212)
#plt.plot(0.03*omega,np.real(epsxx),'g')
plt.plot(0.03*omega,np.imag(epsxx),'g')
plt.xlabel('$\omega(THz)$',fontsize=14,fontname="Times New Roman")
plt.ylabel('$\epsilon_{Bi_{2}Se_{3}} $ (Im)',fontsize=14,fontname="Times New Roman")
plt.legend(['alpha', 'alpha', 'beta'])