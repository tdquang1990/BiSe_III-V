import numpy as np 
import matplotlib.pyplot as plt

dTI=100e-9
c0 = 3e8 #speed of light in vacuum
q0 = 1.602176634e-19 #electron charge
hb = 1.054571817e-34 # reduced Planks constant
h = 2*np.pi*hb # Planks constant
eps0 = 8.8541878128e-12 # free sapce permitivity
mu0 =  4*np.pi*1e-7 # Magnetic free space permeability 
vf = 5*1e5  # Fermi velocity of surface states in Bi2Se3
mu = 0.17*q0 # chemical potential correspondent to 1E+13 cm-2 carrier density
kf = mu/hb/vf #Fermi wavevector
A = eps0*h/vf/kf/q0**2 #constant
# % % % % % % Following Drude-Lorentz Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eps_inf = 1.0
omD_plasma_cm = 908.66 # drude plasma
gamD_cm = 7.43 # Drude damping
omLa_O_cm = 63.03 # Lorentz alpha phonon
omLa_plasma_cm = 675.9  # Lorentz alpha strength
gamLa_cm = 17.5       # Lorentz alpha damping
omLb_O_cm = 126.94    # Lorentz beta phonon
omLb_plasma_cm = 100 # Lorentz beta strength
gamLb_cm = 10 # Lorentz beta damping
omLo_O_cm = 2029.5 # Lorentz band gap (omega)
omLo_plasma_cm = 11249 # Lorentz band gap (omega) strength
gamLo_cm = 3920.5   # Lorentz band gap damping
eps_a = 1.0
n = 500 
omegamin =5
omegamax =300
omega = np.linspace (omegamin,omegamax, n);
k_CLAS_p= np.zeros(n,dtype=complex)

for i in range (0,n):
    eps_sapp = 3.2**2+(3.2**2 -1)*(20.4*10**-4 *omega[i])**2 +1j*0.036*(3.2**2 -1)*(20.4*10**-4 *omega[i]) # Sapphire substrate
    epsD =  -omD_plasma_cm**2/(omega[i]**2 + 1j*omega[i]*gamD_cm)
    epsLa = omLa_plasma_cm**2/(omLa_O_cm**2 - omega[i]**2 - 1j*omega[i]*gamLa_cm)
    epsLb = omLb_plasma_cm**2/(omLb_O_cm**2 - omega[i]**2 - 1j*omega[i]*gamLb_cm)
    epsLo = omLo_plasma_cm**2/(omLo_O_cm**2 - omega[i]**2 - 1j*omega[i]*gamLo_cm)
    eps_bs = eps_inf + 0*epsD + epsLa + epsLb; + 0*epsLo #The Drude term contribution sometimes changes based on material quality
    k_CLAS_p[i]= A*(omega[i]*(2*np.pi*c0*100))**2*(eps_a + eps_sapp)/(1 - A*(omega[i]*(2*np.pi*c0*100))**2*dTI*eps_bs)

plt.plot(k_CLAS_p/1e7,0.03*omega,'c')
plt.xlim([1e-2, 1])
plt.show()