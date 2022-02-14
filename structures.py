import numpy as np
from EMquantities import *
from Physconstant import*
def structure(omega):
    Ef = 0.17*q0
    tau = 1e-12
    N=3
    d = np.array([0,100,0])*1e-9
    materials  = np.array(['Air','BiSe','AlO'])
    mu = np.zeros(N,dtype=complex)
    eps = {}
    for i in range(0,N):
        (eps[i],mu[i]) = EMquantities(materials[i],omega)    
    sig =1j*Ef*q0**2/(4*np.pi*hb**2*(omega*2*np.pi*c0*100+1j*2*np.pi/tau))
    # Extra 2DEG at Bi2Se3/GaAs interface
    N2DE = 1e13   
    G2DE = 1e12
    m2DE = m0
    sig_extra = N2DE*1e4*q0**2/(m2DE*(G2DE-1j*omega*2*np.pi*c0*100)) 
    sigma = np.array([sig, sig])
    return(d,eps,mu,sigma,N)

def structuremeta(omega):
    N=4
    d = np.array([0,100,25,0])*1e-9
    materials  = np.array(['Air','AlGaAs','GaQw','GaAs'])
    mu = np.zeros(N,dtype=complex)
    eps = {}
    for i in range(0,N):
        (eps[i],mu[i]) = EMquantities(materials[i],omega)       
    sigma = np.array([0, 0, 0])
    return(d,eps,mu,sigma,N)
    
def structuregraphene(omega):
    mu_c=q0
    Gamma = 1/(0.5*1e-12)
    T = 300
    N=3
    d = np.array([0,50,0])*1e-9
    materials  = np.array(['Air','GaAs','AlO'])
    mu = np.zeros(N,dtype=complex)
    eps = {}
    for i in range(0,N):
        (eps[i],mu[i]) = EMquantities(materials[i],omega)    
    sig_inter = 1j*q0**2/(4*np.pi*hb)*np.log((2*np.absolute(mu_c)-hb*(omega*2*np.pi*c0*100+1j*Gamma)) /(2*np.absolute(mu_c)+hb*(omega*2*np.pi*c0*100+1j*Gamma)))
    sig_intra = 1j*q0**2*kb*T/(np.pi*hb**2*(omega*2*np.pi*c0*100+1j*Gamma))*(mu_c/(kb*T)+2*np.log(np.exp(-mu_c/(kb*T))+1))
    sig_scatteing = 1j*q0**2*kb*T/(np.pi*hb**2*omega*2*np.pi*c0*100)*(mu_c/(kb*T)+2*np.log(np.exp(-mu_c/(kb*T))+1))
    sigma = np.array([sig_intra+sig_inter, 0])
    return(d,eps,mu,sigma,N)
    
def structure_thickness(omega):
    Ef = 0.17*q0
    tau = 1e-12
    N=3
    materials  = np.array(['Air', 'BiSe','AlO'])
    mu = np.zeros(N,dtype=complex)
    eps = {}
    for i in range(0,N):
        (eps[i],mu[i]) = EMquantities(materials[i],omega)  
    sig =1j*Ef*q0**2/(4*np.pi*hb**2*(omega*2*np.pi*c0*100+1j*2*np.pi/tau))
    # Extra 2DEG at BiSe/GaAs interface
    N2DE = 1e13
    G2DE = 1e12
    m2DE = m0
    sig_extra = N2DE*1e4*q0**2/(m2DE*(G2DE-1j*omega*2*np.pi*c0*100))
    sigma = np.array([sig, sig])
    return(eps,mu,sigma,N)

def structuremqw(omega,Nqw):
    N=2*Nqw+5
    Ef = 0.17*q0
    tau = 1e-12
    mu = np.zeros(N,dtype=complex)      
    sig =1j*Ef*q0**2/(4*np.pi*hb**2*(omega*2*np.pi*c0*100+1j*2*np.pi/tau))
    # Extra 2DEG at BiSe/GaAs interface
    N2DE = 2.6e14
    G2DE = 1e12
    m2DE = m0
    sig_extra = N2DE*1e4*q0**2/(m2DE*(G2DE-1j*omega*2*np.pi*c0*100))
    # Quantum well's parameters
    dqw = np.array([100, 25]*Nqw)
    mqw = np.array(['AlGaAs','GaQw']*Nqw)
    sigqw = np.array([0]*2*Nqw)
    # TI's parameters
    dti = np.array([0,100,10])
    TI  = np.array(['Air', 'BiSe','GaAs'])
    sigti = np.array([sig , sig])
    d = np.concatenate((dti,dqw,np.array([100,0])),axis=0)*1e-9
    materials  = np.concatenate((TI,mqw,np.array(['AlGaAs','GaAs'])),axis=0)
    sigma = np.concatenate((sigti,sigqw,np.array([0,0])),axis=0)
    eps = {} 
    for i in range(0,N):
        (eps[i],mu[i]) = EMquantities(materials[i],omega)
    return(d,eps,mu,sigma,N)

def structuremqw_thickness(omega,Nqw):
    Ef = 0.17*q0
    tau = 1e-12
    N=2*Nqw+5
    # Quantum well's parameters
    mqw = np.array(['AlGaAs','GaQw']*Nqw)
    sigqw = np.array([0]*2*Nqw)
    #TI's parameters
    TI  = np.array(['Air', 'BiSe','GaAs'])
    sig =1j*Ef*q0**2/(4*np.pi*hb**2*(omega*2*np.pi*c0*100+1j*2*np.pi/tau))
    # Extra 2DEG at BiSe/GaAs interface
    N2DE = 1e13
    G2DE = 1e12
    m2DE = m0
    sig_extra = N2DE*1e4*q0**2/(m2DE*(G2DE-1j*omega*2*np.pi*c0*100))
    sigti = np.array([sig , sig])
    materials  = np.concatenate((TI,mqw,np.array(['AlGaAs','GaAs'])),axis=0)
    mu = np.zeros(N,dtype=complex)
    eps = {}
    for i in range(0,N):
        (eps[i],mu[i]) = EMquantities(materials[i],omega)
    sigma = np.concatenate((sigti,sigqw,np.array([0,0])),axis=0)
    return(eps,mu,sigma,N)
def structure_materials(omega,gamma,Ns,Nqw):
    N=2*Nqw+5
    Ef = 0.17*q0
    tau = 1e-12
    mu = np.zeros(N,dtype=complex)      
    sig =1j*Ef*q0**2/(4*np.pi*hb**2*(omega*2*np.pi*c0*100+1j*2*np.pi/tau))
    # Extra 2DEG at BiSe/GaAs interface
    N2DE = 1e13
    G2DE = 1e12
    m2DE = m0
    sig_extra = N2DE*1e4*q0**2/(m2DE*(G2DE-1j*omega*2*np.pi*c0*100))
    # Quantum well's parameters
    dqw = np.array([100, 25]*Nqw)
    mqw = np.array(['AlGaAs','GaQw']*Nqw)
    sigqw = np.array([0]*2*Nqw)
    # TI's parameters
    dti = np.array([0,100,10])
    TI  = np.array(['Air', 'BiSe','GaAs'])
    sigti = np.array([sig , sig])
    d = np.concatenate((dti,dqw,np.array([100,0])),axis=0)*1e-9
    materials  = np.concatenate((TI,mqw,np.array(['AlGaAs','GaAs'])),axis=0)
    sigma = np.concatenate((sigti,sigqw,np.array([0,0])),axis=0)
    eps = {} 
    for i in range(0,N):
        (eps[i],mu[i]) = EMquantity(materials[i],omega,gamma,Ns)     
    return(d,eps,mu,sigma,N)