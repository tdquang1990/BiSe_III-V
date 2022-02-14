from cmath import sqrt
import numpy as np

def EMquantities(materials,omega):
    if materials == 'Air':
        eps=1
        eps = [eps, eps, eps]
        mu=1
    elif materials =='BiSe':
        eps_inf = 1.0
        omD_plasma_cm = 908.66  #drude plasma
        gamD_cm = 7.43          #Drude damping
        omLa_O_cm = 63.03       #Lorentz alpha phonon
        omLa_plasma_cm = 675.9  #Lorentz alpha strength
        gamLa_cm = 17.5         #Lorentz alpha damping
        omLb_O_cm = 126.94      #Lorentz beta phonon
        omLb_plasma_cm = 100    #Lorentz beta strength
        gamLb_cm = 10           #Lorentz beta damping
        omLo_O_cm = 2029.5      #Lorentz band gap (omega)
        omLo_plasma_cm = 11249  #Lorentz band gap (omega) strength
        gamLo_cm = 3920.5       #Lorentz band gap damping
        epsD = -omD_plasma_cm**2/(omega**2 + 1j*omega*gamD_cm)
        epsLa = omLa_plasma_cm**2/(omLa_O_cm**2 - omega**2 - 1j*omega*gamLa_cm)
        epsLb = omLb_plasma_cm**2/(omLb_O_cm**2 - omega**2 - 1j*omega*gamLb_cm)
        epsLo = omLo_plasma_cm**2/(omLo_O_cm**2 - omega**2 - 1j*omega*gamLo_cm)
        eps = eps_inf + 0*epsD + epsLa + epsLb + 0*epsLo     #The Drude term contribution sometimes changes based on material quality
        eps = [eps, eps, eps]
        mu =  1
    elif materials == 'BIS':
        #eps_inf_bis = 7.030       #no drude term for s polarization
        eps_inf_bis = 10.26       #no drude term for P polarization
        omLa_O_cm_bis = 88.84     #Lorentz alpha phonon
        omLa_plasma_cm_bis = 361.936  #Lorentz alpha strength
        gamLa_cm_bis = 60.247  #Lorentz alpha damping
        omLb_O_cm_bis = 142.30  #Lorentz beta phonon
        omLb_plasma_cm_bis = 186.111  #Lorentz beta strength
        gamLb_cm_bis = 22.29  #Lorentz beta damping
        epsLa_bis = omLa_plasma_cm_bis**2/(omLa_O_cm_bis**2 - omega**2 - 1j*omega*gamLa_cm_bis)
        epsLb_bis = omLb_plasma_cm_bis**2/(omLb_O_cm_bis**2 - omega**2 - 1j*omega*gamLb_cm_bis)
        eps = eps_inf_bis + epsLa_bis + epsLb_bis
        eps = [eps, eps, eps]
        mu = 1
    elif materials == 'GaAs':
        eps = 12.9
        eps = [eps, eps, eps]
        mu = 1
    elif materials == 'AlGaAs':
        eps = 12.9
        eps = [eps, eps, eps]
        mu = 1
    elif materials == 'AlO':
        eps = 3.2**2+(3.2**2 -1)*(20.4*10**-4 *omega)**2 +1j*0.036*(3.2**2 -1)*(20.4*10**-4 *omega) #Sapphire substrate
        eps = [eps, eps, eps]
        mu = 1        
    elif materials == 'GaQw':
        from Physconstant import c0,q0,m0,eps0,hb
        omega = 2*np.pi*100*c0*omega 
        m = 0.0665*m0
        g1 = 1e12*2*np.pi  
        g2 = 1e12*2*np.pi  
        f12 = 1   #Dipole matrix element for ISBT
        ns = 1e12   #Doping in QW
        n12 = 1*ns 
        Lqw = 125*1e-9 
        Lmqw = 125*1e-9 
        E12 = 0.0241*q0 
        epsb = 9.88
        epsw = 10.36
        epsy = (1-Lqw/Lmqw)*epsb + epsw*Lqw/Lmqw
        epsz = 1/((1-Lqw/Lmqw)/epsb+Lqw/(epsw*Lmqw)) 
        omega_p = sqrt(n12*1e4*q0**2/(m*eps0*Lmqw)) #
        epsxx = epsy -  (ns*1e4*q0**2/(m*eps0*Lmqw))/(omega**2+1j*omega*g1)
        epsyy = epsy -  (ns*1e4*q0**2/(m*eps0*Lmqw))/(omega**2+1j*omega*g1)
        epszz = 1/(1/epsz-(omega_p**2*f12/(2*omega*g2*epsw))/((E12**2-hb**2*omega**2)/(2*hb**2*g2*omega)-1j))
        eps = epszz
        eps = [epsxx, epsyy, epszz]
        mu = 1
    return(eps,mu)

def EMquantity(materials,omega,gamma,Ns):
    if materials == 'Air':
        eps=1
        eps = [eps, eps, eps]
        mu=1
    elif materials =='BiSe':
        eps_inf = 1.0
        omD_plasma_cm = 908.66  #drude plasma
        gamD_cm = 7.43          #Drude damping
        omLa_O_cm = 63.03       #Lorentz alpha phonon
        omLa_plasma_cm = 675.9  #Lorentz alpha strength
        gamLa_cm = 17.5         #Lorentz alpha damping
        omLb_O_cm = 126.94      #Lorentz beta phonon
        omLb_plasma_cm = 100    #Lorentz beta strength
        gamLb_cm = 10           #Lorentz beta damping
        omLo_O_cm = 2029.5      #Lorentz band gap (omega)
        omLo_plasma_cm = 11249  #Lorentz band gap (omega) strength
        gamLo_cm = 3920.5       #Lorentz band gap damping
        epsD =  -omD_plasma_cm**2/(omega**2 + 1j*omega*gamD_cm)
        epsLa = omLa_plasma_cm**2/(omLa_O_cm**2 - omega**2 - 1j*omega*gamLa_cm)
        epsLb = omLb_plasma_cm**2/(omLb_O_cm**2 - omega**2 - 1j*omega*gamLb_cm)
        epsLo = omLo_plasma_cm**2/(omLo_O_cm**2 - omega**2 - 1j*omega*gamLo_cm)
        eps = eps_inf + 0*epsD + epsLa + epsLb + 0*epsLo     #The Drude term contribution sometimes changes based on material quality
        eps = [eps, eps, eps]
        mu =  1
    elif materials == 'BIS':
        #eps_inf_bis = 7.030       #no drude term for s polarization
        eps_inf_bis = 10.26       #no drude term for P polarization
        omLa_O_cm_bis = 88.84     #Lorentz alpha phonon
        omLa_plasma_cm_bis = 361.936  #Lorentz alpha strength
        gamLa_cm_bis = 60.247  #Lorentz alpha damping
        omLb_O_cm_bis = 142.30  #Lorentz beta phonon
        omLb_plasma_cm_bis = 186.111  #Lorentz beta strength
        gamLb_cm_bis = 22.29  #Lorentz beta damping
        epsLa_bis = omLa_plasma_cm_bis**2/(omLa_O_cm_bis**2 - omega**2 - 1j*omega*gamLa_cm_bis)
        epsLb_bis = omLb_plasma_cm_bis**2/(omLb_O_cm_bis**2 - omega**2 - 1j*omega*gamLb_cm_bis)
        eps = eps_inf_bis + epsLa_bis + epsLb_bis
        eps = [eps, eps, eps]
        mu = 1
    elif materials == 'GaAs':
        eps = 12.9
        eps = [eps, eps, eps]
        mu = 1
    elif materials == 'AlGaAs':
        eps = 12.9
        eps = [eps, eps, eps]
        mu = 1
    elif materials == 'AlO':
        eps = 3.2**2+(3.2**2 -1)*(20.4*10**-4 *omega)**2 +1j*0.036*(3.2**2 -1)*(20.4*10**-4 *omega) #Sapphire substrate
        eps = [eps, eps, eps]
        mu = 1        
    elif materials == 'GaQw':
        from Physconstant import c0,q0,m0,eps0,hb
        omega = 2*np.pi*100*c0*omega 
        m = 0.0665*m0
        g1 = gamma*2e12*np.pi  
        g2 = gamma*2e12*np.pi 
        f12 = 1   #Dipole matrix element for ISBT
        ns = Ns   #Doping in QW
        n12 = 1*ns     
        Lqw = 120*1e-9 
        Lmqw = 120*1e-9 
        E12 = 0.0241*q0 
        epsb = 9.88
        epsw = 10.36
        epsy = (1-Lqw/Lmqw)*epsb + epsw*Lqw/Lmqw
        epsz = 1/((1-Lqw/Lmqw)/epsb+Lqw/(epsw*Lmqw))    
        omega_p = sqrt(n12*1e4*q0**2/(m*eps0*Lmqw)) 
        epsxx = epsy -  (ns*1e4*q0**2/(m*eps0*Lmqw))/(omega**2+1j*omega*g1)
        epsyy = epsy -  (ns*1e4*q0**2/(m*eps0*Lmqw))/(omega**2+1j*omega*g1)
        epszz = 1/(1/epsz-(omega_p**2*f12/(2*omega*g2*epsw))/((E12**2-hb**2*omega**2)/(2*hb**2*g2*omega)-1j))
        eps = epszz
        eps = [epsxx, epsyy, epszz]
        mu = 1
    return(eps,mu)
    
def EMquantity1(materials,omega,gamma,Ns):
    if materials == 'Air':
        eps=1
        eps = [eps, eps, eps]
        mu=1
    elif materials =='BiSe':
        eps_inf = 1.0
        omD_plasma_cm = 908.66  #drude plasma
        gamD_cm = 7.43          #Drude damping
        omLa_O_cm = 63.03       #Lorentz alpha phonon
        omLa_plasma_cm = 675.9  #Lorentz alpha strength
        gamLa_cm =  gamma/0.03#17.5         #Lorentz alpha damping
        omLb_O_cm = 126.94      #Lorentz beta phonon
        omLb_plasma_cm = 100    #Lorentz beta strength
        gamLb_cm = gamma/0.03 #10           #Lorentz beta damping
        omLo_O_cm = 2029.5      #Lorentz band gap (omega)
        omLo_plasma_cm = 11249  #Lorentz band gap (omega) strength
        gamLo_cm = 3920.5       #Lorentz band gap damping
        epsD =  -omD_plasma_cm**2/(omega**2 + 1j*omega*gamD_cm)
        epsLa = omLa_plasma_cm**2/(omLa_O_cm**2 - omega**2 - 1j*omega*gamLa_cm)
        epsLb = omLb_plasma_cm**2/(omLb_O_cm**2 - omega**2 - 1j*omega*gamLb_cm)
        epsLo = omLo_plasma_cm**2/(omLo_O_cm**2 - omega**2 - 1j*omega*gamLo_cm)
        eps = eps_inf + 0*epsD + epsLa + epsLb + 0*epsLo     #The Drude term contribution sometimes changes based on material quality
        eps = [eps, eps, eps]
        mu =  1
    elif materials == 'BIS':
        #eps_inf_bis = 7.030       #no drude term for s polarization
        eps_inf_bis = 10.26       #no drude term for P polarization
        omLa_O_cm_bis = 88.84     #Lorentz alpha phonon
        omLa_plasma_cm_bis = 361.936  #Lorentz alpha strength
        gamLa_cm_bis = 60.247  #Lorentz alpha damping
        omLb_O_cm_bis = 142.30  #Lorentz beta phonon
        omLb_plasma_cm_bis = 186.111  #Lorentz beta strength
        gamLb_cm_bis = 22.29  #Lorentz beta damping
        epsLa_bis = omLa_plasma_cm_bis**2/(omLa_O_cm_bis**2 - omega**2 - 1j*omega*gamLa_cm_bis)
        epsLb_bis = omLb_plasma_cm_bis**2/(omLb_O_cm_bis**2 - omega**2 - 1j*omega*gamLb_cm_bis)
        eps = eps_inf_bis + epsLa_bis + epsLb_bis
        eps = [eps, eps, eps]
        mu = 1
    elif materials == 'GaAs':
        eps = 12.9
        eps = [eps, eps, eps]
        mu = 1
    elif materials == 'AlGaAs':
        eps = 12.9
        eps = [eps, eps, eps]
        mu = 1
    elif materials == 'AlO':
        eps = 3.2**2+(3.2**2 -1)*(20.4*10**-4 *omega)**2 +1j*0.036*(3.2**2 -1)*(20.4*10**-4 *omega) #Sapphire substrate
        eps = [eps, eps, eps]
        mu = 1        
    elif materials == 'GaQw':
        from Physconstant import c0,q0,m0,eps0,hb
        omega = 2*np.pi*100*c0*omega 
        m = 0.0665*m0
        g1 = 2e12*np.pi  
        g2 = 2e12*np.pi 
        f12 = 1   #Dipole matrix element for ISBT
        ns = Ns   #Doping in QW
        n12 = 1*ns     
        Lqw = 120*1e-9 
        Lmqw = 120*1e-9 
        E12 = 0.0241*q0 
        epsb = 9.88
        epsw = 10.36
        epsy = (1-Lqw/Lmqw)*epsb + epsw*Lqw/Lmqw
        epsz = 1/((1-Lqw/Lmqw)/epsb+Lqw/(epsw*Lmqw))    
        omega_p = sqrt(n12*1e4*q0**2/(m*eps0*Lmqw)) 
        epsxx = epsy -  (ns*1e4*q0**2/(m*eps0*Lmqw))/(omega**2+1j*omega*g1)
        epsyy = epsy -  (ns*1e4*q0**2/(m*eps0*Lmqw))/(omega**2+1j*omega*g1)
        epszz = 1/(1/epsz-(omega_p**2*f12/(2*omega*g2*epsw))/((E12**2-hb**2*omega**2)/(2*hb**2*g2*omega)-1j))
        eps = epszz
        eps = [epsxx, epsyy, epszz]
        mu = 1
    return(eps,mu)
