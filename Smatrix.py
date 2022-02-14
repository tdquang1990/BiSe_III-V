from cmath import sqrt
import numpy as np
from Physconstant import*

def Smatrix(dl,epsl,epsr,sigma,omega,kpa,mode):
    omega = omega*(2*np.pi*c0*100) 
    IM = np.zeros((2,2),dtype=complex)
    S_pr = np.eye(2,dtype=complex)
    S = np.zeros((2,2),dtype=complex)
    kzl = sqrt(epsl[0]*omega**2/c0**2 - (epsl[0]/epsl[2])*kpa**2)
    kzr = sqrt(epsr[0]*omega**2/c0**2 - (epsr[0]/epsr[2])*kpa**2)
    P0 = np.exp(-1j*kzl*dl)
    P0t = np.exp(1j*kzl*dl)
    
    if mode =='TM':
        IM[0,0] = 1 + epsr[0]*kzl/(epsl[0]*kzr) + kzl*sigma/(eps0*epsl[0]*omega)
        IM[0,1] = 1 - epsr[0]*kzl/(epsl[0]*kzr) + kzl*sigma/(eps0*epsl[0]*omega)
        IM[1,0] = 1 - epsr[0]*kzl/(epsl[0]*kzr) - kzl*sigma/(eps0*epsl[0]*omega)
        IM[1,1] = 1 + epsr[0]*kzl/(epsl[0]*kzr) - kzl*sigma/(eps0*epsl[0]*omega)
    elif mode == 'TE':
        IM[0,0] = 1 + kzr/kzl + mu0*sigma*omega/kzl
        IM[1,0] = 1 - kzr/kzl + mu0*sigma*omega/kzl
        IM[0,1] = 1 - kzr/kzl - mu0*sigma*omega/kzl
        IM[1,1] = 1 + kzr/kzl - mu0*sigma*omega/kzl
        
    IM = 0.5*IM
    S11 = S_pr[0,0]
    S12 = S_pr[0,1]
    S21 = S_pr[1,0]
    S22 = S_pr[1,1]
    I11 = IM[0,0]
    I12 = IM[0,1]
    I21 = IM[1,0]
    I22 = IM[1,1]
    temp1 = 1/(1 - (S12/I11)*I21)
    temp2 = S11/I11
    S11_new = temp1*temp2
    S12_new=temp1*((S12*I22-I12)/I11)
    S21_new=S22*I21*S11_new+S21
    S22_new=S22*I21*S12_new+S22*I22
    S[0,0] = S11_new
    S[0,1] = S12_new
    S[1,0] = S21_new
    S[1,1] = S22_new
    return(S,P0,P0t)