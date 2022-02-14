import numpy as np
from scipy.linalg import block_diag
from Smatrix import Smatrix

def Gmatrix(d,eps,sigma,N,omega,kpa,mode):
    P0  = {}
    P01 = {}
    S   = {}
    for i in range (0,N-1):
        if (i == 0):
            (Sc,_,_) =Smatrix(d[i],eps[i],eps[i+1],sigma[i],omega,kpa,mode)
        else:
            (Sc,_,P) =Smatrix(d[i],eps[i],eps[i+1],sigma[i],omega,kpa,mode)
            P0[i-1]  = -block_diag(P, 0)
            P01[i-1] = -block_diag(0,P)
        S[i]=np.linalg.inv(Sc)
    if (N > 2):
        D0 = S[0]
        D1 = P0[0]
        D2 = P01[0]
        for i in range(1,N-1):
            D0 = block_diag(D0,S[i])
        for i in range(1,N-2):
            D1 = block_diag(D1,P0[i])
            D2 = block_diag(D2,P01[i])   
        D1  = np.hstack((np.vstack((np.zeros((2,2*(N-2))),D1)),np.zeros((2*(N-1),2))))
        D2  = np.hstack((np.zeros((2*(N-1),2)),np.vstack((D2,np.zeros((2,2*(N-2)))))))
        G = np.linalg.inv(D0+D1+D2)
    else:
        G = S[0]
    t = G[2*(N-2),0]
    r = G[1,0]
    sT = np.multiply(t,np.conj(t))
    sR = np.multiply(r,np.conj(r))
    return(t,r,sT,sR,G)