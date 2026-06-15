# Calculate the Derivative Rx

# Imports -------------------------------------------------------------------------------
import numpy as np
from NON_DIM_HB_AE.NON_DIM_HB_FNL import HB_FNL

def HB_Rx(X, ntot, n, nfft, G, G_inv, qn, FNL, Z, U):
    #Description of inputs-----------------------------------------------------------------
    # X ==> harmonic coefficients of the system
    # ntot ==> total dimension of the system ntot = n*(2*n_hbm+1)
    # n ==> number of DOF
    # nfft ==> number of FFT time points to determine the harmonic coefficients of Fnl
    # G ==> DFT matrix Gamma
    # G_inv ==> inverse DFT matrix Gamma-1
    # knl ==> non-linear stiffness in matrix form [knl1, knl2, ... , knln]
    # qn ==> non-linear force input distribution matrix
    # FNL ==> Fourrier coefficients for nonlinear force
    # Z ==> Linear Forces Z = w NABLA kron In - IH kron Q
    
    epsx = 1e-6; 
    Xp = X.copy()
    dFNLdX = np.zeros((ntot,ntot))

    for ii in range(len(Xp)):
        Xp[ii]= X[ii].copy() + epsx
        FNLP= HB_FNL(n, Xp, nfft, G, G_inv, qn, U)
        dFNLdX[:,ii:ii+1] =(FNLP-FNL)/epsx
        Xp=X.copy()

    #Calculate Residual Derivative Rx
    Rx = Z - dFNLdX # R = Z*X - Fnl = > R' = Z - dFnl/dX

    return Rx