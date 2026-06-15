# HB Newton Raphson scheme function for continuation: here we included the re-calculation of Z since it depends on w

# Imports -------------------------------------------------------------------------------
import numpy as np
import scipy

from NON_DIM_HB_AE.NON_DIM_HB_FNL import HB_FNL
from NON_DIM_AEROELASTIC_MODEL.NON_DIM_make_first_order_matrix import make_first_order_matrix
from NON_DIM_HB_AE.NON_DIM_HB_AE_Rx import HB_Rx
from NON_DIM_HB_AE.NON_DIM_HB_AE_Ru import HB_Ru

np.set_printoptions(suppress=True, precision=4, linewidth=200)
np.set_printoptions(formatter={'float': '{: 0.2e}'.format})

# Start function---------------------------------------------------------------------------
def HB_Newton_Raphson_CONT(eps_NR, iter_NR, n, nmodes, ntot, w, X, nfft, G, G_inv, NABLA, mu, C, D, E, F, W, W1, W2, M_inv, qn, IH, Phase_Condition, dxdw, dU, U, precStab):
    #Description of inputs-----------------------------------------------------------------
    # eps_NR => Newton-Raphson resolution (set by user)
    # iter_NR => Maximum number of iterations of Newton-Raphson (set by user)
    # n ==> number of DOF
    # nmodes ==> number of modes: pitch, plunge, flap-NES
    # ntot ==> total dimension of the system ntot = n*(2*n_hbm+1)
    # w ==> omega (frequency in rad/s)
    # X ==> harmonic coefficients of the system
    # nfft ==> number of FFT time points to determine the harmonic coefficients of Fnl
    # G ==> DFT matrix Gamma
    # G_inv ==> inverse DFT matrix Gamma-1
    # NABLA ==> derivative operator
    # rho ==> density
    # C ==> Damping Matrix
    # D ==> Aerodynamic Dampign Matrix
    # E ==> Stiffness Matrix
    # F ==> Aerodynamic Stiffness Matrix
    # W ==> Aerodynamic Influence Matrix
    # W1, W2 ==> Make Aerodynamic State Equation Matrices
    # M_inv ==> Inv Mass Matrix M=A+rho*B
    # qn ==> non-linear force input distribution matrix
    # IH ==> Identitiy matrix of dimension n*(2*H+1) or ntot
    # Phase_Condition ==> Necessary condition to lock the periodic orbit dimensoin 1xntot
    # dxdw ==> [Delta X Delta w] (continuation scheme direction)
    # dU ==> Delta U (for continuation scheme direction)
    # U ==> Airpseed
    # precStab ==> precision for stability analysis (set in HB_Parameters)
    
    eps = 1 #initial resolution
    iter = 0 #inital number of iterations

    # Estimate the linear state matrix, Z, ------------------------------------------------
    Q = make_first_order_matrix(mu, C, D, E , F , W, W1, W2, M_inv, U, nmodes)
    Z = w*np.kron(NABLA, np.eye(n)) - np.kron(IH,Q)
    #Calculate the non-linear force matrix & derivatives: at an intial condition ----------
    FNL= HB_FNL(n, X, nfft, G, G_inv, qn, U)
        
    # Build the JACOBIAN MATRIX JJ ------------------------------------------------------------
        
    #Calculate Residual Derivative Rx
    Rx = HB_Rx(X, ntot, n, nfft, G, G_inv, qn, FNL, Z,U)
       
    #Calculate the Residual Derivative Rw
    Rw = np.kron(NABLA, np.eye(n)).dot(X) 
       
    #Calculate the Residual Derivative Ru
    Ru = HB_Ru(U, mu, C, D, E, F, W, W1, W2, M_inv, nmodes, w, NABLA, n, IH, Z, X, FNL)

    #Build Jacobian, JJ
    top_row=np.hstack((Rx, Rw, Ru))
    middle_row = np.hstack((dxdw.T, dU))
    bottom_row = np.hstack((Phase_Condition, 0, 0))
    JJ=np.vstack((top_row, middle_row, bottom_row))
    

    #Build RESIDUAL MATRIX RR ---------------------------------------------------------------
    #Calculate Residual R (error)
    R = (w*np.kron(NABLA, np.eye(n))  - np.kron(IH, Q)).dot(X) - FNL # R = Z*X - Fnl & Z = NABLA In - Ih kron Q
    # TOTAL RESIDUAL RR
    RR = np.vstack((R,[[0]],[[Phase_Condition.dot(X).item()]]))

    
    while eps > eps_NR and iter < iter_NR:
        iter = iter+1

        # Solve and update
        correction = -np.linalg.solve(JJ, RR)
        X = X + correction[0:-2]
        w = w + correction[-2].item()
        U = U + correction[-1].item()

        # Estimate the linear state matrix, Z, ------------------------------------------------
        Q = make_first_order_matrix(mu, C, D, E , F , W, W1, W2, M_inv, U, nmodes)
        Z = w*np.kron(NABLA, np.eye(n)) - np.kron(IH,Q)
        #Calculate the non-linear force matrix & derivatives: at an intial condition ----------
        FNL= HB_FNL(n, X, nfft, G, G_inv, qn, U)
        
        # Build the JACOBIAN MATRIX JJ ------------------------------------------------------------
        
        #Calculate Residual Derivative Rx
        Rx = HB_Rx(X, ntot, n, nfft, G, G_inv, qn, FNL, Z,U)
       
        #Calculate the Residual Derivative Rw
        Rw = np.kron(NABLA, np.eye(n)).dot(X) 
       
       #Calculate the Residual Derivative Ru
        Ru = HB_Ru(U, mu, C, D, E, F, W, W1, W2, M_inv, nmodes, w, NABLA, n, IH, Z, X, FNL)


        #Build Jacobian, JJ
        top_row=np.hstack((Rx, Rw, Ru))
        middle_row = np.hstack((dxdw.T, dU))
        bottom_row = np.hstack((Phase_Condition, 0, 0))
        JJ=np.vstack((top_row, middle_row, bottom_row))
    

        #Build RESIDUAL MATRIX RR ---------------------------------------------------------------
        #Calculate Residual R (error)
        R = (w*np.kron(NABLA, np.eye(n))  - np.kron(IH, Q)).dot(X) - FNL # R = Z*X - Fnl & Z = NABLA In - Ih kron Q
        # TOTAL RESIDUAL RR
        RR = np.vstack((R,[[0]],[[Phase_Condition.dot(X).item()]]))

        # Check Convergence
        eps = scipy.linalg.norm(RR)

    # Post-Processing: Stability of the LCO --------------------------------------------
    HILL_matrix =  - Rx
    ValP, VecP=np.linalg.eig(HILL_matrix)

    indices_tries = np.argsort(np.abs((ValP.imag)))  
    ValP_sorted = ValP[indices_tries][0:n]
    VecP_sorted = VecP[:, indices_tries][:,0:n]
    
    # Stability check
    stabi_check = 1
    for iii in range(len(ValP_sorted)):
        if np.real(ValP_sorted[iii]) >= precStab:
            stabi_check = 0


    return X,w, U, eps, iter, stabi_check, HILL_matrix, ValP_sorted, VecP_sorted






        
