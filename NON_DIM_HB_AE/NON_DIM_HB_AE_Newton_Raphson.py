# HB Newton Raphson scheme function
# Adapted for the equations of motion of the aeroelastic equation: xdot = Qx + Fnl
# Resdidual becomes: R = xdot - Qx - Fnl = Z X - Fnl


# Imports -------------------------------------------------------------------------------
import numpy as np
import scipy

from NON_DIM_HB_AE.NON_DIM_HB_FNL import HB_FNL
from NON_DIM_HB_AE.NON_DIM_HB_AE_Rx import HB_Rx

np.set_printoptions(suppress=True, precision=4, linewidth=200)
np.set_printoptions(formatter={'float': '{: 0.2e}'.format})

# Start function---------------------------------------------------------------------------
def HB_Newton_Raphson(eps_NR, iter_NR, n, ntot, X, nfft, G, G_inv, qn, IH, Q, NABLA, Phase_Condition, w, U, precStab):
    #Description of inputs-----------------------------------------------------------------
    # eps_NR => Newton-Raphson resolution (set by user)
    # iter_NR => Maximum number of iterations of Newton-Raphson (set by user)
    # n ==> number of DOF
    # ntot => total size of the system
    # X ==> harmonic coefficients of the system DOF
    # nfft ==> number of FFT time points to determine the harmonic coefficients of Fnl
    # G ==> DFT matrix Gamma
    # G_inv ==> inverse DFT matrix Gamma-1
    # knl ==> non-linear stiffness in matrix form [knlh, knlalpha, knlbeta]
    # qn ==> non-linear force input distribution matrix
    # IH ==> Identitiy matrix of dimension n*(2*H+1) or ntot
    # Q ==> Linear state matrix
    # NABLA ==> block-diagonal derivative operator
    # Phase_Condition ==> Necessary condition to lock the periodic orbit dimensoin 1xntot
    # w ==> frequency (rad/s)
    # precStab ==> precision for stability analysis (set in HB_Parameters)
    
    eps = 1 #initial resolution
    iter = 0 #inital number of iterations
    
    while eps > eps_NR and iter < iter_NR:
        iter = iter+1

        #Calculate the non-linear force matrix: at an intial condition
        FNL= HB_FNL(n, X, nfft, G, G_inv, qn, U)

        # Update the linear matrix Z with the current w
        Z = w*np.kron(NABLA, np.eye(n)) - np.kron(IH,Q)

        #Calculate Residual R (error)
        R = (w*np.kron(NABLA, np.eye(n))  - np.kron(IH, Q)).dot(X) - FNL # R = Z*X - Fnl
        R_aug = np.vstack([R, Phase_Condition.dot(X)])

        #Calculate the derivative of the non-linear force matrix: forward differences method
        Rx = HB_Rx(X, ntot, n, nfft, G, G_inv, qn, FNL, Z,U)

        #Calculate the Residual Derivative Rw
        Rw = np.kron(NABLA, np.eye(n)).dot(X)

        #Build the corresponding Jacobian Matrix
        J  = np.vstack([np.hstack([Rx, Rw]), np.hstack([Phase_Condition, 0])])

        # Solve and update
        dxdw = -np.linalg.solve(J, R_aug)

        X = X + dxdw[0:ntot]
        w = w + dxdw[ntot]

        # 5. Check Convergence
        eps = scipy.linalg.norm(R)
        print("Iter:", iter, "||R|| =", eps, "||dX|| =", scipy.linalg.norm(dxdw[0:ntot]), "dw = ", dxdw[ntot])
    if eps < eps_NR:
        print('Convergence achieved in ', iter, 'iterations')
    else:
        print('Solution not converged and resolution is of', eps)

    
    # Post-Processing: Stability of the LCO --------------------------------------------
    HILL_matrix = - Rx
    ValP, VecP=np.linalg.eig(HILL_matrix)

    indices_tries = np.argsort((ValP.imag))  # Tri sur la partie imaginaire
    ValP_sorted = ValP[indices_tries]
    stabi_check=1
    for iii in range(len(ValP_sorted)):
        if np.real(ValP_sorted[iii])>=precStab:
            stabi_check=0

    return X, w, stabi_check






        
