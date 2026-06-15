#Calculate the Residual Derivative Ru (relative to airspeed U)

# Imports -------------------------------------------------------------------------------
import numpy as np
from NON_DIM_AEROELASTIC_MODEL.NON_DIM_make_first_order_matrix import make_first_order_matrix

def HB_Ru(U, mu, C, D, E, F, W, W1, W2, M_inv, nmodes, w, NABLA, n, IH, Z, X, FNL):
    #Description of inputs-----------------------------------------------------------------
    # U ==> Airpseed
    # mu ==> nondimensional coefficient: mu = m/(pi*rho*b^2)
    # C ==> Damping Matrix
    # D ==> Aerodynamic Dampign Matrix
    # E ==> Stiffness Matrix
    # F ==> Aerodynamic Stiffness Matrix
    # W ==> Aerodynamic Influence Matrix
    # W1, W2 ==> Make Aerodynamic State Equation Matrices
    # M_inv ==> Inv Mass Matrix M=A+rho*B
    # nmodes ==> number of modes: pitch, plunge, flap-NES
    # w ==> omega (frequency in rad/s)
    # NABLA ==> derivative operator
    # n ==> number of DOF
    # IH ==> Identitiy matrix of dimension n*(2*H+1) or ntot
    # Z ==> Linear Forces Z = w NABLA kron In - IH kron Q
    # X ==> harmonic coefficients of the system

    epsx = 1e-6 
    Up = U + epsx
    Qp = make_first_order_matrix(mu, C, D, E , F , W, W1, W2, M_inv, Up, nmodes)
    Zp = w*np.kron(NABLA, np.eye(n)) - np.kron(IH,Qp)
    
    # 1. The linear part of the derivative
    dZdU_X = ((Zp - Z)/epsx).dot(X)
    
    # 2. The exact analytical derivative of the non-linear force w.r.t U
    dFNLdU = -(2.0 / U) * FNL
    
    # Total Residual derivative: R' = d(ZX)/dU - d(Fnl)/dU
    Ru = dZdU_X - dFNLdU

    return Ru