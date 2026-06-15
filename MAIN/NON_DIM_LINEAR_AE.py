# ===================================================================================================================================
# Linear Stability Function: Hopf Bifurcation Detection and Futter Diagrams for each set of Flap-NES Parameters
# ====================================================================================================================================


# Imports ------------------------------------------------------------------------------------
import numpy as np
from scipy.optimize import linear_sum_assignment
import os
import sys

# Dynamic Path Linking
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

from NON_DIM_Linear_Analysis.flutter_calc_test import flutter_calc_test
from NON_DIM_AEROELASTIC_MODEL.NON_DIM_make_first_order_matrix import make_first_order_matrix
from NON_DIM_AEROELASTIC_MODEL.NON_DIM_Aero_Model_Pitch_Plunge_Control_Unsteady import Model_NON_DIM
from NON_DIM_HB_AE.NON_DIM_Hopf_Bifurcation import flutter_calc

np.set_printoptions(suppress=True, precision=4, linewidth=200)
np.set_printoptions(formatter={'float': '{: 0.2e}'.format})


#####################################################################################
# START FUNCTION
#####################################################################################

def NON_DIM_LINEAR_AE(mu, OMEGA1, ah, ch, xa, ra, Zz, Za, c_beta, OMEGA2, xbeta, rbeta, k_h, k_a, GAMMA_beta, Umax, nUmax, nmodes,nstates):
    #####################################################################################
    # SYSTEM PARAMETERS: Unsteady Aerodynamics (Theodorsen)
    #####################################################################################

    # Obtain the Aeroelastic Matrices ---------------------------------------------------
    # NON-DIM MODEL:
    model = Model_NON_DIM(mu=mu , OMEGA1 = OMEGA1, ah=ah, ch=ch, xa=xa, ra=ra, Zz=Zz, Za= Za, c_beta = c_beta, OMEGA2=OMEGA2, xbeta=xbeta, rbeta=rbeta, k_h = k_h, k_a = k_a, GAMMA_beta = GAMMA_beta)

    A = model.make_mass_matrix()
    B = model.make_aero_mass_matrix()
    C = model.make_damping_matrix()
    D = model.make_aerodynamic_damping_matrix()
    E = model.make_stiffness_matrix()
    F = model.make_aerodynamic_stiffness_matrix()
    W = model.make_aerodynamic_influence_matrix()
    W1, W2 = model.make_aerodynamic_state_equation_matrices()

    M = A + B/mu
    M_inv = np.linalg.inv(M)

    #####################################################################################
    # LINEAR VELOCITY SWEEP (Calculate Eigenvalues across Uvec)
    #####################################################################################
    Uvec =  np.linspace(0.5, Umax, nUmax)
    Zeta_matrix = np.zeros((3,nUmax))
    Freq_matrix =np.zeros((3,nUmax))

    tau_flut = np.zeros_like(Uvec)
    tau_div = np.zeros_like(Uvec)
    ncomp = np.zeros_like(Uvec)

    U_divergence = None
    divergence_found = False


    for i in range(nUmax):
        Q = make_first_order_matrix(mu, C, D, E , F , W, W1, W2, M_inv, Uvec[i], nmodes)

        tau_flut[i], tau_div[i], comp_eigs, eigis = flutter_calc_test(Q) 
        ncomp[i] = len(comp_eigs)

        if len(comp_eigs) >= 2 * nmodes:
            # Group conjugates by sorting imaginary parts descending
            iko = np.argsort(np.imag(comp_eigs))[::-1]
            sorted_comp_eigs = comp_eigs[iko]
            
            #Extract 3 modes
            unique_modes = sorted_comp_eigs[0:3]
            
            for k  in range(nmodes):
                omega = np.imag(unique_modes[k])
                zeta = -np.real(unique_modes[k]) / omega
                freq = omega / (2 * np.pi)*  Uvec[i]

                Zeta_matrix[k,i]=zeta
                Freq_matrix[k,i]=freq

    ####################################################################################
    # LOCATE HOPF SPEED WITH BISECTION ALGORITHM
    ####################################################################################
    nUmax= 500 #number of points

    Uflut, omega_flut, Q = flutter_calc(Umax, nUmax, nstates, nmodes, mu, C, D, E , F , W, W1, W2, M_inv)

    return Uvec, Zeta_matrix, Freq_matrix, Uflut, omega_flut