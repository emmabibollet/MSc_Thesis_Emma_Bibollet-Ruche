# ===================================================================================================================================
# HBM & Stability Function: Calculate the Post-Flutter Dynamic Response and Stability Analysis for each set of Aeroelastic Parameters
# ====================================================================================================================================

# Imports ------------------------------------------------------------------------------------
import numpy as np
from ypstruct import struct
import scipy
from tqdm import tqdm
import os
import sys
import numpy as np

# Dynamic Path Linking
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))
if parent_dir not in sys.path:
    sys.path.append(parent_dir)


from NON_DIM_AEROELASTIC_MODEL.NON_DIM_Aero_Model_Pitch_Plunge_Control_Unsteady import Model_NON_DIM 
from NON_DIM_HB_AE.NON_DIM_Hopf_Bifurcation import flutter_calc
from NON_DIM_HB_AE.NON_DIM_HB_AE_Newton_Raphson import HB_Newton_Raphson
from NON_DIM_HB_AE.NON_DIM_HB_AE_Newton_Raphson_Cont import HB_Newton_Raphson_CONT
from NON_DIM_AEROELASTIC_MODEL.NON_DIM_make_first_order_matrix import make_first_order_matrix
from NON_DIM_HB_AE.NON_DIM_HB_FNL import HB_FNL
from NON_DIM_HB_AE.NON_DIM_HB_AE_Rx import HB_Rx
from NON_DIM_HB_AE.NON_DIM_HB_AE_Ru import HB_Ru

np.set_printoptions(suppress=True, precision=4, linewidth=200)
np.set_printoptions(formatter={'float': '{: 0.2e}'.format})

def NON_DIM_HBM_AE(mu, OMEGA1, ah, ch, xa, ra, Zz, Za, c_beta, OMEGA2, xbeta, rbeta, k_h, k_a, GAMMA_beta, n, nmodes, nstates, n_hbm, nfft, eps_NR, iter_NR, ds, dir_cont, avg_iter, ds_min, ds_max, Uf, max_steps, dU, dX, DOF_i, precStab, HB):

    #####################################################################################
    # SYSTEM PARAMETERS: Unsteady Aerodynamics (Theodersen)
    #####################################################################################

    # Obtain the Aeroelastic Matrices ---------------------------------------------------------------------------------------------------
    model = Model_NON_DIM(mu=mu , OMEGA1 = OMEGA1, ah=ah, ch=ch, xa=xa, ra=ra, Zz=Zz, Za= Za, c_beta = c_beta, OMEGA2=OMEGA2, xbeta=xbeta, rbeta=rbeta, k_h = k_h, k_a = k_a, GAMMA_beta = GAMMA_beta)

    A = model.make_mass_matrix()
    B = model.make_aero_mass_matrix()
    C = model.make_damping_matrix()
    D = model.make_aerodynamic_damping_matrix()
    E = model.make_stiffness_matrix()
    K = model.make_CUBIC_stiffness_matrix()
    F = model.make_aerodynamic_stiffness_matrix()
    W = model.make_aerodynamic_influence_matrix()
    W1,W2 = model.make_aerodynamic_state_equation_matrices()

    M=A+B/mu
    M_inv=np.linalg.inv(M)

    qn = -M_inv@K

    #####################################################################################
    # HB PARAMETERS: Defined in the HB_Parameters.py file
    #####################################################################################
    ntot = n*(2*n_hbm+1) #total dimension of the system

    # Storage of solutions -------------------------------------------------------------
    SOL = struct()

    #####################################################################################
    #CREATE DFT MATRIX: G, GAMMA
    #####################################################################################

    G = np.zeros((nfft, 2*n_hbm+1)) #empty gamma

    for i in range(nfft):
        wt = 2* np.pi * (i+1) /nfft #omega * ti
        G[i,0] = 1 #constant terms in the first column at all rows
        for j in range(1, 2*n_hbm, 2): #cos terms at from the 2nd column every two colums
            Hj = (j+1)/2 #harmonic order of each column, j
            G[i,j] = np.cos(Hj*wt)
        for k in range (2, 2*n_hbm +1, 2): #sin terms at from the 2nd column every two colums
            Hk = k/2 #harmonic order of each column, k
            G[i,k] = np.sin(Hk*wt)

    # Now for the inv Gamma matrix------------------------------------
    G_inv = np.zeros((2*n_hbm+1, nfft))  #empty gamma_inv

    for i in range(nfft):
        wt = 2* np.pi * (i+1) /nfft #omega * ti
        G_inv[0,i] = 1 #constant terms in the first column at all rows
        for j in range(1, 2*n_hbm, 2): #cos terms at from the 2nd column every two colums
            Hj = (j+1)/2 #harmonic order of each column, j
            G_inv[j,i] = 2*np.cos(Hj*wt)
        for k in range (2, 2*n_hbm +1, 2): #sin terms at from the 2nd column every two colums
            Hk = k/2 #harmonic order of each column, k
            G_inv[k,i] = 2*np.sin(Hk*wt)
    G_inv = G_inv/nfft


    # Modify G and G_inv when we have multiple DOF
    if (n>1):
        G = np.kron(G, np.eye(n))
        G_inv = np.kron(G_inv, np.eye(n))

    #####################################################################################
    # Other Useful Derivatives
    #####################################################################################

    # Define the Derivative Operator: NABLA ----------------------------------------------
    NABLA = np.zeros((2*n_hbm+2,2*n_hbm+2))
    for i in range(1,n_hbm+1):
        NABLA[2*i, 2*i+1] = i
        NABLA[2*i+1, 2*i] = -i
    NABLA= NABLA[1:2*n_hbm+2,1:2*n_hbm+2] #crop the extra row & column to get the correct size (2*n_hbm+1)

    #Define Identity Matrix, IH ----------------------------------------------------------
    IH = np.eye(2*n_hbm+1)

    #####################################################################################
    # INTIAL CONDITIONS: Hopf Bifurcation
    #####################################################################################
    nUmax= 300 #number of points

    Uflut, omega_flut, Q = flutter_calc(Uf, nUmax, nstates, nmodes, mu, C, D, E , F , W, W1, W2, M_inv)
    
    # If there are multiple: HB updates it represents the current Hopf bifurcation under study.
    Uflut = Uflut[HB]
    omega_flut = omega_flut[HB]

    # Store the exact Hopf Bifurcation point (The LCO Origin)
    SOL.w = [omega_flut.item()] 
    SOL.X = [np.zeros((ntot, 1))]  # The amplitude is strictly zero at the linear boundary
    SOL.U = [Uflut.item()] 
    SOL.norm = [0.0]
    SOL.stabi  = [0]

    # ESTIMATION OF SECOND POINT ON THE LIMIT CYCLE BRANCH -------------------------------------
    #Initial Airspeed
    Ui = Uflut + dir_cont*dU #Define the second point on the limit cycle branche at a slightly higher airspeed
    U = Ui

    # Initial Frequency
    w = omega_flut #Our initial guess for the second point frequency is the Hopf Bifurcation frequency

    # Constant Fourrier Coefficients
    X = np.zeros((ntot, 1)) #Initial =0
    X[n + DOF_i] = dX # non-zero starting guess X[n* constant terms [hdot alphadot betadot h alpha beta w1..w6] sin etc..]
    #X[:] = dX # non-zero starting guess X[n* constant terms [hdot alphadot betadot h alpha beta w1..w6] sin etc..]
    #Phase constraint
    Phase_Condition = np.zeros((ntot))
    Phase_Condition[n*2 + DOF_i] = 1 # Phase condition 1st cos of heave

    #Initial State Matrix Q at the flutter velocity
    Q = make_first_order_matrix(mu, C, D, E , F , W, W1, W2, M_inv, Ui, nmodes)

    # Correct the X and w value with the Newton-Raphson Scheme-----------------------------------
    [X, w, stabi_check] = HB_Newton_Raphson(eps_NR, iter_NR, n, ntot, X, nfft, G, G_inv, qn, IH, Q, NABLA, Phase_Condition, w, U, precStab)

    # Store intial conditions solutions ---------------------------------------------------
    SOL.w.append(w.item())
    SOL.X.append(X.copy())
    SOL.U.append(U.item())
    SOL.norm.append(float(scipy.linalg.norm(X)))
    SOL.stabi.append(stabi_check)
    SOL.bifurcations = []

    #####################################################################################
    # START THE CONTINUATION SCHEME
    #####################################################################################

    # Find the tangent vector, t, to the intial solution of X,w,U ------------------------------------------
    # Calculate the non-linear force matrix: at an intial condition
    FNL= HB_FNL(n, X, nfft, G, G_inv, qn, U)

    # Update the linear matrix Z with the current w
    Z = w*np.kron(NABLA, np.eye(n)) - np.kron(IH,Q)

    #Calculate the derivative of the non-linear force matrix: forward differences method
    Rx = HB_Rx(X, ntot, n, nfft, G, G_inv, qn, FNL, Z, U)
    #Calculate the Residual Derivative Rw
    Rw = np.kron(NABLA, np.eye(n)).dot(X)
    #Calculate the Residual Derivative Ru: forward differences method
    Ru = HB_Ru(U, mu, C, D, E, F, W, W1, W2, M_inv, nmodes, w, NABLA, n, IH, Z, X, FNL)

    #Build the corresponding Jacobian Matrix
    J  = np.vstack([np.hstack([Rx, Rw]), np.hstack([Phase_Condition, 0])])

    # Solve and update
    dxdw = -np.linalg.solve(J, np.vstack([Ru, [[0]]]))

    a: float = 1/np.sqrt(1 + dxdw.T.dot(dxdw)) #normalization condition to make the tangent vector unitary & we move 1 unit in U direction
    # t is the tangent vector under the form [dx; dw]
    dxdw1 = a*dir_cont*dxdw

    dU1 = dir_cont*a
    t = np.concatenate((dxdw1, dU1), axis=0)
    t = t / np.linalg.norm(t)

    p = 0
    iiter = 0

    # --- ADD PROGRESS BAR INITIALIZATION HERE ---
    # Calculate total airspeed distance to cover
    pbar = tqdm(total=max_steps, desc="Continuation Steps", bar_format="{l_bar}{bar}| {n:.2f}/{total:.2f} m/s [{elapsed}]")

    while p < max_steps:
        iter = 0
        eps = 1
        p += 1 

        # Record Previous DATA ------------------------------------------------------------------
        Xw = np.vstack([X, [[w.item()]]])
        Xwp = Xw
        Up = U

        # --- ADD PROGRESS BAR UPDATE HERE ---
        # Advance the bar by the exact amount U changed in this step
        pbar.update(1)
        # ------------------------------------
        
        # Estimation at distance ds ---------------------------------------------------------------
        Xw = Xw + ds*dxdw1
        U = U + ds*dU1

        # Correct Xw and U with a Newton-Raphson scheme ---------------------------------------------------------------
        X, w, U, eps, iter, stabi_check, HILL_matrix, ValP_sorted, VecP_sorted = HB_Newton_Raphson_CONT(eps_NR, iter_NR, n, nmodes, ntot, Xw[-1], Xw[:-1], nfft, G, G_inv, NABLA, mu, C, D, E, F, W, W1, W2, M_inv, qn, IH, Phase_Condition, dxdw1, dU1, U, precStab)

        # THE SOLUTION HAS CONVERGED :) -------------------------------------------------------------------------------
        if eps < eps_NR:

            # Reset iiter
            iiter = 0

            #Store the solutions
            SOL.w.append(w.item())
            SOL.X.append(X.copy())
            SOL.U.append(U.item())
            SOL.norm.append(float(scipy.linalg.norm(X)))
            SOL.stabi.append(stabi_check)

            # If there is a change in stability between two steps: Bifurcation Identification
            #if len(SOL.stabi) > 1 and SOL.stabi[-1] != SOL.stabi[-2]:
            if SOL.stabi[-1] != SOL.stabi[-2]:
                print("FLAG_Change_In_Stab")
                #Find the critical eigenvalue: which just crossed the Imaginary axis (smallest real part)
                # If it goes from stable to unstable: current step is unstable therefore Re(Lambda) > 0
                if SOL.stabi[-1] == 0:
                    critical_eig = min((eig for eig in ValP_sorted if np.real(eig) > precStab), key=np.real, default=None)
                    critical_vec = VecP_sorted[:, np.where(ValP_sorted == critical_eig)[0][0]] if critical_eig is not None else None 
                # Else if the reponse goes from unstable to stable: current step is unstable therefore Re(Lambda) < 0
                elif SOL.stabi[-1] == 1:
                    critical_eig = max((eig for eig in ValP_sorted if np.real(eig) < -precStab), key=np.real, default=None)
                    critical_vec = VecP_sorted[:, np.where(ValP_sorted == critical_eig)[0][0]] if critical_eig is not None else None
                # DISTINGUISH BIFURCATION TYPE
                if np.abs(np.imag(critical_eig)) <= precStab:
                    bif_type = "LP"
                    print("LP Bifurcation DETECTED !")
                elif np.abs(np.imag(critical_eig) - (w / 2.0)) <= precStab:
                    bif_type = "PD"
                    print("PD Bifurcation DETECTED!")
                else:
                    bif_type = "NS"
                    print("NS Bifurcation DETECTED !")

                # Storage  
                xt_bif = G.dot(X).flatten() 
                
                # Extract max amplitudes (Heave = index 3, Pitch = index 4, Flap = index 5)
                xt_heave = xt_bif[3::n]
                heave_bif_amp = float((np.max(xt_heave) - np.min(xt_heave)) / 2)
                xt_pitch = xt_bif[4::n]
                pitch_bif_amp = float(np.rad2deg((np.max(xt_pitch) - np.min(xt_pitch)) / 2)) # Converted to deg
                xt_flap = xt_bif[5::n]
                flap_bif_amp = float(np.rad2deg((np.max(xt_flap) - np.min(xt_flap)) / 2))  # Converted to deg

                SOL.bifurcations.append({
                    'type': bif_type,
                    'Airspeed': U.item(),
                    'heave_amp': heave_bif_amp,
                    'pitch_amp': pitch_bif_amp,
                    'flap_amp': flap_bif_amp
                })

            # Prepare for the next point -----------------------------------------------------------
            # Calculate the non-linear force matrix: at an intial condition
            FNL= HB_FNL(n, X, nfft, G, G_inv, qn, U)

            # Update the linear matrix Z with the current w, U
            Q = make_first_order_matrix(mu, C, D, E , F , W, W1, W2, M_inv, U, nmodes)
            Z = w*np.kron(NABLA, np.eye(n)) - np.kron(IH,Q)

            #Calculate the derivative of the non-linear force matrix: forward differences method
            Rx = HB_Rx(X, ntot, n, nfft, G, G_inv, qn, FNL, Z, U)
            #Calculate the Residual Derivative Rw
            Rw = np.kron(NABLA, np.eye(n)).dot(X)
            #Calculate the Residual Derivative Ru: forward differences method
            Ru = HB_Ru(U, mu, C, D, E, F, W, W1, W2, M_inv, nmodes, w, NABLA, n, IH, Z, X, FNL)

            #Build the corresponding Jacobian Matrix
            J  = np.vstack([np.hstack([Rx, Rw]), np.hstack([Phase_Condition, 0])])

            # Solve and update
            dxdw = -np.linalg.solve(J, np.vstack([Ru, [[0]]]))

            a: float = 1/np.sqrt(1 + dxdw.T.dot(dxdw)) #normalization condition to make the tangent vector unitary & we move 1 unit in U direction
            
            # To ensure a forward progress calculate the dot product of the previous tangent and the new candidate tangent (dx)
            if p!=1:
                a = a * np.sign(tgp.T.dot(np.concatenate(((dxdw,[[1]])),axis=0)))
            else:
                a = a *dir_cont
        
            #Update the size of the continuation step ds
            update_coef = 2**((avg_iter - iter)/2)
            ds = ds *update_coef
            if ds > ds_max:
                ds = ds_max
            if ds < ds_min:
                ds = ds_min

            # t is the tangent vector under the form [dx; dw]
            dxdw1 = a*dxdw
            dU1 = a
            tgp = np.concatenate((dxdw1, dU1), axis=0)
            tgp = tgp / np.linalg.norm(tgp)

            #print('Step:', p, ', w=', SOL.w[p], ', amp=', SOL.norm[p], ', U=', SOL.U[p], ', res=', eps, ', iter=', iter, ', ds=', ds)

        # THE SOLUTION HAS NOT CONVERGED :( ----------------------------------------------------------
        else :
            iiter  = iiter + 1
            if iiter > 3:
                print('REPEATED NON-CONVERGED SOLUTION NEED TO STOP')
                break
            if eps > eps_NR:
                ds = ds/4
                print('The continuation step is reduced as the Nb of maximum iteration of the NR scheme has been reached')
                if (ds < ds_min) :
                    ds = ds_min
            Xw = Xwp
            U = Up
            p = p-1
            iter = 0

    pbar.close()
    print('Continuation Scheme Finished :)')    

    #####################################################################################
    # SAVE_FINAL_DATA
    #####################################################################################

    # EXTRACT the LCO amplitude of each DOF & Phase Angles ----------------------------------------------
    A = [[] for _ in range(n)]

    for X_sol in SOL.X:
        #Go back into time domain: get x for each DOF over an entire cycle
        xt = G.dot(X_sol).flatten() #xt has a size (n*nfft,1)

        # Extract for each DOF, the maximum amplitude during the cycle 
        for k in range(n):
            xt_n = xt[k::n]
            An = (np.max(xt_n) - np.min(xt_n)) / 2
            A[k].append(An)

    # Clean the common X-axis data (Airspeed) to prevent array shape errors
    U_plot = [float(u) for u in SOL.U]

    # Extract and convert the exact data arrays
    heave_amp = [float(val) for val in A[3]]          # Non-dimensional
    pitch_amp = [np.rad2deg(float(val)) for val in A[4]] # Converted to degrees
    flap_amp  = [np.rad2deg(float(val)) for val in A[5]] # Converted to degrees
    stabi_arr = np.array(SOL.stabi)

    # EXTRACT BIFURCATION POINTS VIA INDEX MAPPING --------------------------------------
    
    U_NS, heave_NS, pitch_NS, flap_NS = [], [], [], []
    U_LP, heave_LP, pitch_LP, flap_LP = [], [], [], []

    for bif in SOL.bifurcations:
        
        if bif['type'] == 'NS':
            U_NS.append(bif['Airspeed'])
            heave_NS.append(bif['heave_amp'])
            pitch_NS.append(bif['pitch_amp'])
            flap_NS.append(bif['flap_amp'])
            
        elif bif['type'] == 'LP':
            U_LP.append(bif['Airspeed'])
            heave_LP.append(bif['heave_amp'])
            pitch_LP.append(bif['pitch_amp'])
            flap_LP.append(bif['flap_amp'])

    w_plot = [float(w) / (2 * np.pi) for w in SOL.w]

# USE THIS FOR COMPLETE BIFURCATION DIAGRAM (1): to include frequency
#    return U_plot, heave_amp, pitch_amp, flap_amp, stabi_arr, U_NS, heave_NS, pitch_NS, flap_NS, U_LP, heave_LP, pitch_LP, flap_LP, w_plot

# USE THIS FOR PARAMETRIC SWEEP (2): to exclude frequency
    return U_plot, heave_amp, pitch_amp, flap_amp, stabi_arr, U_NS, heave_NS, pitch_NS, flap_NS, U_LP, heave_LP, pitch_LP, flap_LP

