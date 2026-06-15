# NON DIMENSIONAL FUNCTION TO CALCULATE HOPF (flutter)  AIRSPEED FOR THE 2D UNSTEADY PITCH-PLUNGE-CONTROL WING USING SUCCESSIVE BISECTION METHOD
import numpy as np
from NON_DIM_Linear_Analysis.flutter_calc_test import flutter_calc_test
from NON_DIM_AEROELASTIC_MODEL.NON_DIM_make_first_order_matrix import make_first_order_matrix


def flutter_calc(Umax, nUmax, nstates, nmodes, mu, C, D, E , F , W, W1, W2, M_inv):
    
    #Select airspeed values in a range of interest
    Uvec =  np.linspace(0.5,Umax,nUmax)
    #Set up array for Hopf test function
    tau_flut = np.zeros_like(Uvec)
    #Set up array for pitchfork test function
    tau_div = np.zeros_like(Uvec)
    #Set up array for the eigenvalues
    eigis = np.zeros((nstates, nUmax))
    #Set up array for the number of complex eigenvalues
    ncomp = np.zeros_like(Uvec)

    for i in range(nUmax):
            
        Q = make_first_order_matrix(mu, C, D, E , F , W, W1, W2, M_inv, Uvec[i], nmodes)
        
        tau_flut[i], tau_div[i], comp_eigs, eigis = flutter_calc_test(Q) 
        ncomp[i] = len(comp_eigs)


    # Hopf Condition: Sign change indicates a crossing; ncomp check ensures the mode is still oscillatory. 
    # This ignores 'jumps' where eigenvalues collide on the real axis and transition from complex to real.
    iflut = np.where(((tau_flut[1:] * tau_flut[:-1]) < 0) & (ncomp[1:] == ncomp[:-1]))[0]

    # Static Divergenve Condition: sign change of tau_div indicates a purely real eigenvalue crossing zero
    idiv = np.where((tau_div[1:] * tau_div[:-1]) < 0)[0]

    ####################################################################################
    #Locate Hopf speed(s) and frequency(ies)
    ####################################################################################
    tol = 1e-10
    Uflut = np.zeros(np.size(iflut)) #Set array for all Hopf airspeeds
    omega_flut = np.zeros(np.size(iflut)) #Set array for all flutter frequencies


    for i in range(len(iflut)):
        #Calculate airspeed increment
        dU = (Uvec[iflut[i]+1] - Uvec[iflut[i]])/2
        Unow = Uvec[iflut[i]]
        while abs(dU) > tol:
            #Increase airspeed
            Unow = Unow+dU
            #New test function
            Q = make_first_order_matrix(mu, C, D, E , F , W, W1, W2, M_inv, Unow, nmodes)
            tau_flut_loc, tau_div_loc, comp_eigs_loc, eigis = flutter_calc_test(Q) 
            #Test for convergence
            if tau_flut_loc*tau_flut[iflut[i]] < 0:
                #If we have not converged decrease airspeed
                Unow = Unow - dU
                #Decrease airspeed step
                dU = dU/2
        
        #Store Hopf speed
        Uflut[i] = Unow


        #Store eigenvalues at airspeed
        Q = make_first_order_matrix(mu, C, D, E , F , W, W1, W2, M_inv, Unow, nmodes)
        tau_flut_loc, tau_div_loc, comp_eigs_loc, eigis = flutter_calc_test(Q) 
        criteigs = eigis[np.abs(np.real(eigis)) < 1e-4]
        omega_flut[i] = np.abs(criteigs[1])

    ####################################################################################
    # Locate Divergence speed(s)
    ####################################################################################
    Udiv = np.zeros(np.size(idiv)) # Set array for all Divergence airspeeds

    for i in range(len(idiv)):
        # Calculate airspeed increment
        dU = (Uvec[idiv[i]+1] - Uvec[idiv[i]])/2
        Unow = Uvec[idiv[i]]

        while abs(dU) > tol:
            # Increase airspeed
            Unow = Unow + dU
            # New test function
            Q = make_first_order_matrix(mu, C, D, E , F , W, W1, W2, M_inv, Unow, nmodes)
            tau_flut_loc, tau_div_loc, comp_eigs_loc, eigis = flutter_calc_test(Q) 

            # Test for convergence using the divergence tracking variable
            if tau_div_loc * tau_div[idiv[i]] < 0:
                # If we have not converged, decrease airspeed
                Unow = Unow - dU
                # Decrease airspeed step
                dU = dU / 2

        # Store Divergence speed
        Udiv[i] = Unow
        print(f"\n[Bisection Alert] Precise Static Divergence located at U* = {Unow:.5f}")
        
         
    return Uflut, omega_flut, Q

