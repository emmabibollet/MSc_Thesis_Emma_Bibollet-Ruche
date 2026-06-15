# Calculate the Harmonic Coefficients of the non-linear force:
# The equation of motion of the aeroelastic equation is: xdot = Qx + Fnl and Fnl = qn * f(hdot, alphadot, betadot, h, alpha, beta)

#Method: the harmonic coefficients of the non-linear force are computed by evaluating the system state x, xdot in the time domain and applying an FFT to the resulting force signal.
# Imports -------------------------------------------------------------------------------
import numpy as np

np.set_printoptions(suppress=True, precision=4, linewidth=200)
np.set_printoptions(formatter={'float': '{: 0.2e}'.format})

# Start function---------------------------------------------------------------------------
def HB_FNL(n, X, nfft, G, G_inv, qn, U):
    #Description of inputs-----------------------------------------------------------------
    # n ==> total size of the system (n_modes*2 + lag states = 12)
    # X ==> harmonic coefficients of the system DOF
    # nfft ==> number of FFT time points to determine the harmonic coefficients of Fnl
    # G ==> DFT matrix Gamma
    # G_inv ==> inverse DFT matrix Gamma-1
    # knl ==> non-linear stiffness in matrix form [knlh, knlalpha, knlbeta]
    # qn ==> non-linear force input distribution matrix

    # Back in time domain give the displacements (x):
    x = G.dot(X)
  
    idx = 0 
    fnl = np.zeros((len(x),1))
    
    # Loop through every time step
    while idx < nfft*n:
        xn= x[idx:idx+n] #extract x for each DOF at time i [hdot, alphadot, betadot, h, alpha, beta, w1 .... w6]
        xn_disp = xn[3:6] #extract the displacements [h, alpha, beta]
       
        #Calculate the f(hdot, alphadot, betadot, h, alpha, beta), in time domain x
        f = 1/(U**2)* (xn_disp**3)

        #Calculate Fnl in time domain
        fnl[idx:idx+n] = np.vstack([qn @ f , np.zeros((9,1))])
            
        #update index
        idx = idx + n

    #Convert back to frequency domain: inv DFT
    FNL = G_inv.dot(fnl)
    return FNL






