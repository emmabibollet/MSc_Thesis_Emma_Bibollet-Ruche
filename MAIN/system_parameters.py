# =============================================================================
# AEROELASTIC SYSTEM PARAMETERS 
# =============================================================================

# --- System Dimensions ---
nmodes = 3       # pitch, plunge, flap
nstates = 6      # Wagner lag states
n = nmodes * 2 + nstates  # total number of degrees of freedom

# --- Mass and Geometry ---
mu = 100       # Mass ratio: m / (pi * rho * b^2)
ah = -0.5       # Elastic axis location
ch = [0.2,0.3, 0.4, 0.5]         # Hinge location
xa = 0.25        # CG offset Pitch
xbeta = [0.069785, 0.047836,0.03091, 0.018438]  # CG offset Flap
ra = 0.5         # Radius of gyration Pitch
rbeta = [0.170774, 0.132249,0.098459, 0.069488]   # Radius of gyration Flap
   
# --- Frequency Ratios ---
OMEGA1 = 0.5         # omega_h / omega_alpha (Plunge/Pitch)
OMEGA2 = [0, 0.1, 0.2, 0.3]       # omega_beta / omega_alpha (Flap/Pitch)
        
# --- Structural Damping ---
Zz = 0.002      # Damping Plunge
Za = 0.002      # Damping Pitch
c_beta = [0.02, 0.04, 0.06, 0.08] # Damping Flap

# --- Non-Linear Stiffnesses ---
k_h = 0
k_a = 0
GAMMA_beta = [10, 15, 20,25]