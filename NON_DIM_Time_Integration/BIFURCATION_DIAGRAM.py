#Time-Integration: Bifurcation Diagram

# system parameters are defined in the MAIN folder under the system_parameters file
import sys
import os
import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp
import time
from scipy.signal import find_peaks

# Dynamic Path Linking
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

from MAIN.system_parameters import *
from NON_DIM_AEROELASTIC_MODEL.NON_DIM_Aero_Model_Pitch_Plunge_Control_Unsteady import Model_NON_DIM
from NON_DIM_AEROELASTIC_MODEL.NON_DIM_make_first_order_matrix import make_first_order_matrix
from Plot_Configuration import*

np.set_printoptions(suppress=True, precision=4, linewidth=200)
np.set_printoptions(formatter={'float': '{: 0.2e}'.format})

plt.rcParams.update({
    'font.size': 14,          # Global text size
    'axes.titlesize': 16,     # Title size
    'axes.labelsize': 16,     # Axis label size
    'xtick.labelsize': 14,    # X-tick size
    'ytick.labelsize': 14,    # Y-tick size
    'legend.fontsize': 14,    # Legend text size
    'font.family': 'serif',   # Professional serif font (like LaTeX)
    'lines.markersize': 10    # Make markers visible
})
#####################################################################################
# TIME INTEGRATION PARAMETERS
#####################################################################################
# Extract the airspeeds and initial conditions from the HBM simulation (copy-paste the resulting csv file here): 
#filename = os.path.join(current_dir, "NON-DIM-HB_Bifurcation_X0_B1.csv")
filename = os.path.join(current_dir, "NON-DIM-HB_Bifurcation_X0_B2.csv")

data = np.loadtxt(filename,  skiprows=1)

U_range = data[:, 0]        # Airspeed is the 1st column (Index 0)
X0_matrix = data[:, 5:]     # Initial conditions start at the 6th column (Index 5)

dt = 0.001 #Non-Dimensional Time
tfinal = 1000 #Non-Dimensional Time

# Record DATA
scatter_U_h = []
scatter_U_alpha = []
scatter_U_beta = []

scatter_h = []
scatter_alpha = []
scatter_beta = []

#####################################################################################
# SYSTEM PARAMETERS: Unsteady Aerodynamics (Theodorsen)
#####################################################################################

# Obtain the Aeroelastic Matrices ---------------------------------------------------------------------------------------------------
model = Model_NON_DIM(mu=mu , OMEGA1 = OMEGA1, ah=ah, ch=ch[0], xa=xa, ra=ra, Zz=Zz, Za= Za, c_beta = c_beta[0], OMEGA2=OMEGA2[0], xbeta=xbeta[0], rbeta=rbeta[0], k_h = k_h, k_a = k_a, GAMMA_beta = GAMMA_beta[0])

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

start_time = time.time()

for i, U in enumerate(U_range):


    X0 = X0_matrix[i, :]

    Q = make_first_order_matrix(mu, C, D, E , F , W, W1, W2, M_inv, U, nmodes)
    IE = model.make_initial_excitation(X0) #Initial Excitation 

    #####################################################################################
    # Time Integration: RK4
    #####################################################################################

    # Define the first-order ODE system
    def ODE(t, X, U):
        dX = np.dot(Q, X) + np.concatenate([qn @ (X[3:6]**3)*1/U**2, np.zeros(9)]) + IE(t)
        return dX

    # Set Time Span
    t_span = (0, tfinal)  
    t_eval = np.linspace(t_span[0], t_span[1], int((t_span[1]-t_span[0])/dt)+1)

    # Integrate using RK45
    sol = solve_ivp(ODE, t_span, X0, t_eval=t_eval, method='RK45', args=(U,), rtol=1e-9, atol=1e-9)

    #####################################################################################
    # Record Maximum Amplitudes (Steady-State)
    #####################################################################################
    ss_idx = int(0.8 * len(sol.t))
    
    # Extract steady-state signals
    ss_h = sol.y[3, ss_idx:]
    ss_alpha = sol.y[4, ss_idx:]
    ss_beta = sol.y[5, ss_idx:]

    # Find peaks independently for each DOF to capture their true phase-shifted maximums
    peaks_h_idx, _ = find_peaks(np.abs(ss_h))
    peaks_alpha_idx, _ = find_peaks(np.abs(ss_alpha))
    peaks_beta_idx, _ = find_peaks(np.abs(ss_beta))
    
    # Extract the true local maxima values
    vals_h = np.abs(ss_h[peaks_h_idx])
    vals_alpha = np.degrees(np.abs(ss_alpha[peaks_alpha_idx]))
    vals_beta = np.degrees(np.abs(ss_beta[peaks_beta_idx]))

    # Append to scatter lists (using .extend for cleaner array population)
    scatter_U_h.extend([U] * len(vals_h))
    scatter_h.extend(vals_h)

    scatter_U_alpha.extend([U] * len(vals_alpha))
    scatter_alpha.extend(vals_alpha)

    scatter_U_beta.extend([U] * len(vals_beta))
    scatter_beta.extend(vals_beta)


print(f"Total simulation time: {(time.time() - start_time) / 60:.2f} minutes")

#####################################################################################
# Extract and Save Raw Data
#####################################################################################

heave_data = np.column_stack((scatter_U_h, scatter_h, np.ones(len(scatter_h))))
pitch_data = np.column_stack((scatter_U_alpha, scatter_alpha, np.full(len(scatter_alpha), 2)))
flap_data  = np.column_stack((scatter_U_beta, scatter_beta, np.full(len(scatter_beta), 3)))

saved_data = np.vstack((heave_data, pitch_data, flap_data))

# Define the filename and column headers
filename = os.path.join(current_dir, "NON_DIM_Bifurcation_Raw_Data.txt")
header_titles = "U_star, Amplitude, DOF_ID (1=Heave 2=Pitch 3=Flap)"

# Save to a .txt file
np.savetxt(filename, saved_data, delimiter=",", header=header_titles, comments='', fmt='%.6f')

print(f"Raw data successfully saved to: {filename}")

#####################################################################################
# Plotting: Combined Bifurcation Diagram
#####################################################################################

fig1, ax1 = plt.subplots(figsize=(10, 6))
ax2_twin = ax1.twinx()  

# Define scatter plot styling parameters (no connecting lines, small transparent dots)
scatter_style = {'linestyle': 'none', 'marker': '.', 'markersize': 4, 'alpha': 0.6}

# Heave (Left Axis)
line1, = ax1.plot(scatter_U_h, scatter_h, 
                 color=AERO_STYLE_NON_DIM['heave']['color'], 
                 label=AERO_STYLE_NON_DIM['heave']['label'],
                 **scatter_style)

# Pitch & Flap (Right Axis)
line2, = ax2_twin.plot(scatter_U_alpha, scatter_alpha, 
                      color=AERO_STYLE_NON_DIM['pitch']['color'], 
                      label=AERO_STYLE_NON_DIM['pitch']['label'],
                      **scatter_style)
line3, = ax2_twin.plot(scatter_U_beta, scatter_beta, 
                      color=AERO_STYLE_NON_DIM['flap']['color'], 
                      label=AERO_STYLE_NON_DIM['flap']['label'],
                      **scatter_style)

# Axis Limits
ax1.set_ylim(bottom=0)
ax2_twin.set_ylim(bottom=0)

ax1.set_title("Brute-Force Bifurcation Diagram (Time-Integration Peaks)")
ax1.set_xlabel(LABEL_U_NONDIM) 

# Format Left Y-Axis (Heave)
ax1.set_ylabel(LABEL_AMP_HEAVE_NONDIM)
ax1.grid(True, linestyle='--', alpha=0.7)

# Format Right Y-Axis (Pitch & Flap)
ax2_twin.set_ylabel(LABEL_AMP_ANGLE)

# Unified Legend (Reconstruct lines for a cleaner legend appearance)
import matplotlib.lines as mlines
leg_heave = mlines.Line2D([], [], color=AERO_STYLE_NON_DIM['heave']['color'], marker='.', linestyle='none', markersize=10, label=AERO_STYLE_NON_DIM['heave']['label'])
leg_pitch = mlines.Line2D([], [], color=AERO_STYLE_NON_DIM['pitch']['color'], marker='.', linestyle='none', markersize=10, label=AERO_STYLE_NON_DIM['pitch']['label'])
leg_flap  = mlines.Line2D([], [], color=AERO_STYLE_NON_DIM['flap']['color'], marker='.', linestyle='none', markersize=10, label=AERO_STYLE_NON_DIM['flap']['label'])

ax1.legend(handles=[leg_heave, leg_pitch, leg_flap], loc='best')

fig1.tight_layout()
plt.show()

