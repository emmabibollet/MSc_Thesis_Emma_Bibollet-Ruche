# Time-Integration for the NON-DIMENSIONAL AEROELASTIC MODEL: used to verify the HB results + observe the quasi-steady behaviour
import sys
import os
import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp

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
U =3.99 #Non Dimensional Airspeed

X0 = np.zeros(n)
X0[5] = np.deg2rad(5)

dt = 0.001 #Non-Dimensional Time
tfinal = 20000 #Non-Dimensional Time

#####################################################################################
# SYSTEM PARAMETERS: Unsteady Aerodynamics (Theodorsen)
#####################################################################################

# Obtain the Aeroelastic Matrices ---------------------------------------------------------------------------------------------------
model = Model_NON_DIM(mu=mu , OMEGA1 = OMEGA1, ah=ah, ch=ch[0], xa=xa, ra=ra, Zz=Zz, Za= Za, c_beta = c_beta[0], OMEGA2=1.2, xbeta=xbeta[0], rbeta=rbeta[0], k_h = k_h, k_a = k_a, GAMMA_beta = GAMMA_beta[-1])

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

# Define the steady-state region: 20% of the time vector (ignore the transient)
ss_idx = int(0.2 * len(sol.t))

# Calculate absolute maximums in the steady-state region
amp_h = np.max(np.abs(sol.y[3, ss_idx:]))               # Heave amplitude
amp_alpha = np.max(np.abs(np.degrees(sol.y[4, ss_idx:])))           # Pitch amplitude (radians)
amp_beta = np.max(np.abs(np.degrees(sol.y[5, ss_idx:])))            # Flap amplitude (radians)

# Convert rotations to degrees for output
amp_alpha_deg =amp_alpha
amp_beta_deg = amp_beta

print("-" * 50)
print(f"Steady-State LCO Amplitudes at U* = {U}")
print("-" * 50)
print(f"Plunge (h) Amplitude:     {amp_h:.4e}")
print(f"Pitch (\u03b1) Amplitude: {amp_alpha_deg:.4f}°")
print(f"Flap (\u03b2) Amplitude:  {amp_beta_deg:.4f}°")
print("-" * 50)

#####################################################################################
# Plotting
#####################################################################################

fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)

# Prepare the grouped data for the loop
t_plot = sol.t[ss_idx:]
y_plot = [
    sol.y[3, ss_idx:],             # Plunge (non-dimensional)
    np.degrees(sol.y[4, ss_idx:]), # Pitch (degrees)
    np.degrees(sol.y[5, ss_idx:])  # Flap (degrees)
]

# Formatting values for the text boxes
amp_text = [
    f'Max Amp: {amp_h:.4e}', 
    f'Max Amp: {amp_alpha:.4f}°', 
    f'Max Amp: {amp_beta:.4f}°'
]

# Keys mapping to the AERO_STYLE dictionary
dof_keys = ['heave', 'pitch', 'flap']
props = dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray')

# Execute the plotting loop
for i, ax in enumerate(axes):
    key = dof_keys[i]
    
    # Plot the line with the label assigned for the legend
    ax.plot(t_plot, y_plot[i], color=AERO_STYLE_NON_DIM[key]['color'], 
            linewidth=2, label=AERO_STYLE_NON_DIM[key]['label'])
    
    # Apply global units strictly to the Y-axis
    ax.set_ylabel(AERO_STYLE_NON_DIM[key]['unit'])
    ax.grid(True, linestyle='--', alpha=0.5)
    
    # Apply text box (top left)
    ax.text(0.02, 0.85, amp_text[i], transform=ax.transAxes, 
            fontsize=12, verticalalignment='top', bbox=props)
    
    # Generate the legend (top right to avoid text box)
    ax.legend(loc='upper right')
    
    # Apply symmetric limits with a 20% margin
    y_min, y_max = ax.get_ylim()
    max_y = max(abs(y_min), abs(y_max))
    ax.set_ylim(-max_y * 1, max_y * 1)

# Finalize the outer figure formatting
axes[-1].set_xlabel(LABEL_TIME_NONDIM) 
fig.suptitle(f'Non-Dimensional Time Integration (RK45) at U* = {U}')
fig.tight_layout()
plt.show()