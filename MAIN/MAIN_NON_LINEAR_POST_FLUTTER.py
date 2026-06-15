# ====================================================================================
# MAIN EXECUTION FILE: Non-Linear Post-Flutter Response (Bifurcation Diagram)
# ====================================================================================

# !! This file must be run with the input system parameters defined in system_parameters.py and HBM parameters in HB_Parameters.py in this same folde !!

# File Structure:
# 1 - Obtain complete bifurcation diagram for 1 parameter combination of the Flap-NES (l.48)
# 2 - Parameter sweep (l.176):
    # a) SWEEP FLAP-NES Nonlinear Stiffness: GAMMA_beta (l.182)
    # b) SWEEP FLAP-NES Linear Damping: c_beta (l.211)
    # c) SWEEP FLAP-NES size of Flap: ch, xbeta, rbeta (l.241)
    # d) SWEEP FLAP-NES Linear Stiffness: OMEGA2 (l.271)


# Imports ---------------------------------------------------------------------
import numpy as np
from matplotlib import pyplot as plt
import os
import sys

# Dynamic Path Linking
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

from Plot_Configuration import*
from HB_Parameters import*
from system_parameters import*
from PLOT_HBM_Sweep import PLOT_HBM_Sweep # used for the parametric sweep
from NON_DIM_HBM_AE import NON_DIM_HBM_AE

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



################################################################################################
#  1- PLOT ENTIRE BIFURCATION DIAGRAM FOR A SINGLE SPECIFIC PARAMETRIC COMBINATION OF FLAP-NES
#################################################################################################

# !! NOTE !!: use return line 392 in NON_DIM_HBM_AE.py

# Create a figure with 2 subplots (Amplitude top, Frequency bottom), sharing the x-axis
fig, (ax1, ax_freq) = plt.subplots(nrows=2, ncols=1, figsize=(10, 10), sharex=True)
ax2_twin = ax1.twinx()  # Create secondary Y-axis for the top plot

boundary_u = 4.76 # 2-DOF Baseline linear flutter speed to beat

# The HBM scheme launches for every Hopf bifurcation: below define the number of Hopf bifurcations you wish to study (here 2)
for HB in range(2):
    if HB == 1:
        max_steps = 10 #increase or decrease the number of steps for each bifurcation branch under study
    
    # Launch the HBM scheme and chose the Flap-NES parameters below (either write the parameters or select the position in the system_parameters.py file):
    outputs = NON_DIM_HBM_AE(mu, OMEGA1, ah, ch[-1], xa, ra, Zz, Za, c_beta[0], 1.2 , xbeta[-1], rbeta[-1], 
                             k_h, k_a, GAMMA_beta[-1], n, nmodes, nstates, n_hbm, nfft, eps_NR, iter_NR, 
                             ds, dir_cont, avg_iter, ds_min, ds_max, Uf, max_steps, dU, dX, DOF_i, precStab, HB)
    
    # Unpack the specific data
    U_plot_t, heave_amp, pitch_amp, flap_amp, stabi_arr, U_NS, heave_NS, pitch_NS, flap_NS, U_LP, heave_LP, pitch_LP, flap_LP, freq_arr = outputs

    U_plot = np.array(U_plot_t)
    U_NS = np.array(U_NS) if len(U_NS) > 0 else []
    U_LP = np.array(U_LP) if len(U_LP) > 0 else []

    # ---------------------------------------------------------
    # LCO AMPLITUDE PLOT (ax1 / ax2_twin)
    # ---------------------------------------------------------

    l_style = '-'

    # Plot Stable Branches
    ax1.plot(U_plot, heave_amp, color=AERO_STYLE_NON_DIM['heave']['color'], linewidth=2, linestyle=l_style)
    ax2_twin.plot(U_plot, pitch_amp, color=AERO_STYLE_NON_DIM['pitch']['color'], linewidth=2, linestyle=l_style)
    ax2_twin.plot(U_plot, flap_amp, color=AERO_STYLE_NON_DIM['flap']['color'], linewidth=2, linestyle=l_style)
    
    # Plot Stability: Red Overlay for Unstable
    unstable_heave = np.where(stabi_arr == 0, heave_amp, np.nan)
    unstable_pitch = np.where(stabi_arr == 0, pitch_amp, np.nan)
    unstable_flap  = np.where(stabi_arr == 0, flap_amp, np.nan)
    
    ax1.plot(U_plot, unstable_heave, color='red', linewidth=2, linestyle=l_style)
    ax2_twin.plot(U_plot, unstable_pitch, color='red', linewidth=2, linestyle=l_style)
    ax2_twin.plot(U_plot, unstable_flap, color='red', linewidth=2, linestyle=l_style)
        
    # Plot Bifurcations
    if len(U_NS) > 0: 
        ax1.plot(U_NS, heave_NS, color='purple', marker='s', linestyle='None')
        ax2_twin.plot(U_NS, pitch_NS, color='purple', marker='s', linestyle='None')
        ax2_twin.plot(U_NS, flap_NS, color='purple', marker='s', linestyle='None')
    
    if len(U_LP) > 0: 
        ax1.plot(U_LP, heave_LP, color='black', marker='s', linestyle='None')
        ax2_twin.plot(U_LP, pitch_LP, color='black', marker='s', linestyle='None')
        ax2_twin.plot(U_LP, flap_LP, color='black', marker='s', linestyle='None')

    # ---------------------------------------------------------
    # FREQUENCY PLOT (ax_freq)
    # ---------------------------------------------------------
    unstable_freq = np.where(stabi_arr == 0, freq_arr, np.nan)
    ax_freq.plot(U_plot, freq_arr, color='blue', linewidth=2, linestyle=l_style)
    ax_freq.plot(U_plot, unstable_freq, color='red', linewidth=2, linestyle=l_style)

# =============================================================================
# Formatting Amplitude Plot (Top)
# =============================================================================
ax1.axvline(x=boundary_u, color='red', linestyle=':', linewidth=1.5)
ax1.axvspan(boundary_u, 5.0, color='red', alpha=0.1)

ax1.grid(True, linestyle='--', alpha=0.7)
ax1.set_ylabel(LABEL_AMP_HEAVE_NONDIM)
ax2_twin.set_ylabel(LABEL_AMP_ANGLE)
ax1.set_ylim(bottom=0)
ax2_twin.set_ylim(bottom=0)
ax1.set_xlim(right=5.0) 

# Simplified LEGEND for Amplitude
custom_lines = [
    plt.Line2D([0], [0], color=AERO_STYLE_NON_DIM['heave']['color'], lw=2),
    plt.Line2D([0], [0], color=AERO_STYLE_NON_DIM['pitch']['color'], lw=2),
    plt.Line2D([0], [0], color=AERO_STYLE_NON_DIM['flap']['color'], lw=2),
    plt.Line2D([0], [0], color='red', lw=2),
    plt.Line2D([0], [0], color='purple', marker='s', linestyle='None'),
    plt.Line2D([0], [0], color='black', marker='s', linestyle='None')
]
custom_labels = [
    AERO_STYLE_NON_DIM['heave']['label'],
    AERO_STYLE_NON_DIM['pitch']['label'],
    AERO_STYLE_NON_DIM['flap']['label'],
    'Unstable',
    'Neimark-Sacker (NS)',
    'Limit Point (LP)'
]
ax1.legend(custom_lines, custom_labels, loc='best')

# =============================================================================
# Formatting Frequency Plot (Bottom)
# =============================================================================
ax_freq.axvline(x=boundary_u, color='red', linestyle=':', linewidth=1.5)
ax_freq.axvspan(boundary_u, 5.0, color='red', alpha=0.1)
ax_freq.annotate(f'2-DOF Baseline\nU* = {boundary_u}', 
                 xy=(boundary_u, 0.95), xycoords=('data', 'axes fraction'),
                 xytext=(-10, 0), textcoords='offset points',
                 ha='right', va='top', color='red',
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="red", alpha=0.8))

ax_freq.grid(True, linestyle='--', alpha=0.7)
ax_freq.set_xlabel(LABEL_U_NONDIM)
ax_freq.set_ylabel(LABEL_FREQ_NONDIM)
ax_freq.set_ylim(bottom=0)
ax_freq.set_xlim(right=5.0) 

# Simplified Legend for frequency 
freq_legend_lines = [
    plt.Line2D([0], [0], color='blue', lw=2),
    plt.Line2D([0], [0], color='red', lw=2)
]
freq_legend_labels = ['Stable', 'Unstable']
ax_freq.legend(freq_legend_lines, freq_legend_labels, loc='best')

# Show/Render
plt.tight_layout()


# #############################################################################
# # 2 - Flap-NES Parametric Sweep
# #############################################################################

# The number of bifurcation branches studied is set to 1: therefore HB = 0 in all the following sweeps

# #=============================================================================
# # a) SWEEP FLAP-NES Nonlinear Stiffness: GAMMA_beta
# #=============================================================================
# # # !! NOTE !!: use return line 395 in NON_DIM_HBM_AE.py. This excludes frequency we only plot LCO amplitude evolution here

# sweep_results = {}
# for gb in GAMMA_beta:
#    print(f"\n{'='*60}")
#    print(f"EXECUTING CONTINUATION FOR GAMMA_beta = {gb}")
#    print(f"{'='*60}")
    
# #    Execute the function with the current gb value
#    outputs = NON_DIM_HBM_AE(mu, OMEGA1, ah, ch[1], xa, ra, Zz, Za, c_beta[0], OMEGA2[1], xbeta[1], rbeta[1], 
#                             k_h, k_a, gb, n, nmodes, nstates, n_hbm, nfft, eps_NR, iter_NR, 
#                             ds, dir_cont, avg_iter, ds_min, ds_max, Uf, max_steps, dU, dX, DOF_i, precStab, 0)
    
# #    Store the returned tuple in the dictionary
#    sweep_results[gb] = outputs

# PLOT_HBM_Sweep( 
#    sweep_results=sweep_results, 
#    param_label=r'$\Gamma_\beta$', 
#     AERO_STYLE_NON_DIM=AERO_STYLE_NON_DIM, 
#     LABEL_U_NONDIM_NORM=LABEL_U_NONDIM, 
#     LABEL_AMP_HEAVE_NONDIM=LABEL_AMP_HEAVE_NONDIM, 
#    LABEL_AMP_ANGLE=LABEL_AMP_ANGLE
# )


# #=============================================================================
# # b) SWEEP FLAP-NES Linear Damping: c_beta
# #=============================================================================

# # !! NOTE !!: use return line 395 in NON_DIM_HBM_AE.py. This excludes frequency we only plot LCO amplitude evolution here

# sweep_results = {}
# for cb in c_beta:
#     print(f"\n{'='*60}")
#     print(f"EXECUTING CONTINUATION FOR c_beta = {cb}")
#     print(f"{'='*60}")
    
#     # Execute the function with the current gb value
#     outputs = NON_DIM_HBM_AE(mu, OMEGA1, ah, ch[2], xa, ra, Zz, Za, cb, OMEGA2[1], xbeta[2], rbeta[2], 
#                              k_h, k_a, GAMMA_beta[0], n, nmodes, nstates, n_hbm, nfft, eps_NR, iter_NR, 
#                             ds, dir_cont, avg_iter, ds_min, ds_max, Uf, max_steps, dU, dX, DOF_i, precStab, 0)
    
#     # Store the returned tuple in the dictionary
#     sweep_results[cb] = outputs

# PLOT_HBM_Sweep(
#     sweep_results=sweep_results, 
#     param_label=r'$c_\beta$', 
#     AERO_STYLE_NON_DIM=AERO_STYLE_NON_DIM, 
#     LABEL_U_NONDIM_NORM=LABEL_U_NONDIM, 
#     LABEL_AMP_HEAVE_NONDIM=LABEL_AMP_HEAVE_NONDIM, 
#     LABEL_AMP_ANGLE=LABEL_AMP_ANGLE
# )


# # =============================================================================
# # c) SWEEP FLAP-NES size of Flap: ch, xbeta, rbeta
# # =============================================================================

# # !! NOTE !!: use return line 395 in NON_DIM_HBM_AE.py. This excludes frequency we only plot LCO amplitude evolution here

# sweep_results = {}
# for i, ch_val in enumerate(ch):
#    print(f"\n{'='*60}")
#    print(f"EXECUTING CONTINUATION FOR ch = {ch_val}")
#    print(f"{'='*60}")
    
# #    Execute the function with the current gb value
#    outputs = NON_DIM_HBM_AE(mu, OMEGA1, ah, ch_val, xa, ra, Zz, Za, c_beta[0], OMEGA2[1], xbeta[i], rbeta[i], 
#                             k_h, k_a, GAMMA_beta[0], n, nmodes, nstates, n_hbm, nfft, eps_NR, iter_NR, 
#                             ds, dir_cont, avg_iter, ds_min, ds_max, Uf, max_steps, dU, dX, DOF_i, precStab, 0)
    
# #   Store the returned tuple in the dictionary
#    sweep_results[ch_val] = outputs

# PLOT_HBM_Sweep(
#    sweep_results=sweep_results, 
#    param_label=r'$c_h$', 
#    AERO_STYLE_NON_DIM=AERO_STYLE_NON_DIM, 
#    LABEL_U_NONDIM_NORM=LABEL_U_NONDIM, 
#    LABEL_AMP_HEAVE_NONDIM=LABEL_AMP_HEAVE_NONDIM, 
#    LABEL_AMP_ANGLE=LABEL_AMP_ANGLE
# )


# #=============================================================================
# # d) SWEEP FLAP-NES Linear Stiffness: OMEGA2
# #=============================================================================

# # !! NOTE !!: use return line 395 in NON_DIM_HBM_AE.py. This excludes frequency we only plot LCO amplitude evolution here

# sweep_results = {}
# for O2 in OMEGA2:
#    print(f"\n{'='*60}")
#    print(f"EXECUTING CONTINUATION FOR OMEGA2 = {O2}")
#    print(f"{'='*60}")
    
# #    Execute the function with the current O2 value
#    outputs = NON_DIM_HBM_AE(mu, OMEGA1, ah, ch[0], xa, ra, Zz, Za, c_beta[0], O2, xbeta[0], rbeta[0], 
#                             k_h, k_a, GAMMA_beta[0], n, nmodes, nstates, n_hbm, nfft, eps_NR, iter_NR, 
#                             ds, dir_cont, avg_iter, ds_min, ds_max, Uf, max_steps, dU, dX, DOF_i, precStab, 0)
    
# #    Store the returned tuple in the dictionary
#    sweep_results[O2] = outputs

# PLOT_HBM_Sweep( 
#    sweep_results=sweep_results, 
#    param_label=r'$\Omega_2$', 
#     AERO_STYLE_NON_DIM=AERO_STYLE_NON_DIM, 
#     LABEL_U_NONDIM_NORM=LABEL_U_NONDIM, 
#     LABEL_AMP_HEAVE_NONDIM=LABEL_AMP_HEAVE_NONDIM, 
#    LABEL_AMP_ANGLE=LABEL_AMP_ANGLE
# )


plt.show()

