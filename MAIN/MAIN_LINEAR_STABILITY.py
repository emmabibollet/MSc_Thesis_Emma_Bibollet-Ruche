# ====================================================================================
# MAIN EXECUTION FILE: Linear Stability Analysis (Flutter Diagrams)
# ====================================================================================

# !! This file must be run with the input system parameters defined in system_parameters.py and HBM parameters in HB_Parameters.py in this same folde !!

# File Structure:
# 1 - Flutter diagram of the baseline 2-DOF (locked Flap-NES)(l.50)
# 2 - Futter diagram for 1 parameter combination of the Flap-NES (l.150)
# 3 - Flutter diagram parameter sweep (l.255):
    # a) SWEEP FLAP-NES Linear Damping: c_beta (l.259)
    # b) SWEEP FLAP-NES Linear Stiffness: OMEGA2 (l.282)
    # c) SWEEP FLAP-NES size of Flap: ch, xbeta, rbeta (l.305)


# Imports ---------------------------------------------------------------------
import numpy as np
from matplotlib import pyplot as plt
import os
import sys
import matplotlib.lines as mlines

# Dynamic Path Linking
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

from Plot_Configuration import*
from HB_Parameters import*
from system_parameters import*
from PLOT_Linear_Sweep import PLOT_Linear_Sweep
from NON_DIM_LINEAR_AE_2DOF import NON_DIM_LINEAR_AE_2DOF
from NON_DIM_LINEAR_AE import NON_DIM_LINEAR_AE


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

# ################################################################################
# # 1 - Flutter Diagram: Baseline 2-DOF system (locked Flap-NES)
# ###############################################################################
# print(f"\n{'='*60}")
# print("CALCULATING BASELINE 2-DOF WING FLUTTER SPEED")
# print(f"{'='*60}")
# nmodes = 2
# nstates = 4

# outputs = NON_DIM_LINEAR_AE_2DOF(mu, OMEGA1, ah, xa, ra, Zz, Zz, k_h, k_a, Uf, 100, nmodes,nstates)

# # =============================================================================
# # UNPACK & PREPARE DATA
# # =============================================================================
# # Unpack the returned tuple from NON_DIM_LINEAR_AE
# Uvec, Zeta_matrix, Freq_matrix, Uflut, omega_flut = outputs

# boundary_u = 4.76  # 2-DOF Baseline boundary

# # Mapping style dictionary entries to match mode order [Heave, Pitch, Flap]
# mode_keys = ['pitch', 'heave']

# custom_lines = []
# custom_labels = []

# # =============================================================================
# # INITIALIZE COMBINED FIGURE (Stacked Subplots Sharing X-Axis)
# # =============================================================================
# fig, (ax_damp, ax_freq) = plt.subplots(nrows=2, ncols=1, figsize=(10, 12), sharex=True)

# # Plot each degree of freedom / mode
# for k in range(nmodes):
#     key = mode_keys[k % len(mode_keys)]
#     mc = AERO_STYLE_NON_DIM[key]['color']
#     label = AERO_STYLE_NON_DIM[key]['label']
    
#     # Assign markers sequentially
#     ms = ['o', 's', '^'][k % 3]
    
#     # Top Plot: Damping points
#     ax_damp.plot(Uvec, Zeta_matrix[k, :], color=mc, marker=ms, linestyle='None', markersize=5)
#     # Bottom Plot: Frequency points
#     ax_freq.plot(Uvec, Freq_matrix[k, :], color=mc, marker=ms, linestyle='None', markersize=5)
    
#     # Save elements for a clean unified legend
#     custom_lines.append(mlines.Line2D([0], [0], color=mc, marker=ms, linestyle='None'))
#     custom_labels.append(label)

# # Mark Flutter crossings on the Damping plot (where zeta crosses 0)
# if len(Uflut) > 0:
#     for u_val in np.atleast_1d(Uflut):
#         if u_val > 0:
#             ax_damp.plot(u_val, 0, color='red', marker='x', markeredgewidth=3, zorder=5)


# # =============================================================================
# # FORMATTING: DAMPING PLOT (TOP AXIS)
# # =============================================================================
# ax_damp.axhline(0, color='black', linewidth=1)

# # Baseline Boundary Indicators
# ax_damp.axvline(x=boundary_u, color='red', linestyle=':', linewidth=1.5)
# ax_damp.axvspan(boundary_u, Uvec[-1], color='red', alpha=0.1)

# ax_damp.set_ylabel(r'Damping Ratio $\zeta$')
# ax_damp.grid(True, linestyle='--', alpha=0.7)
# ax_damp.set_ylim(-0.2, 0.5)
# ax_damp.set_xlim(Uvec[0], Uvec[-1])
# ax_damp.legend(custom_lines, custom_labels, loc='upper right', bbox_to_anchor=(0.98, 0.82))

# # =============================================================================
# # FORMATTING: FREQUENCY PLOT (BOTTOM AXIS)
# # =============================================================================
# ax_freq.axhline(0, color='black', linewidth=1)

# # Baseline Boundary Indicators
# ax_freq.axvline(x=boundary_u, color='red', linestyle=':', linewidth=1.5)
# ax_freq.axvspan(boundary_u, Uvec[-1], color='red', alpha=0.1)
# # ax_freq.annotate(f'2-DOF Baseline\n{LABEL_U_NONDIM} = {boundary_u}', 
# #                  xy=(boundary_u, 0.90), xycoords=('data', 'axes fraction'),
# #                  xytext=(-10, 0), textcoords='offset points',
# #                  ha='right', va='top', color='red',
# #                  bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="red", alpha=0.8))

# ax_freq.set_xlabel(LABEL_U_NONDIM)
# ax_freq.set_ylabel(LABEL_FREQ_NONDIM)
# ax_freq.grid(True, linestyle='--', alpha=0.7)
# ax_freq.set_ylim(bottom=0)
# ax_freq.set_xlim(Uvec[0], 5)

# # =============================================================================
# # MAIN FIGURE SETTINGS & DISPLAY
# # =============================================================================
# # Apply main global title across both figures without subplot subtitles
# # fig.suptitle(r"2-DOF Baseline: Flap-NES Locked")

# plt.tight_layout()
# # Adjust layout to cleanly incorporate the main suptitle space
# fig.subplots_adjust(top=0.93)


################################################################################
# 2 - Futter diagram for 1 parameter combination of the Flap-NES 
###############################################################################

# # Select the parameters of the system under study:
# outputs = NON_DIM_LINEAR_AE(mu, OMEGA1, ah, ch[-1], xa, ra, Zz, Za, 
#                                 c_beta[0], 1.2, xbeta[-1],rbeta[-1], 
#                                 k_h, k_a, GAMMA_beta[-1], Uf, 100, nmodes, nstates)
    

# # =============================================================================
# # UNPACK & PREPARE DATA
# # =============================================================================
# # Unpack the returned tuple from NON_DIM_LINEAR_AE
# Uvec, Zeta_matrix, Freq_matrix, Uflut, omega_flut = outputs

# boundary_u = 4.76  # 2-DOF Baseline boundary

# # Mapping style dictionary entries to match mode order [Heave, Pitch, Flap]
# mode_keys = ['pitch', 'heave', 'flap']

# custom_lines = []
# custom_labels = []

# # =============================================================================
# # INITIALIZE COMBINED FIGURE (Stacked Subplots Sharing X-Axis)
# # =============================================================================
# fig, (ax_damp, ax_freq) = plt.subplots(nrows=2, ncols=1, figsize=(10, 12), sharex=True)

# # Plot each degree of freedom / mode
# for k in range(nmodes):
#     key = mode_keys[k % len(mode_keys)]
#     mc = AERO_STYLE_NON_DIM[key]['color']
#     label = AERO_STYLE_NON_DIM[key]['label']
    
#     # Assign markers sequentially
#     ms = ['o', 's', '^'][k % 3]
    
#     # Top Plot: Damping points
#     ax_damp.plot(Uvec, Zeta_matrix[k, :], color=mc, marker=ms, linestyle='None', markersize=5)
#     # Bottom Plot: Frequency points
#     ax_freq.plot(Uvec, Freq_matrix[k, :], color=mc, marker=ms, linestyle='None', markersize=5)
    
#     # Save elements for a clean unified legend
#     custom_lines.append(mlines.Line2D([0], [0], color=mc, marker=ms, linestyle='None'))
#     custom_labels.append(label)

# # Mark Flutter crossings on the Damping plot (where zeta crosses 0)
# if len(Uflut) > 0:
#     for u_val in np.atleast_1d(Uflut):
#         if u_val > 0:
#             ax_damp.plot(u_val, 0, color='red', marker='x', markeredgewidth=3, zorder=5)

# # Generate Text Box Info for Flutter Speed
# flutter_text_lines = []
# if len(Uflut) > 0:
#     u_str = ", ".join([f"{u:.3f}" for u in Uflut if u > 0])
#     flutter_text_lines.append(f"Hopf Bifurcation: {LABEL_U_NONDIM}$_f = {u_str}$")
# else:
#     flutter_text_lines.append("System: Stable")

# flutter_text_block = "\n".join(flutter_text_lines)
# box_props = dict(boxstyle='square,pad=0.5', facecolor='white', alpha=0.9, edgecolor='black')

# # Place Flutter Box in upper right corner of the top panel
# ax_damp.text(0.95, 0.95, flutter_text_block, transform=ax_damp.transAxes, 
#              verticalalignment='top', horizontalalignment='right', bbox=box_props, zorder=10)

# # =============================================================================
# # FORMATTING: DAMPING PLOT (TOP AXIS)
# # =============================================================================
# ax_damp.axhline(0, color='black', linewidth=1)

# # Baseline Boundary Indicators
# ax_damp.axvline(x=boundary_u, color='red', linestyle=':', linewidth=1.5)
# ax_damp.axvspan(boundary_u, 5, color='red', alpha=0.1)

# ax_damp.set_ylabel(r'Damping Ratio $\zeta$')
# ax_damp.grid(True, linestyle='--', alpha=0.7)
# ax_damp.set_ylim(-0.2, 0.5)
# ax_damp.set_xlim(Uvec[0], 5)
# ax_damp.legend(custom_lines, custom_labels, loc='upper right', bbox_to_anchor=(0.98, 0.82))

# # =============================================================================
# # FORMATTING: FREQUENCY PLOT (BOTTOM AXIS)
# # =============================================================================
# ax_freq.axhline(0, color='black', linewidth=1)

# # Baseline Boundary Indicators
# ax_freq.axvline(x=boundary_u, color='red', linestyle=':', linewidth=1.5)
# ax_freq.axvspan(boundary_u, Uvec[-1], color='red', alpha=0.1)
# ax_freq.annotate(f'2-DOF Baseline\n{LABEL_U_NONDIM} = {boundary_u}', 
#                  xy=(boundary_u, 0.90), xycoords=('data', 'axes fraction'),
#                  xytext=(-10, 0), textcoords='offset points',
#                  ha='right', va='top', color='red',
#                  bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="red", alpha=0.8))

# ax_freq.set_xlabel(LABEL_U_NONDIM)
# ax_freq.set_ylabel(LABEL_FREQ_NONDIM)
# ax_freq.grid(True, linestyle='--', alpha=0.7)
# ax_freq.set_ylim(bottom=0)
# ax_freq.set_xlim(Uvec[0], 5)

# plt.tight_layout()

################################################################################
# 3 - Flutter diagram parameter sweep
###############################################################################

# # =============================================================================
# # a ) SWEEP FLAP-NES Linear Damping: c_beta
# # =============================================================================

# sweep_results = {}
# for cb in c_beta:
#     print(f"\n{'='*60}")
#     print(f"EXECUTING Linear Stability Analysis FOR c_beta = {cb}")
#     print(f"{'='*60}")
    
#     # Execute the function with the current gb value
#     outputs = NON_DIM_LINEAR_AE(mu, OMEGA1, ah, ch[2], xa, ra, Zz, Za, cb, OMEGA2[1], xbeta[2], rbeta[2], k_h, k_a, GAMMA_beta,Uf, 100, nmodes,nstates)
    
#     # Store the returned tuple in the dictionary
#     sweep_results[cb] = outputs

# PLOT_Linear_Sweep(
#     sweep_results=sweep_results, 
#     param_label=r'$c_\beta$',  
#     LABEL_U_NONDIM=LABEL_U_NONDIM,
#     nmodes=nmodes
# )

# # =============================================================================
# # b) SWEEP FLAP-NES Linear Stiffness: OMEGA2
# # =============================================================================

# sweep_results = {}
# for OMEGA2_i in OMEGA2:
#    print(f"\n{'='*60}")
#    print(f"EXECUTING CONTINUATION FOR OMEGA2 = {OMEGA2_i}")
#    print(f"{'='*60}")
    
# #    Execute the function with the current gb value
#    outputs = NON_DIM_LINEAR_AE(mu, OMEGA1, ah, ch[0], xa, ra, Zz, Za, c_beta[0], OMEGA2_i, xbeta[0], rbeta[0], k_h, k_a, GAMMA_beta[1], Uf, 100, nmodes,nstates)
    
#     # Store the returned tuple in the dictionary
#    sweep_results[OMEGA2_i] = outputs

# PLOT_Linear_Sweep(
#    sweep_results=sweep_results, 
#    param_label=r'$\Omega_2$', 
#    LABEL_U_NONDIM=LABEL_U_NONDIM,
#    nmodes=nmodes
# )

# # =============================================================================
# # c ) SWEEP FLAP-NES Linear Stiffness: Size of the Flap
# # =============================================================================

# sweep_results = {}
# for i, ch_val in enumerate(ch):
#    print(f"\n{'='*60}")
#    print(f"EXECUTING Linear Stability Analysis FOR ch = {ch_val}")
#    print(f"{'='*60}")
   

#    outputs = NON_DIM_LINEAR_AE(mu, OMEGA1, ah, ch_val, xa, ra, Zz, Za, 
#                               c_beta[0], OMEGA2[1], xbeta[i], rbeta[i], 
#                                k_h, k_a, GAMMA_beta[0], Uf, 100, nmodes, nstates)
    
#    sweep_results[ch_val] = outputs

# PLOT_Linear_Sweep(
#    sweep_results=sweep_results, 
#    param_label=r'$c_h$', 
#    LABEL_U_NONDIM=LABEL_U_NONDIM,
#    nmodes=nmodes
# )


# Do not comment
plt.show()