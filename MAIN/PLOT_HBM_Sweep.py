# =============================================================================
# Plotting Function: Plot Post-Futter Reponse for each set of parameters
# =============================================================================

# Imports ----------------------------------------------------------------------
import numpy as np
from matplotlib import pyplot as plt


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



def PLOT_HBM_Sweep(sweep_results, param_label, AERO_STYLE_NON_DIM, LABEL_U_NONDIM_NORM, LABEL_AMP_HEAVE_NONDIM, LABEL_AMP_ANGLE):
    fig1, ax1 = plt.subplots(figsize=(10, 6))
    ax2_twin = ax1.twinx()  # Create secondary Y-axis
    
    # Define distinct line styles for the sweeps
    line_styles = ['-', '--', ':', '-.']
    
    # Store proxy artists for the custom legend
    custom_lines = []
    custom_labels = []

    # Plot 3-DOF PARAMETER SWEEP -----------------------------------------------------------
    # Iterate through the stored results
    for idx, (param_val, data) in enumerate(sweep_results.items()):
        # Unpack the specific data
        U_plot_t, heave_amp, pitch_amp, flap_amp, stabi_arr, U_NS, heave_NS, pitch_NS, flap_NS, U_LP, heave_LP, pitch_LP, flap_LP = data

        U_plot = np.array(U_plot_t)
        U_NS = np.array(U_NS)if len(U_NS) > 0 else []
        U_LP = np.array(U_LP) if len(U_LP) > 0 else []

        # Normalize the airspeeds wrt the linear flutter speed (Hopf Bifurcation)
        #U_f = U_plot_t[0] #Flutter Speed
        #U_plot = np.array(U_plot_t) / U_f
        #U_NS = np.array(U_NS) / U_f if len(U_NS) > 0 else []
        #U_LP = np.array(U_LP) / U_f if len(U_LP) > 0 else []
        
        ls = line_styles[idx % len(line_styles)] # Assign line style for this sweep
        
        # Plot Stable Branches (Continuous Base)
        ax1.plot(U_plot, heave_amp, 
                 color=AERO_STYLE_NON_DIM['heave']['color'], 
                 linewidth=2, linestyle=ls)
    
        ax2_twin.plot(U_plot, pitch_amp, 
                      color=AERO_STYLE_NON_DIM['pitch']['color'], 
                      linewidth=2, linestyle=ls)
    
        ax2_twin.plot(U_plot, flap_amp, 
                      color=AERO_STYLE_NON_DIM['flap']['color'], 
                      linewidth=2, linestyle=ls)
    
        # Plot Stability: Red Overlay for Unstable
        unstable_heave = np.where(stabi_arr == 0, heave_amp, np.nan)
        unstable_pitch = np.where(stabi_arr == 0, pitch_amp, np.nan)
        unstable_flap  = np.where(stabi_arr == 0, flap_amp, np.nan)
    
        ax1.plot(U_plot, unstable_heave, color='red', linewidth=2)
        ax2_twin.plot(U_plot, unstable_pitch, color='red', linewidth=2)
        ax2_twin.plot(U_plot, unstable_flap, color='red', linewidth=2)
        
        # Plot Bifurcations
        if len(U_NS) > 0: # Neimark-Sacker Bifurcation
            ax1.plot(U_NS, heave_NS, color='purple', marker='s', linestyle='None', markersize=6)
            ax2_twin.plot(U_NS, pitch_NS, color='purple', marker='s', linestyle='None', markersize=6)
            ax2_twin.plot(U_NS, flap_NS, color='purple', marker='s', linestyle='None', markersize=6)
    
        if len(U_LP) > 0: # Limit Point Bifurcation
            ax1.plot(U_LP, heave_LP, color='black', marker='s', linestyle='None', markersize=6)
            ax2_twin.plot(U_LP, pitch_LP, color='black', marker='s', linestyle='None', markersize=6)
            ax2_twin.plot(U_LP, flap_LP, color='black', marker='s', linestyle='None', markersize=6)
    
        # Add a black line to the legend proxy to indicate which line style belongs to which parameter value
        custom_lines.append(plt.Line2D([0], [0], color='black', linewidth=2, linestyle=ls))
        custom_labels.append(f'{param_label} = {param_val}')
    
    # Formatting and Labels
    ax1.grid(True, linestyle='--', alpha=0.7)
    ax1.set_title("LCO Amplitude vs Airspeed")
    ax1.set_xlabel(LABEL_U_NONDIM_NORM) 
    ax1.set_ylabel(LABEL_AMP_HEAVE_NONDIM)
    ax2_twin.set_ylabel(LABEL_AMP_ANGLE)
    ax1.set_ylim(bottom=0)
    ax2_twin.set_ylim(bottom=0)
    
    # LEGEND ------------------------------------------------------------------------------------------
    custom_lines.extend([
        plt.Line2D([0], [0], color=AERO_STYLE_NON_DIM['heave']['color'], lw=2),
        plt.Line2D([0], [0], color=AERO_STYLE_NON_DIM['pitch']['color'], lw=2),
        plt.Line2D([0], [0], color=AERO_STYLE_NON_DIM['flap']['color'], lw=2),
        plt.Line2D([0], [0], color='red', lw=2),
        plt.Line2D([0], [0], color='purple', marker='s', linestyle='None', markersize=6),
        plt.Line2D([0], [0], color='black', marker='s', linestyle='None', markersize=6)
    ])
    custom_labels.extend([
        AERO_STYLE_NON_DIM['heave']['label'],
        AERO_STYLE_NON_DIM['pitch']['label'],
        AERO_STYLE_NON_DIM['flap']['label'],
        'Unstable',
        'Neimark-Sacker (NS)',
        'Limit Point (LP)'
    ])
    
    ax1.legend(custom_lines, custom_labels, loc='best')

    plt.tight_layout() 
    
