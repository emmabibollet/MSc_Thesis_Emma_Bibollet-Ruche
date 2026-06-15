# =============================================================================
# Plotting Function: Plot Flutter Diagrams for each set of Flap-NES parameters
# =============================================================================

# Imports ----------------------------------------------------------------------
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.lines as mlines

plt.rcParams.update({
    'font.size': 16,          # Global text size
    'axes.titlesize': 16,     # Title size
    'axes.labelsize': 16,     # Axis label size
    'xtick.labelsize': 16,    # X-tick size
    'ytick.labelsize': 16,    # Y-tick size
    'legend.fontsize': 16,    # Legend text size
    'font.family': 'serif',   # Professional serif font (like LaTeX)
    'lines.markersize': 10    # Make markers visible
})    

def PLOT_Linear_Sweep(sweep_results, param_label, LABEL_U_NONDIM, nmodes):

    fig_damp, ax_damp = plt.subplots(figsize=(10, 6))
    fig_freq, ax_freq = plt.subplots(figsize=(10, 6))

    marker_colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray']
    marker_styles = ['o', 's', '^', 'D', 'v']

    custom_lines_sweep = []
    custom_labels_sweep = []
    flutter_text_lines = []

    # Iterate through the sweep results
    for idx, (param_val, data) in enumerate(sweep_results.items()):
        # Unpack the returned tuple from NON_DIM_LINEAR_AE
        Uvec, Zeta_matrix, Freq_matrix, Uflut, omega_flut = data
        
        ms = marker_styles[idx % len(marker_styles)]
        mc = marker_colors[idx % len(marker_colors)]
        
        if len(Uflut) > 0:
            # Format flutter speeds for the text box
            u_str = ", ".join([f"{u:.3f}" for u in Uflut if u > 0])
            flutter_text_lines.append(f"{param_label} = {param_val}: $U_f = {u_str}$")
            #flutter_text_lines.append(f"3-DOF: $U_f = {u_str}$")
        else:
            flutter_text_lines.append(f"{param_label} = {param_val}: Stable")

        for k in range(nmodes):
            # Plot Damping
            ax_damp.plot(Uvec, Zeta_matrix[k, :], color = mc , marker=ms, linestyle='None', markersize=4)
            # Plot Frequency
            ax_freq.plot(Uvec, Freq_matrix[k, :], color = mc, marker=ms, linestyle='None', markersize=4)

        # Mark Flutter crossings on the Damping plot (where zeta crosses 0)
        if len(Uflut) > 0:
            for u_val in np.atleast_1d(Uflut):
                if u_val > 0:
                    ax_damp.plot(u_val, 0, color='red', marker='x', markersize=8, markeredgewidth=2, zorder=5)

        # Record the line style for the custom legend
        custom_lines_sweep.append(mlines.Line2D([0], [0], color = mc, marker=ms, linestyle='None', markersize=4))
        custom_labels_sweep.append(f'{param_label} = {param_val}')


    # ---------------------------------------------------------
    # Render the Flutter Speed Text Box
    # ---------------------------------------------------------
    flutter_text_block = "\n".join(flutter_text_lines)
    box_props = dict(boxstyle='square,pad=0.5', facecolor='white', alpha=0.9, edgecolor='black')
    
    # Place text box in the upper right corner of both plots
    ax_damp.text(0.95, 0.95, flutter_text_block, transform=ax_damp.transAxes, 
                 fontsize=16, verticalalignment='top', horizontalalignment='right', bbox=box_props, zorder=10)
    ax_freq.text(0.95, 0.95, flutter_text_block, transform=ax_freq.transAxes, 
                 fontsize=16, verticalalignment='top', horizontalalignment='right', bbox=box_props, zorder=10)

    # ---------------------------------------------------------
    # Formatting Damping Plot
    # ---------------------------------------------------------
    ax_damp.axhline(0, color='black', linewidth=1)
    ax_damp.set_xlabel(LABEL_U_NONDIM)
    ax_damp.set_ylabel(r'Damping Ratio $\zeta$')
    ax_damp.grid(True, linestyle='--', alpha=0.7)
    ax_damp.set_ylim(-0.2, 0.5)
    ax_damp.set_xlim(Uvec[0], Uvec[-1])

    # ---------------------------------------------------------
    # Formatting Frequency Plot
    # ---------------------------------------------------------
    ax_freq.axhline(0, color='black', linewidth=1)
    ax_freq.set_xlabel(LABEL_U_NONDIM)
    ax_freq.set_ylabel(r'Frequency Ratio ($\omega / \omega_\alpha$)')
    ax_freq.grid(True, linestyle='--', alpha=0.7)
    ax_freq.set_ylim(bottom=0)
    ax_freq.set_xlim(Uvec[0], Uvec[-1])

    # ---------------------------------------------------------
    # Unified Legend Construction
    # ---------------------------------------------------------
    # Append sweep parameters
    legend_elements = custom_lines_sweep.copy()
    labels = custom_labels_sweep.copy()

    ax_damp.legend(legend_elements, labels, loc='lower left')
    ax_freq.legend(legend_elements, labels, loc='center left')

    fig_damp.tight_layout()
    fig_freq.tight_layout()
    plt.show()


def PLOT_Linear_2DOF(outputs, LABEL_U_NONDIM, nmodes):

    fig_damp, ax_damp = plt.subplots(figsize=(10, 6))
    fig_freq, ax_freq = plt.subplots(figsize=(10, 6))

    # Unpack the returned tuple directly
    Uvec, Zeta_matrix, Freq_matrix, Uflut, omega_flut = outputs
    
    # Establish a unified color for the baseline system
    mc = 'tab:blue'
    ms = 'o'

    flutter_text_lines = []

    # 1. Format flutter speeds for the text box
    if len(Uflut) > 0:
        u_str = ", ".join([f"{u:.3f}" for u in Uflut if u > 0])
        flutter_text_lines.append(f"2-DOF Baseline: $U_f = {u_str}$")
    else:
        flutter_text_lines.append("2-DOF Baseline: Stable")

    # 2. Iterate strictly up to nmodes (2) to avoid plotting empty zero-rows
    for k in range(nmodes):
        # Plot discrete points from the 1D matrix slices
        ax_damp.plot(Uvec, Zeta_matrix[k, :], color=mc, marker=ms, linestyle='None', markersize=4)
        ax_freq.plot(Uvec, Freq_matrix[k, :], color=mc, marker=ms, linestyle='None', markersize=4)

    # 3. Mark Flutter crossings on the Damping plot
    if len(Uflut) > 0:
        for u_val in np.atleast_1d(Uflut):
            if u_val > 0:
                ax_damp.plot(u_val, 0, color='red', marker='x', markersize=10, markeredgewidth=3, zorder=5)

    # ---------------------------------------------------------
    # Render the Flutter Speed Text Box
    # ---------------------------------------------------------
    flutter_text_block = "\n".join(flutter_text_lines)
    box_props = dict(boxstyle='square,pad=0.5', facecolor='white', alpha=0.9, edgecolor='black')
    
    ax_damp.text(0.95, 0.95, flutter_text_block, transform=ax_damp.transAxes, 
                 fontsize=16, verticalalignment='top', horizontalalignment='right', bbox=box_props, zorder=10)
    ax_freq.text(0.95, 0.95, flutter_text_block, transform=ax_freq.transAxes, 
                 fontsize=16, verticalalignment='top', horizontalalignment='right', bbox=box_props, zorder=10)

    # ---------------------------------------------------------
    # Formatting Damping Plot
    # ---------------------------------------------------------
    ax_damp.axhline(0, color='black', linewidth=1)
    ax_damp.set_xlabel(LABEL_U_NONDIM)
    ax_damp.set_ylabel(r'Damping Ratio $\zeta$')
    ax_damp.grid(True, linestyle='--', alpha=0.7)
    ax_damp.set_ylim(-0.2, 0.5)
    ax_damp.set_xlim(Uvec[0], Uvec[-1])

    # ---------------------------------------------------------
    # Formatting Frequency Plot
    # ---------------------------------------------------------
    ax_freq.axhline(0, color='black', linewidth=1)
    ax_freq.set_xlabel(LABEL_U_NONDIM)
    ax_freq.set_ylabel(r'Frequency Ratio ($\omega / \omega_\alpha$)')
    ax_freq.grid(True, linestyle='--', alpha=0.7)
    ax_freq.set_ylim(bottom=0)
    ax_freq.set_xlim(Uvec[0], Uvec[-1])

    # ---------------------------------------------------------
    # Legend Construction
    # ---------------------------------------------------------
    legend_line = mlines.Line2D([0], [0], color=mc, marker=ms, linestyle='None', markersize=8)
    
    ax_damp.legend([legend_line], ['2-DOF Eigenvalues'], loc='lower left')
    ax_freq.legend([legend_line], ['2-DOF Eigenvalues'], loc='center left')

    fig_damp.tight_layout()
    fig_freq.tight_layout()
    plt.show()