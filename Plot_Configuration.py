# -----------------------------------------------------------------------------------
# GLOBAL AEROELASTIC PLOT CONFIGURATION
# -----------------------------------------------------------------------------------
AERO_STYLE_NON_DIM = {
    'heave': {'label': 'Heave ($\\xi$)', 'color': 'tab:blue', 'unit': '[-]'},
    'pitch': {'label': 'Pitch ($\\alpha$)', 'color': 'tab:orange', 'unit': '[°]'},
    'flap':  {'label': 'Flap ($\\beta$)',  'color': 'tab:green', 'unit': '[°]'}
}

AERO_STYLE = {
    'heave': {'label': 'Heave ($h$)', 'color': 'tab:blue', 'unit': '[m]'},
    'pitch': {'label': 'Pitch ($\\alpha$)', 'color': 'tab:orange', 'unit': '[°]'},
    'flap':  {'label': 'Flap ($\\beta$)',  'color': 'tab:green', 'unit': '[°]'}
}

# Dimensional: Standardized Axis Labels
LABEL_TIME = 'Time [s]'
LABEL_U_DIM = '$U$ [m/s]'
LABEL_AMP_HEAVE_DIM = '$h$ Amplitude [m]'
LABEL_FREQ_DIM = '$\\omega$ [Hz]'
LABEL_AMP_ANGLE = '$\\alpha$ / $\\beta$ Amplitude [°]'

# NON-Dimensional: Standardized Axis Labels
LABEL_U_NONDIM = '$U^*$'
LABEL_U_NONDIM_NORM = '$U^*/U_L^*$'
LABEL_AMP_HEAVE_NONDIM = '$\\xi$ Amplitude [-]'
LABEL_TIME_NONDIM = '$\\tau$ [-]'
LABEL_FREQ_NONDIM = '$\\omega / \\omega_\\alpha$ [-]'