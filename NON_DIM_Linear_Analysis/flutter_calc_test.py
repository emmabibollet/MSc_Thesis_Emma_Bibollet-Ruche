import numpy as np

def flutter_calc_test (Q):

    eigis = np.linalg.eigvals(Q)
    
    # Determine the complex eigenvalues (masking)
    comp_eigs = eigis[np.abs(np.imag(eigis)) > 1e-6]

    # Sort the complex eigenvalues by descending imaginary part
    iko = np.argsort(np.abs(np.imag(comp_eigs)))[::-1]
    comp_eigs = comp_eigs[iko]

    # Detects Hopf bifurcation by monitoring the real parts of complex eigenvalue pairs (λ = σ ± iω). 
    # Since (σ + iω) + (σ - iω) = 2σ, the product 'tau_flut' will change sign when the leading real part (σ) crosses from stable (negative) to unstable (positive).
    tau_flut = np.real(np.prod(comp_eigs[0::2] + comp_eigs[1::2]) )

    # Pitchfork (static divergence) test function
    real_eigs = eigis[np.abs(np.imag(eigis)) < 1e-4]
    tau_div = np.prod(np.real(real_eigs))

    return tau_flut, tau_div, comp_eigs, eigis