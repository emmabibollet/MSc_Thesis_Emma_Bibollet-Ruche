# =============================================================================
# HARMONIC BALANCE PARAMETERS: NON-DIMENSIONAL
# =============================================================================

# Included Harmonics
n_hbm =5 #number of harmonics included
nfft = 2**8 #number of FFT time points to determine the harmonic coefficients of Fnl

#Newton-Raphson ---------------------------------------------------------------
eps_NR = 1e-9 #Newton-Raphson precision
iter_NR = 30 # Maximum number of iterations for Newton Raphson scheme

#Continuation -----------------------------------------------------------------
ds = 1e-4 #initial length of the step ds
dir_cont = -1 #direction of the continuation scheme (either 1 or -1)
avg_iter = 3 # to help adapt the size of ds
ds_min = 1e-4
ds_max = 0.1

Uf =8 # Maximum airspeed
max_steps = 4 #Maximum number of steps for continuation

#Initial Conditions ----------------------------------------------------------
dU = 0.0001 #slightly higher or lover velocity
dX = 0.01 #give a little initial ampltiude to kick off of the stable branch
DOF_i =5 # DOF to which the initial amplitude is given (h = 3, alpha =4 and beta = 5)

#Stability -------------------------------------------------------------------
precStab = 1e-6