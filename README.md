This GitHub contains the Python code used for my MSc Thesis titled "Wing Flutter Control Using Passive Devices". The final report is included. For any questions, please contact me at emma@bibollet.fr

MAIN: This is the one you want to open :)
=> Linear stability analysis execution file: MAIN_LINEAR_STABILITY.py
=> Post-Flutter analysis execution file: MAIN_NON_LINEAR_POST_FLUTTER.py
=> System parameters: system_parameters.py
=> Harmonic balance parameters: HB_Parameters.py
=> Plotting functions for parametric sweep: PLOT_Linear_Sweep.py & PLOT_HBM_Sweep.py

The following folders contain the functions necessary for running the linear and post-flutter analysis
=> NON_DIM_AEROELASTIC MODEL: building the aeroelastic matrices from the system parameters
=> Non_DIM_HB_AE:
a) NON_DIM_HB_AE_Newton_Raphson.py = correcting the first point on the bifurcation branch using Newton-Raphson
b) NON_DIM_HB_Newton_Raphson_Cont.py = correcting the points along the bifurcation branch during the pseudo-arc length continuation scheme
c) NON_DIM_HB_AE_Ru.py & NON_DIM_HB_AE_Rx.py = calculating the residual derivatives with respect to airspeed (U) and displacements (X)
d) NON_DIM_HB_FNL.py = calculating the Harmonic Coefficients of the non-linear force
e) NON_DIM_Hopf_Bifurcation.py = locating the Hopf bifurcation airspeed and frequency
=> NON_DIM_Linear Analysis: contains the flutter_calc_test.py necessary to locate the Hopf bifurcation

To perform time-integration analysis, open the NON_DIM_Time_Integration folder
=> Bifurcation Diagram.py = obtain the complete bifurcation diagram. This requires precise initial conditions (refer to Appendix B of the report for more information)
=> Time Response.py = observe the motion of the system from an initial condition

Finally, the file Plot_Configurations.py centralizes the plot axes, legends, etc.
