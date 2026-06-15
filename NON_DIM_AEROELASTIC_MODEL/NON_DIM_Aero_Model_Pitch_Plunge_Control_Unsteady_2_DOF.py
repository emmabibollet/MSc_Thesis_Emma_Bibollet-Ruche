# NON DIMENSIONAL AEROELASTIC MODEL: 2D PITCH-PLUNGE-CONTROL WING SECTION WITH UNSTEADY AERODYNAMICS (JESUS PhD: o)

#######################################################################
# Imports
#######################################################################
import numpy as np

#######################################################################
# Model
#######################################################################

class Model_NON_DIM: 

    ########################################################################
    # Initialize Parameters
    ########################################################################

    # Define the list of critical parameters you expect to receive
    CRITICAL_PARAMS = ['mu', 'OMEGA1', 'ah', 'xa','Zz', 'Za', 'k_h', 'k_a']

    def __init__(self, **kwargs):

        # DEFAULT VALUES -------------------------------------------------------------
        self.mu = 100
        self.OMEGA1 = 0.5
        self.ah = -0.5
        self.xa = 0.25
        self.ra = 0.5
        self.Zz = 0
        self.Za = 0.02
        self.k_h = 0
        self.k_a = 0

        # OVERWRITE WITH USER AGRUMENTS (from MAIN_AEROELASTIC CALCULATION.py) --------
        for key, value in kwargs.items():
            setattr(self, key, value)

        # VERIFICATION: Inform if there are no user arguments --------------------------
        if not kwargs: # 0 user arguments, kwargs is empty
            print("WARNING: DEFAULT VALUES FOR PARAMETERS USED") 

        # Find which parameters are missing from user arguments
        missing_params = [ param for param in self.CRITICAL_PARAMS 
            if param not in kwargs ]
        if missing_params:
            print(f"WARNING: The following CRITICAL parameters were missing from input and are using DEFAULT values: {', '.join(missing_params)}")

        
        # EXPONENTIAL APPROXIMATION PARAMETERS OF WAGNER FUNCTION ----------------------
        self.psi1 = 0.165
        self.psi2 = 0.335  
        self.eps1= 0.0455
        self.eps2= 0.3

        #WAGNER FUNCTION AT t=0
        self.phi0 = 1-self.psi1 - self.psi2

        #WAGNER DERIVATIVE AT t = 0
        self.phi0dot = self.psi1*self.eps1+self.psi2*self.eps2


    ########################################################################
    # STRUCTURAL Mass Matrix (A) : Intertial Coupling
    ######################################################################## 
    def make_mass_matrix(self):

        A=np.array([[1 ,self.xa],
                   [self.xa/(self.ra**2),1]])
        return A
    
    ########################################################################
    # STRUCTURAL Stiffness Matrix (E) : Uncoupled
    ########################################################################
    def make_stiffness_matrix(self):
        E=np.array([[self.OMEGA1**2,0],
                   [0,1]])
        
        return E 
    
    ########################################################################
    # STRUCTURAL NON_LINEAR Stiffness Matrix (K_nl) : Uncoupled
    ########################################################################
    def make_CUBIC_stiffness_matrix(self):
        K_nl=np.array([[self.OMEGA1**2*self.k_h,0],
                   [0,1*self.k_a]])
        
        return K_nl 

    
    ########################################################################
    # STRUCTURAL Damping Matrix (C) : Uncoupled
    ########################################################################
    def make_damping_matrix(self):
        C=np.array([[2*self.Zz*self.OMEGA1,0],
                   [0,2*self.Za]])
        return C

    ########################################################################
    # AERO MASS MATRIX (B) : Non-Circulatory, No lag, acceleration terms
    ########################################################################
    def make_aero_mass_matrix(self):
        B=np.array([[1, -self.ah],
                   [-self.ah/(self.ra**2), (1/8+self.ah**2)/(self.ra**2)]])
        return B
    
    ########################################################################
    # AERO DAMPING MATRIX (D): velocity terms
    ########################################################################
    def make_aerodynamic_damping_matrix(self):
        #Non Circulatory
        D1=np.array([[0, 1],
                   [0, (1/2-self.ah)/(self.ra**2)]])
        
        #Circulatory
        D2=np.array([[2, 2*(1/2-self.ah)],
                   [-2*(self.ah+1/2)/(self.ra**2), -2*(self.ah+1/2)*(1/2-self.ah)/(self.ra**2)]])

        return D1 + self.phi0*D2

    ########################################################################
    # AERO STIFFNESS MATRIX (F)
    ########################################################################
    def make_aerodynamic_stiffness_matrix(self):
        #Non Circulatory
        F1=np.array([[0,0],
                   [0,0]])
        
        #Circulatory (factor wagner function at t= 0)
        F2=np.array([[0,2],
                   [0,-2*(self.ah+1/2)/(self.ra**2)]])
        
        #Circulatory (factor of derivative of wagner function at t = 0)
        F3=np.array([[2,2*(1/2-self.ah)],
                   [-2*(self.ah+1/2)/(self.ra**2), -2*(self.ah+1/2)*(1/2-self.ah)/(self.ra**2)]])
        
        return F1+self.phi0*F2+self.phi0dot*F3

    ########################################################################
    # AERO STATE INFLUENCE MATRIX (W)
    ########################################################################
    def make_aerodynamic_influence_matrix(self):
        W0=np.array([[-self.psi1*(self.eps1)**2],
                   [-self.psi2*(self.eps2)**2],
                   [self.psi1*self.eps1*(1-self.eps1*(1/2-self.ah))],
                   [self.psi2*self.eps2*(1-self.eps2*(1/2-self.ah))]])
        
        W=np.hstack((2*W0, -2*(self.ah+1/2)*W0/self.ra**2)).T
        return W
    
    ########################################################################
    # AERO STATE EQUATION MATRICES (W1 and W2)
    ########################################################################
    def make_aerodynamic_state_equation_matrices(self):
        W1=np.array([[1,0],
                     [1,0],
                     [0,1],
                     [0,1]])
        
        W2=np.array([[-self.eps1,0,0,0],
                     [0,-self.eps2,0,0],
                     [0,0,-self.eps1,0],
                     [0,0,0,-self.eps2]])
        return W1,W2

    ########################################################################
    # INITIAL CONDITION EXCITATION VECTOR (g)
    ########################################################################
    def compute_initial_excitation(self, X0):
        g=(X0[2]+(1/2-self.ah)*X0[3])*np.array([[2],
                                                [-2*(self.ah+1/2)/(self.ra**2)]])
    
        return g
    
    ########################################################################
    # WAGNER FUNCTION DERIVATIVE (t is non dimensional time)
    ########################################################################
    def Wagner_function_derivative(self,t):
        return self.psi1 * self.eps1 * np.exp(-self.eps1 * t) + self.psi2 * self.eps2 * np.exp(-self.eps2 * t)
    

    ########################################################################
    # INITIAL EXCITATION
    ########################################################################    
    def make_initial_excitation(self, X0):
    
        A=self.make_mass_matrix()
        B=self.make_aero_mass_matrix()
        M=A+B/self.mu
        M_inv=np.linalg.inv(M)
        g=self.compute_initial_excitation(X0)

        q = np.vstack((np.dot(M_inv, g)/self.mu, np.zeros((3, 1)), np.zeros((3, 1)), np.zeros((3, 1)))).reshape(12)
        return lambda t: q * self.Wagner_function_derivative(t)


   
    
