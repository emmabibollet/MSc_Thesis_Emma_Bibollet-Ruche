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
    CRITICAL_PARAMS = ['mu', 'OMEGA1', 'ah', 'ch', 'xa', 'ra', 'Zz', 'Za', 'c_beta', 'OMEGA2', 'xbeta', 'rbeta', 'k_h', 'k_a', 'GAMMA_beta']

    def __init__(self, **kwargs):

        # DEFAULT VALUES -------------------------------------------------------------
        self.mu = 100
        self.OMEGA1 = 0.5
        self.ah = -0.5
        self.ch = 0.5
        self.xa = 0.25
        self.ra = 0.5
        self.Zz = 0
        self.Za = 0.02
        self.OMEGA2 = 3.5
        self.xbeta = 0.015
        self.rbeta = 0.08
        self.c_beta = 1.2
        self.k_h = 0
        self.k_a = 0
        self.GAMMA_beta = 100

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

        # CALCULATE DERIVED CONSTANTS --------------------------------------------------
        self.nu = np.arccos(self.ch)
        
        # EXPONENTIAL APPROXIMATION PARAMETERS OF WAGNER FUNCTION ----------------------
        self.psi1 = 0.165
        self.psi2 = 0.335  
        self.eps1= 0.0455
        self.eps2= 0.3

        # THEODORSEN COEFFICIENTS (T1 to T14)--------------------------------------------
        self.make_T_variables()

        #WAGNER FUNCTION AT t=0
        self.phi0 = 1-self.psi1 - self.psi2

        #WAGNER DERIVATIVE AT t = 0
        self.phi0dot = self.psi1*self.eps1+self.psi2*self.eps2


    ########################################################################
    # STRUCTURAL Mass Matrix (A) : Intertial Coupling
    ######################################################################## 
    def make_mass_matrix(self):

        A=np.array([[1 ,self.xa ,self.xbeta],
                   [self.xa/(self.ra**2),1, (self.rbeta**2+(self.ch-self.ah)*self.xbeta)/(self.ra**2)],
                   [self.xbeta/(self.rbeta**2),(self.rbeta**2+(self.ch-self.ah)*self.xbeta)/(self.rbeta**2), 1]])
        return A
    
    ########################################################################
    # STRUCTURAL Stiffness Matrix (E) : Uncoupled
    ########################################################################
    def make_stiffness_matrix(self):
        E=np.array([[self.OMEGA1**2,0,0],
                   [0,1,0],
                  [0,0,self.OMEGA2**2]]) #With linear stiffness
        
        return E 
    
    ########################################################################
    # STRUCTURAL NON_LINEAR Stiffness Matrix (K_nl) : Uncoupled
    ########################################################################
    def make_CUBIC_stiffness_matrix(self):
        K_nl=np.array([[self.OMEGA1**2*self.k_h,0,0],
                   [0,1*self.k_a,0],
                  [0,0, self.GAMMA_beta]]) #With linear stiffness
        
        return K_nl 

    
    ########################################################################
    # STRUCTURAL Damping Matrix (C) : Uncoupled
    ########################################################################
    def make_damping_matrix(self):
        C=np.array([[2*self.Zz*self.OMEGA1,0,0],
                   [0,2*self.Za,0],
                  [0,0,self.c_beta]])
        return C

    ########################################################################
    # AERO MASS MATRIX (B) : Non-Circulatory, No lag, acceleration terms
    ########################################################################
    def make_aero_mass_matrix(self):
        B=np.array([[1, -self.ah, -self.T1/np.pi],
                   [-self.ah/(self.ra**2), (1/8+self.ah**2)/(self.ra**2),-(self.T7+(self.ch-self.ah)*self.T1)/(np.pi*(self.ra**2))],
                   [-self.T1/(np.pi*self.rbeta**2),2*self.T13/(np.pi*self.rbeta**2),-self.T3/(np.pi*self.rbeta)**2]])
        return B
    
    ########################################################################
    # AERO DAMPING MATRIX (D): velocity terms
    ########################################################################
    def make_aerodynamic_damping_matrix(self):
        #Non Circulatory
        D1=np.array([[0, 1 ,-self.T4/np.pi],
                   [0, (1/2-self.ah)/(self.ra**2),(self.T1-self.T8-(self.ch-self.ah)*self.T4+self.T11/2)/(np.pi*(self.ra**2))],
                   [0,(-2*self.T9-self.T1-self.T4*(0.5-self.ah))/(np.pi*self.rbeta**2),-self.T4*self.T11/(2*np.pi**2*self.rbeta**2)]])
        
        #Circulatory
        D2=np.array([[2, 2*(1/2-self.ah),self.T11/np.pi],
                   [-2*(self.ah+1/2)/(self.ra**2),-2*(self.ah+1/2)*(1/2-self.ah)/(self.ra**2),-(self.ah+1/2)*self.T11/(np.pi*(self.ra**2))],
                   [self.T12/(np.pi*self.rbeta**2),self.T12/np.pi*(1/2-self.ah)/(self.rbeta**2),self.T12*self.T11/(2*np.pi**2*self.rbeta**2)]])

        return D1 + self.phi0*D2

    ########################################################################
    # AERO STIFFNESS MATRIX (F)
    ########################################################################
    def make_aerodynamic_stiffness_matrix(self):
        #Non Circulatory
        F1=np.array([[0,0,0],
                   [0,0,(self.T4+self.T10)/(np.pi*(self.ra**2))],
                   [0,0,(self.T5-self.T4*self.T10)/(np.pi**2*self.rbeta**2)]])
        
        #Circulatory (factor wagner function at t= 0)
        F2=np.array([[0,2,2*self.T10/np.pi],
                   [0,-2*(self.ah+1/2)/(self.ra**2),-2*(self.ah+1/2)*self.T10/(np.pi*(self.ra**2))],
                   [0, self.T12/(np.pi*self.rbeta**2) , self.T12*self.T10/(np.pi**2*self.rbeta**2)]])
        
        #Circulatory (factor of derivative of wagner function at t = 0)
        F3=np.array([[2,2*(1/2-self.ah), self.T11/np.pi],
                   [-2*(self.ah+1/2)/(self.ra**2),-2*(self.ah+1/2)*(1/2-self.ah)/(self.ra**2),-(self.ah+1/2)*self.T11/(np.pi*(self.ra**2))],
                   [self.T12/(np.pi*self.rbeta**2), self.T12/np.pi*(1/2-self.ah)/(self.rbeta**2), self.T12*self.T11/(2*np.pi**2*self.rbeta**2)]])
        
        return F1+self.phi0*F2+self.phi0dot*F3

    ########################################################################
    # AERO STATE INFLUENCE MATRIX (W)
    ########################################################################
    def make_aerodynamic_influence_matrix(self):
        W0=np.array([[-self.psi1*(self.eps1)**2],
                   [-self.psi2*(self.eps2)**2],
                   [self.psi1*self.eps1*(1-self.eps1*(1/2-self.ah))],
                   [self.psi2*self.eps2*(1-self.eps2*(1/2-self.ah))],
                   [self.psi1*self.eps1*(self.T10-self.eps1*self.T11/2)/(np.pi)],
                   [self.psi2*self.eps2*(self.T10-self.eps2*self.T11/2)/(np.pi)]])
        
        W=np.hstack((2*W0, -2*(self.ah+1/2)*W0/self.ra**2, self.T12*W0/(np.pi*self.rbeta**2))).T
        return W
    
    ########################################################################
    # AERO STATE EQUATION MATRICES (W1 and W2)
    ########################################################################
    def make_aerodynamic_state_equation_matrices(self):
        W1=np.array([[1,0,0],
                     [1,0,0],
                     [0,1,0],
                     [0,1,0],
                     [0,0,1],
                     [0,0,1]])
        W2=np.array([[-self.eps1,0,0,0,0,0],
                     [0,-self.eps2,0,0,0,0],
                     [0,0,-self.eps1,0,0,0],
                     [0,0,0,-self.eps2,0,0],
                     [0,0,0,0,-self.eps1,0],
                     [0,0,0,0,0,-self.eps2]])
        return W1,W2

    ########################################################################
    # INITIAL CONDITION EXCITATION VECTOR (g)
    ########################################################################
    def compute_initial_excitation(self, X0):
        g=(X0[3]+(1/2-self.ah)*X0[4]+self.T11/(2*np.pi)*X0[5])*np.array([[2],
                                                                        [-2*(self.ah+1/2)/(self.ra**2)],
                                                                        [self.T12/(np.pi*self.rbeta**2)]])
    
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

    ########################################################################
    # UNSTEADY AERODYNAMICS : Theodorsen Coefficients (T1 to T14)
    ########################################################################
    def make_T_variables(self):
        self.T1=-1/3*np.sqrt(1-self.ch**2)*(2+self.ch**2)+self.ch*self.nu
        self.T2=self.ch*(1-self.ch**2)-np.sqrt(1-self.ch**2)*(1+self.ch**2)*self.nu+self.ch*self.nu**2
        self.T3=-(1/8+self.ch**2)*self.nu**2+1/4*self.ch*np.sqrt(1-self.ch**2)*self.nu*(7+2*self.ch**2)-1/8*(1-self.ch**2)*(5*self.ch**2+4)
        self.T4=-self.nu+self.ch*np.sqrt(1-self.ch**2)
        self.T5=-(1-self.ch**2)-self.nu**2+2*self.ch*np.sqrt(1-self.ch**2)*self.nu
        self.T6=self.T2
        self.T7=-(1/8+self.ch**2)*self.nu+1/8*self.ch*np.sqrt(1-self.ch**2)*(7+2*self.ch**2)
        self.T8=-1/3*np.sqrt(1-self.ch**2)*(2*self.ch**2+1)+self.ch*self.nu
        self.T9=1/2*(1/3*(1-self.ch**2)**(3/2)+self.ah*self.T4)
        self.T10=np.sqrt(1-self.ch**2)+self.nu
        self.T11=self.nu*(1-2*self.ch)+np.sqrt(1-self.ch**2)*(2-self.ch)
        self.T12=np.sqrt(1-self.ch**2)*(2+self.ch)-self.nu*(2*self.ch+1)
        self.T13=1/2*(-self.T7-(self.ch-self.ah)*self.T1)
        self.T14=1/16+1/2*self.ah*self.ch

   
    
