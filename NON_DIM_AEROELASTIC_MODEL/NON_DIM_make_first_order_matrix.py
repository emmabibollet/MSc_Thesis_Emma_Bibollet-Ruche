import numpy as np

#Non-dimensional Aeroelastic model: Build the first order matrix, Q

def make_first_order_matrix(mu, C, D, E , F , W, W1, W2, M_inv, U, nmodes):

    l1=np.hstack((-np.dot(M_inv,C/U+D/mu),-np.dot(M_inv,E/(U**2)+F/mu), -np.dot(M_inv,W/(mu))))
    l2=np.hstack((np.eye(nmodes),np.zeros((nmodes,3*nmodes))))
    l3=np.hstack((np.zeros((2*nmodes,nmodes)),W1, W2))
    Q=np.vstack((l1,l2,l3))
    
    return Q