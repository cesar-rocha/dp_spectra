# Python functions for computing pressure and density vertical modes
#   and vertical structure of SQG solutions
# Cesar B Rocha
# SIO, Summer 2014

import numpy as np
import scipy as sp
import scipy.linalg

def Dm(n):
    """ it creates forward difference matrix
            n: number of rows/columns """

    a = [1.,0.]
    b = sp.zeros(n-2)
    row1 = np.concatenate([a,b])
    a = [1.,-1.]
    col1 = np.concatenate([a,b])

    D = sp.linalg.toeplitz(col1,row1)

    return D

# four special matrices
def kctb(n):
    """ it creates four 'special' matrices for
            second order 1d problem.
            n: number of rows/columns """

    a = [2.,-1.]
    b = sp.zeros(n-2)
    row1 = np.concatenate([a,b])

    # Kn: Stiffness matrix (second order difference)
    K = sp.linalg.toeplitz(row1)

    # Cn: Circulant  matrix 
    C = sp.copy(K); C[0,n-1] = -1; C[n-1,0] = -1

    # Tn: Kn with changed the upper bc
    T = sp.copy(K); T[0,0] = 1

    # Bn: Tn with changed the lower bc
    B = sp.copy(T); B[n-1,n-1] = 1

    return K, C, T, B

# orthonormal eigenvector
def normal_evector(E,z):
    """ it normalizes eigenvectors such that 
           sum(Ei^2 x dz/H) = 1    """
    ix,jx =  E.shape
    dz = np.float(np.abs(z[1]-z[0]))
    H = dz*ix
    for j in range(jx):
        s = np.sqrt( (E[:,j]**2).sum()*(dz/H) )
        E[:,j] = E[:,j]/s
    return E

# compute vertical pressure modes
def pmodes(N2,z,lat,nm):

    ''' Compute vertical modes '''

    # Settings
    dz = np.float(np.abs(z[1])-np.abs(z[0]))
    f2 = (2.*(7.29e-5)*np.sin(lat*(np.pi/180.)))**2

    # Assembling matrices
    C = np.matrix( np.diag(f2/N2) )
    D = Dm(N2.size)
    
    K = (D.T*C*D)/(dz**2)
    
    # Enforce upper BC
    K[0,0] = -K[0,1]

    w,v = np.linalg.eigh(K)

    v = v[:,w.argsort()]
    w = np.sort(w)
    w = np.array(w[:nm])
    v = np.array(v[:,:nm])
    
    # Here the normalization const.
    #   is ALWAYS H.
    v = normal_evector(v,z)

    v[:,0] = v[:,0]*0. + 1. # set barotropic = 1 for better accuracy
    w[0] = 0.  # zeroth eigenvalue is always zero 
               #    within machine precision

    return v,w

# vertical structure of SQG solution
def sqgz(N2,z,lat,k,norm='False'):

    """ Compute the SQG vertical structure 
        numerically by solving the BVP
       d/dz( f2/N2 * d/dz )F - k^2 F = 0
       subject to dF/dz  = N2/f2 @ z = 0
       and dF/dz = 0 @ z = -H
     N2 = stratification squared [(cps)^2]
     lat = local latitude
     k = wavenumber [cpm] """

    # Settings
    dz = np.float(np.abs(z[1])-np.abs(z[0]))
    f2 = (2.*(7.29e-5)*np.sin(lat*(np.pi/180.)))**2

    # Assembling matrices
    C = np.matrix(np.diag(f2/N2) )
    D = Dm(N2.size)

    M1 = -(D.T*C*D)
    M1[0,0],M1[0,1] = dz,-dz   # Enforce boundary conditions
    M1[-1,-2],M1[-1,-1] = dz,-dz

    M1 = M1/(dz**2)

    M2 = np.matrix(np.eye(N2.size))*(k**2)
    M2[0,0],M2[-1,-1] = 0.,0.

    # Point load (surface buoyancy)
    f = np.zeros(N2.size)
    f[0] = N2[0]/f2

    v = sp.linalg.solve(M1-M2,f)

    # normalize to such that v[0] = 1
    if norm == True:
        v = v/v[0]

    return v        

