# Auxiliary functions for analysis of
#   ll_4320 mitgcm simulation
# crocha, sio summer 2014

import numpy as np

def rmean(A):
    """ Removes time-mean of llc_4320 3d fields; axis=2 is time"""
    ix,jx,kx = A.shape    
    Am = np.repeat(A.mean(axis=2),kx)
    Am = Am.reshape(ix,jx,kx)
    return A-Am

def spec_est_meridional(U,dx):
    """ Computes 1d (meridional) spectral estimates of 3d llc_4320 fields"""
    ix,jx,kx = U.shape
    N = ix          # record length
    df = 1./(N*dx)  # frequency resolution [cycles / (unit dx)]
    fNy = 1./(2*dx) # Nyquist frequency
    an = np.fft.fft(U,axis=0)

    an = an[1:N/2-1,:,:]
    E = 2*(an*an.conj())/df/(N**2)  # spectral estimate
    f = np.arange(1,N/2-1)*df

    return E.mean(axis=2),f,df,fNy 

def spec_est_zonal(U,dx):
    """ Computes 1d (zonal) spectral estimates of 3d llc_4320 fields"""
    ix,jx,kx = U.shape
    N = jx          # record length
    df = 1./(N*dx)  # frequency resolution [cycles / (unit dx)]
    fNy = 1./(2*dx) # Nyquist frequency
    an = np.fft.fft(U,axis=1)

    an = an[:,1:N/2-1,:]
    E = 2*(an*an.conj())/df/(N**2)  # spectral estimate
    f = np.arange(1,N/2-1)*df

    return E.mean(axis=2),f,df,fNy 

def spec_est_time(U,dt):
    """ Computes spectral estimate in time (axis=2) """
    ix,jx,kx = U.shape
    N = kx          # record length
    df = 1./(N*dt)  # frequency resolution [cycles / (unit time)]
    fNy = 1./(2*dt) # Nyquist frequency
    an = np.fft.fft(U,axis=2)

    an = an[:,:,1:N/2-1]
    E = 2*(an*an.conj())/df/(N**2)  # spectral estimate
    f = np.arange(1,N/2-1)*df
    return E.mean(axis=0),f,df,fNy 

def spec_error(E,sn,ci):

    """ Computes confidence interval for spectral 
        estimate E.
           sn is the number of spectral realizations (dof/2)
           ci = .95 for 95 % confidence interval 
        returns lower (El) and upper (Eu) bounds on E
        as well as pdf and cdf used to estimate errors """

    ## params
    dbin = .001
    yN = np.arange(0,5.+dbin,dbin)
    dof = 2.*sn # DOF = 2 x # of spectral estimates

    ## PDF for E/E0, where E (E0) is the estimate (true) 
    ##  process spectrum (basically a chi^2 distribution)
    C = dof / ( (2**sn) * np.math.gamma(sn) )  # constant
    pdf_yN = C * ( (dof*yN)**(sn-1) ) * np.exp( -(sn*yN) ) # chi^2(E/E0)

    ## CDF
    cdf_yN = np.cumsum(pdf_yN*dbin)  # trapezoidal-like integration

    ## compute confidence limits

    # lower
    el = ci
    fl = np.where( np.abs(cdf_yN - el) ==  np.abs(cdf_yN - el).min())
    El = E/yN[fl]

    # upper 
    eu = 1 - ci
    fu = np.where( np.abs(cdf_yN - eu) ==  np.abs(cdf_yN - eu).min())
    Eu = E/yN[fu]

    return El, Eu, cdf_yN, pdf_yN

# if sn larger than 150, assume it is normally-distributed (e.g., Bendat and Piersol) 
def spec_error2(E,sn):
    std_E = (1/np.sqrt(sn))
    El = E/(1 + 2*std_E)
    Eu = E/(1 - 2*std_E)
    return El, Eu, std_E

def spectral_slope(k,E,kmin,kmax,stdE):
    ''' compute spectral slope in log space in
        a wavenumber subrange [kmin,kmax],
        m: spectral slope; mm: uncertainty'''

    fr = np.where((k>=kmin)&(k<=kmax))

    ki = np.matrix((np.log10(k[fr]))).T
    Ei = np.matrix(np.log10(np.real(E[fr]))).T
    dd = np.matrix(np.eye(ki.size)*((np.abs(np.log10(stdE)))**2))

    G = np.matrix(np.append(np.ones((ki.size,1)),ki,axis=1))
    Gg = ((G.T*G).I)*G.T
    m = Gg*Ei
    mm = np.sqrt(np.array(Gg*dd*Gg.T)[1,1])
    yfit = np.array(G*m)
    m = np.array(m)[1]

    return m, mm

def leg_width(lg,fs):
    """"  Sets the linewidth of each legend object """
    for legobj in lg.legendHandles:
        legobj.set_linewidth(fs)

def auto_corr(x):
    """ Computes auto-correlation of 1d array """
    a = np.correlate(x,x,mode='full')
    a = a[x.size-1:]
    a = a/a[0]
    return a

def fit_gauss(x,y):
    """ Estimate characteristic scale of a auto-correlation function
            by fitting a Gaussian to auto_corr""" 
    y = np.matrix(np.log(y)).T
    A1 =  np .matrix(np.ones((x.size,1)))
    A2 = np.matrix(-(x**2)).T
    A = A2

    xmax = 650
    we = np.float(xmax) - x
    we =  (we/(xmax))**2
    We = np.matrix(np.diag(we))

    Gg = ((A.T*We*A).I)*A.T
    c = Gg*y

    Lfit = np.sqrt(1/c[-1])

    yfit = np.exp( A*c )

    return Lfit

def block_ave(dist,U,dx):
    """ Block-averages the array u onto grid with resolution dx  """
    ix,jx,kx = U.shape
    disti = np.arange(dx/2.,dist[-1]+dx/2.,dx)
    Ui = np.zeros((disti.size,jx,kx))
    for i in range(0,disti.size):
        fn = ((dist >= disti[i]-dx/2.) & (dist <= disti[i]+dx/2.))
        fns = np.sum(fn)
        if fns>0:
            Ui[i,:,:] = np.nansum(U[fn,:,:],axis=0)/fns
        else:
            Ui[i,:,:] = np.nan
    return Ui

