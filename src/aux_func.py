import numpy as np
import scipy as sp

def spec_est2(A,d1,d2,win=True):

    """    computes 2D spectral estimate of A
           obs: the returned array is fftshifted
           and consistent with the f1,f2 arrays
           d1,d2 are the sampling rates in rows,columns   """
    
    import numpy as np

    l1,l2,l3 = A.shape
    df1 = 1./(l1*d1)
    df2 = 1./(l2*d2)
    f1Ny = 1./(2*d1)
    f2Ny = 1./(2*d2)

    f1 = np.arange(-f1Ny,f1Ny,df1)
    f2 = np.arange(-f2Ny,f2Ny,df2)
    
    if win == True:
        wx = np.matrix(np.hanning(l1))
        wy =  np.matrix(np.hanning(l2))
        window_s = np.repeat(np.array(wx.T*wy),l3).reshape(l1,l2,l3)
    else:
        window_s = np.ones((l1,l2,l3))

    an = np.fft.fft2(A*window_s,axes=(0,1))
    E = (an*an.conjugate()) / (df1*df2) / ((l1*l2)**2)
    E = np.fft.fftshift(E)
    E = E.mean(axis=2)

    return np.real(E),f1,f2,df1,df2,f1Ny,f2Ny

def ps(u,v,dx,dy):

    """ decompose the vector field (u,v) into potential (up,vp)
        and solenoidal (us,vs) fields using 2D FT a la Smith JPO 2008 """

    ix,jx,kx = u.shape
    dl = 1./(ix*dy)
    dk = 1./(jx*dx)
    kNy = 1./(2*dx)
    lNy = 1./(2*dy)
    k = np.arange(-kNy,kNy,dk)
    k = np.fft.fftshift(k)
    l = np.arange(-lNy,lNy,dl)
    l = np.fft.fftshift(l)
    K,L = np.meshgrid(k,l)
    THETA = (np.arctan2(L,K))
    THETA = np.repeat(THETA,kx).reshape(ix,jx,kx)

    U = np.fft.fft2(u,axes=(0,1))
    V = np.fft.fft2(v,axes=(0,1))

    P = U*np.cos(THETA) + V*np.sin(THETA)
    S = -U*np.sin(THETA) + V*np.cos(THETA)

    # back to physical space
    up = np.real(np.fft.ifft2(P*np.cos(THETA),axes=(0,1)))
    vp = np.real(np.fft.ifft2(P*np.sin(THETA),axes=(0,1)))

    us = np.real(np.fft.ifft2(-S*np.sin(THETA),axes=(0,1)))
    vs = np.real(np.fft.ifft2(S*np.cos(THETA),axes=(0,1)))

    return up,vp,us,vs

def bcf(k,EU,EV):
    """ Compute Helmholtz decomposition in wavenumber space
            a la Buhler et al JFM in press
             obs: notation different than BCF: Here u,v is across-,along-track"""

    ki=np.linspace(k[0],k[-1],1000)
    interpu,interpv=sp.interpolate.interp1d(k,EU),sp.interpolate.interp1d(k,EV)
    
    EU,EV,k=interpu(ki),interpv(ki),ki
    

    dk=k[1]-k[0]
    kp=np.arange(k[-1]+dk,k[-1]+4*dk,dk)
    Ep=1./(kp**2)
    Epu=(Ep/Ep[0])*EU[-1]
    Epv=(Ep/Ep[0])*EV[-1]
    k=np.append(k,kp)
    EU=np.append(EU,Epu)
    EV=np.append(EV,Epv)

    s=np.log(k)

    Dpsi,Dphi=np.zeros(s.size),np.zeros(s.size)

    for i in range(s.size-1):
        ds=np.diff(s[i:])
        Dpsi[i] = ((EV[i+1:]*np.sinh(s[i]-s[i+1:]) + 
            EU[i+1:]*np.cosh(s[i]-s[i+1:]))*ds).sum()
        Dphi[i] = ((EV[i+1:]*np.cosh(s[i]-s[i+1:]) + 
            EU[i+1:]*np.sinh(s[i]-s[i+1:]))*ds).sum()

    Dphi[(Dpsi<=0.)]=0.

    # total wave energy
    Dw = (Dphi[1:]+Dphi[:-1])/2.
    Dpsi2=(Dpsi[1:]+Dpsi[:-1])/2.
    kw=(k[1:]+k[:-1])/2.
    Kpsi= Dpsi2 - kw*(np.diff(Dpsi)/np.diff(k))
    Kphi= Dw - kw*(np.diff(Dphi)/np.diff(k))
    kw = (k[1:]+k[:-1])/2.
    Ew = Dw - kw*(np.diff(Dphi)/np.diff(k))
    Dwpsi=Ew-Dw

    # 10 km cutoff
    l=1./k
    fn=(l>=1.)
    k,Dpsi,Dphi = k[fn],Dpsi[fn],Dphi[fn] #,((EV+EU+EB)[fn])
    Dpsi,Dphi=(Dpsi[1:]+Dpsi[0:-1])/2.,(Dphi[1:]+Dphi[0:-1])/2. #,(E[1:]+E[0:-1])/2.
    k=(k[1:]+k[0:-1])/2.
    l=1./kw
    fn=(l>=1.)
    Kpsi,Kphi,kw=Kpsi[fn],Kphi[fn],kw[fn]

    return Kpsi,Kphi,kw

def spec_est2(A,d1,d2,win=True):

    """    computes 2D spectral estimate of A
           obs: the returned array is fftshifted
           and consistent with the f1,f2 arrays
           d1,d2 are the sampling rates in rows,columns   """
    
    import numpy as np

    l1,l2,l3 = A.shape
    df1 = 1./(l1*d1)
    df2 = 1./(l2*d2)
    f1Ny = 1./(2*d1)
    f2Ny = 1./(2*d2)

    f1 = np.arange(-f1Ny,f1Ny,df1)
    f2 = np.arange(-f2Ny,f2Ny,df2)
    
    if win == True:
        wx = np.matrix(np.hanning(l1))
        wy =  np.matrix(np.hanning(l2))
        window_s = np.repeat(np.array(wx.T*wy),l3).reshape(l1,l2,l3)
    else:
        window_s = np.ones((l1,l2,l3))

    an = np.fft.fft2(A*window_s,axes=(0,1))
    E = (an*an.conjugate()) / (df1*df2) / ((l1*l2)**2)
    E = np.fft.fftshift(E)
    E = E.mean(axis=2)

    return np.real(E),f1,f2,df1,df2,f1Ny,f2Ny

