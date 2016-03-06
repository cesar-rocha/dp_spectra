import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LogNorm

import seawater as sw

plt.rcParams.update({'font.size': 12
    , 'legend.markerscale': 1., 'axes.titlesize': 12, 'axes.labelsize' : 12,
      'legend.fontsize' : 8,'legend.handlelength': 3})

plt.rc('xtick', labelsize=12) 
plt.rc('ytick', labelsize=12)


plt.close('all')

def calc_ispec(k,l,E):
    """ calculates isotropic spectrum from 2D spectrum """

    dk,dl = k[1,]-k[0],l[1]-l[0]
    l,k = np.meshgrid(l,k)
    wv = np.sqrt(k**2 + l**2)

    if k.max()>l.max():
        kmax = l.max()
    else:
        kmax = k.max()

    # create radial wavenumber
    dkr = np.sqrt(dk**2 + dl**2)
    kr =  np.arange(dkr/2.,kmax+dkr,dkr)
    ispec = np.zeros(kr.size)

    for i in range(kr.size):
        fkr =  (wv>=kr[i]-dkr/2) & (wv<=kr[i]+dkr/2)
        dth = pi / (fkr.sum()-1)
        ispec[i] = E[fkr].sum() * kr[i] * dth

    return kr, ispec

def add_second_axis(ax1):
    """ Add a x-axis at the top of the spectra figures """
    ax2 = ax1.twiny() 
    ax2.set_xscale('log')
    ax2.set_xlim(ax1.axis()[0], ax1.axis()[1])
    kp = 1./np.array([200.,100.,40.,20.,10.,5.])
    lp=np.array([200,100,40,20,10,5])
    ax2.set_xticks(kp)
    ax2.set_xticklabels(lp)
    plt.xlabel('Wavelength [km]')


#wv = np.load('wavenumber_frequency_spec.npz')
wv = np.load('wavenumber_frequency_spec_full.npz')
wvf = np.load('frequency_spec_full.npz')

omg = wv['omg']
E = wv['E']
k = wv['k']
l = wv['l']

# isotropic
Eiso = np.empty((292,omg.size))

for i in range(omg.size):
    kiso, Eiso[:,i] = calc_ispec(k,l,E[:,:,i])

# linear dispersion relationship
kr = 2*pi*kiso*1.e-3
kr2 = kr**2

m = np.logspace(-3, 0., 500)

#omgr = 2*pi*omg/8600.

N2 = (1.7e-5)
f2 = sw.f(59.247177)**2

b = 1.e3
jmax = 100

#for j in range(jmax):
for j in range(m.size):
    #m = (pi*j)/b
    m2 = m[i]**2
    omgr = np.sqrt( (f2*m2 + N2*kr2)/(kr2 + m2) )
    krp = 1.e3*kr/(2*pi)
    
    if j == 0:
        omgrp = 3600*omgr/(2*pi)
    else:
        omgrp = np.vstack([omgrp,3600*omgr/(2*pi)])

omg_int = np.empty(kr.size)
for i in range(kr.size):
    omg_int[i] = np.trapz(omgrp[:,i], x=m)


# auxiliary frequencies
omg_m2 = 1./12.4
omg_f  = 1./13.9
omg_k1 = 1./23.93
omg_o1 = 1./25.82 
ks = np.array([1.e-3,1.])

fig = plt.figure(figsize=(8.27/2+1.5,11.69/3+1.))
ax1 = fig.add_subplot(111)
plt.pcolormesh(kiso[1:],omg[1:],Eiso.T[1:,1:],shading='flat',
                cmap="viridis",norm = LogNorm())
ax1.plot(kiso,omg_int,linewidth=1.75,color='k')
ax1.plot(ks,[omg_m2,omg_m2],'w--',linewidth=1.)
ax1.plot(ks,[omg_f,omg_f],'w--',linewidth=1.)
ax1.plot(ks,[omg_k1,omg_k1],'w--',linewidth=1.)
ax1.plot(ks,[omg_o1,omg_o1],'w--',linewidth=1.)

ax1.text(0.1,omg_m2+.0075,'M2',color='w')
ax1.text(0.1,omg_f-.013,'f',color='w')
ax1.text(0.1,omg_k1+.002,'K1',color='w')
ax1.text(0.1,omg_o1-.0075,'O1',color='w')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_ylim(0.005,1.)
ax1.set_xlim(1/400.,.2)
plt.clim([1.e-5,1.e0])
cb = plt.colorbar(label=r"KE density [m$^2$ s$^{-2}$/(cpkm $\times$ cph)]")
plt.xlabel(r'Horizontal wavenumber [cpkm]')
plt.ylabel(r'Frequency [cph]')
add_second_axis(ax1)
plt.savefig('figs/wavenumber_freq_spec_100m')
plt.savefig('figs/wavenumber_freq_spec_100m.eps',rasterized=True)
plt.savefig('figs/wavenumber_freq_spec_100m.pdf',rasterized=True)

# now integrate the spectrum
domg = omg[1]
dk = kiso[1]
Eiso2 = Eiso.copy()
Etot = Eiso2.sum()*domg*dk

f1 = omg>omg_f
#fn1 = omg<0.9*omg_f
#f1 = omg>=1./24

E1 = (Eiso2[:,f1]).sum()*domg*dk

f2 = (kiso>=1./40)&(kiso<1./10)
Eiso3 = Eiso2[:,f1]   # wavelegnths smaller than 40 km
Eiso4 = Eiso2[:,~f1]  # wavelengths larger than 40 km
E2 = (Eiso3[f2,:]).sum()*domg*dk
E3 = (Eiso4[f2,:]).sum()*domg*dk


#import seaborn as sns
#sns.set_style("whitegrid")
#plt.rcParams.update({'font.size': 25, 'legend.handlelength'  : 1.5
#    , 'legend.markerscale': 1., 'axes.titlesize': 25, 'axes.labelsize' : 25,
#      'legend.fontsize' : 35})
#
#plt.rc('xtick', labelsize=25) 
#plt.rc('ytick', labelsize=25)

plt.rcParams.update({'font.size': 10
    , 'legend.markerscale': 1., 'axes.titlesize': 10, 'axes.labelsize' : 10,
      'legend.fontsize' : 6,'legend.handlelength': 3})


def add_second_axis2(ax1):
    """ Add a x-axis at the top of the spectra figures """
    ax2 = ax1.twiny() 
    ax2.set_xscale('log')
    ax2.set_xlim(ax1.axis()[0], ax1.axis()[1])
    kp = 1./np.array([800,250.,100.,50.,25.,10.,5.])
    lp = np.array([800,250,100,50,25,10,5])
    ax2.set_xticks(kp)
    ax2.set_xticklabels(lp)
    plt.xlabel('Wavelength [km]')

ymin,ymax = 1.e-8,1.e1

Eiso2 = 4*Eiso2

klarge = np.array([1./1100.,1./245.])
kmeso = np.array([1./245.,1./75.])
ksubmeso = np.array([1./75.,1./7.5])
kdiss = np.array([1./7.5,2.])

fig = plt.figure(figsize=(8.27/2+1.5,11.69/3+.5))
ax1 = fig.add_subplot(111)


ax1.fill_between(x=klarge,y1=ymin,y2=ymax,color='m',alpha=.2)
ax1.text(1./980,.5e-4,'Large scale')
ax1.text(1./700,.25e-4,'range')

ax1.fill_between(x=kmeso,y1=ymin,y2=ymax,color='r',alpha=.2)
ax1.text(1./250,.5e-4,' Mesoscale')
ax1.text(1./195,.25e-4,'range')

ax1.fill_between(x=ksubmeso,y1=ymin,y2=ymax,color='g',alpha=.2)
ax1.text(1./50,.5e-4,'Submesoscale')
ax1.text(1./29.5,.25e-4,'range')

ax1.fill_between(x=kdiss,y1=ymin,y2=ymax,color='b',alpha=.2)
ax1.text(1./4.75,.5e-4,'Dissipation')
ax1.text(1./3.5,.25e-4,'range')

ax1.loglog(kiso,Eiso2.sum(axis=1)*domg,linewidth=3,color='0.4')

plt.xlabel(r'Wavenumber [cpkm]')
plt.ylabel(r'KE density [m$^2$ s$^{-2}$]')
plt.xlim(1.e-3,1.)
plt.ylim(1.e-6,1.e1)
add_second_axis2(ax1)
plt.savefig('figs/wavenumber_ke_spec_ranges.pdf')


fig = plt.figure(figsize=(8.27/2+1.5,11.69/3+.5))
ax1 = fig.add_subplot(111)


#ax1.fill_between(x=klarge,y1=ymin,y2=ymax,color='m',alpha=.2)
#ax1.text(1./980,.5e-4,'Large scale')
#ax1.text(1./700,.25e-4,'range')

ax1.fill_between(x=kmeso,y1=ymin,y2=ymax,color='r',alpha=.2)
ax1.text(1./250,.5e-4,' Mesoscale')
ax1.text(1./195,.25e-4,'range')

#ax1.fill_between(x=ksubmeso,y1=ymin,y2=ymax,color='g',alpha=.2)
#ax1.text(1./50,.5e-4,'Submesoscale')
#ax1.text(1./29.5,.25e-4,'range')

#ax1.fill_between(x=kdiss,y1=ymin,y2=ymax,color='b',alpha=.2)
#ax1.text(1./4.75,.5e-4,'Dissipation')
#ax1.text(1./3.5,.25e-4,'range')

ax1.loglog(kiso,Eiso2.sum(axis=1)*domg,linewidth=3,color='0.4')

plt.xlabel(r'Wavenumber [cpkm]')
plt.ylabel(r'KE density [(m$^2$ s$^{-2}$)/cpkm]')


plt.xlim(1.e-3,1.)
plt.ylim(1.e-6,1.e1)
add_second_axis2(ax1)
plt.savefig('figs/wavenumber_ke_spec_mesoscale.pdf')


fig = plt.figure(figsize=(8.27/2+1.5,11.69/3+.5))
ax1 = fig.add_subplot(111)


ax1.fill_between(x=klarge,y1=ymin,y2=ymax,color='m',alpha=.2)
ax1.text(1./980,.5e-4,' Large scale')
ax1.text(1./700,.25e-4,'range')

ax1.fill_between(x=kmeso,y1=ymin,y2=ymax,color='r',alpha=.2)
ax1.text(1./250,.5e-4,' Mesoscale')
ax1.text(1./195,.25e-4,'range')

#ax1.fill_between(x=ksubmeso,y1=ymin,y2=ymax,color='g',alpha=.2)
#ax1.text(1./50,.5e-4,'Submesoscale')
#ax1.text(1./29.5,.25e-4,'range')

#ax1.fill_between(x=kdiss,y1=ymin,y2=ymax,color='b',alpha=.2)
#ax1.text(1./4.75,.5e-4,'Dissipation')
#ax1.text(1./3.5,.25e-4,'range')

ax1.loglog(kiso,Eiso2.sum(axis=1)*domg,linewidth=3,color='0.4')

plt.xlabel(r'Wavenumber [cpkm]')
plt.ylabel(r'KE density [(m$^2$ s$^{-2}$)/cpkm]')


plt.xlim(1.e-3,1.)
plt.ylim(1.e-6,1.e1)
add_second_axis2(ax1)
plt.savefig('figs/wavenumber_ke_spec_mesoscale_large.pdf')


fig = plt.figure(figsize=(8.27/2+1.5,11.69/3+.5))
ax1 = fig.add_subplot(111)


ax1.fill_between(x=klarge,y1=ymin,y2=ymax,color='m',alpha=.2)
ax1.text(1./980,.5e-4,' Large scale')
ax1.text(1./700,.25e-4,'range')

ax1.fill_between(x=kmeso,y1=ymin,y2=ymax,color='r',alpha=.2)
ax1.text(1./250,.5e-4,' Mesoscale')
ax1.text(1./195,.25e-4,'range')

ax1.fill_between(x=ksubmeso,y1=ymin,y2=ymax,color='g',alpha=.2)
ax1.text(1./50,.5e-4,'Submesoscale')
ax1.text(1./29.5,.25e-4,'range')

#ax1.fill_between(x=kdiss,y1=ymin,y2=ymax,color='b',alpha=.2)
#ax1.text(1./4.75,.5e-4,'Dissipation')
#ax1.text(1./3.5,.25e-4,'range')

ax1.loglog(kiso,Eiso2.sum(axis=1)*domg,linewidth=3,color='0.4')

plt.xlabel(r'Wavenumber [cpkm]')
plt.ylabel(r'KE density [(m$^2$ s$^{-2}$)/cpkm]')



plt.xlim(1.e-3,1.)
plt.ylim(1.e-6,1.e1)
add_second_axis2(ax1)
plt.savefig('figs/wavenumber_ke_spec_mesoscale_large_sub.pdf')


fig = plt.figure(figsize=(8.27/2+1.5,11.69/3+.5))
ax1 = fig.add_subplot(111)

ax1.fill_between(x=klarge,y1=ymin,y2=ymax,color='m',alpha=.2)
ax1.text(1./980,.5e-4,' Large scale')
ax1.text(1./700,.25e-4,'range')

ax1.fill_between(x=kmeso,y1=ymin,y2=ymax,color='r',alpha=.2)
ax1.text(1./250,.5e-4,' Mesoscale')
ax1.text(1./195,.25e-4,'range')
ax1.fill_between(x=ksubmeso,y1=ymin,y2=ymax,color='g',alpha=.2)
ax1.text(1./50,.5e-4,'Submesoscale')
ax1.text(1./32.5,.25e-4,'range')
ax1.fill_between(x=kdiss,y1=ymin,y2=ymax,color='b',alpha=.2)
ax1.text(1./3.75,.5e-4,'Dissipation')
ax1.text(1./3.,.25e-4,'range')
ax1.loglog(kiso,Eiso2.sum(axis=1)*domg,linewidth=3,color='0.4')
plt.xlabel(r'Wavenumber [cpkm]')
plt.ylabel(r'KE density [(m$^2$ s$^{-2}$)/cpkm]')


plt.xlim(1.e-3,1.)
plt.ylim(1.e-6,1.e1)
add_second_axis2(ax1)
plt.savefig('figs/wavenumber_ke_spec_mesoscale_large_sub_diss.pdf')


fig = plt.figure(figsize=(8.27/2+1.5,11.69/3+.5))
ax1 = fig.add_subplot(111)
ax1.loglog(kiso,Eiso2.sum(axis=1)*domg,linewidth=3,color='0.4')
plt.xlabel(r'Wavenumber [cpkm]')
plt.ylabel(r'KE density [(m$^2$ s$^{-2}$)/cpkm]')


plt.xlim(1.e-3,1.)
plt.ylim(1.e-6,1.e1)
add_second_axis2(ax1)
plt.savefig('figs/wavenumber_ke_noshade.pdf')


