import numpy as np
import netCDF4 as nc
from lorenz import lz_rk4, lz_rk4_tlm
import matplotlib.pyplot as plt

if __name__=='__main__':

    nx=40
    tlength=10
    dt=0.0125
    nt=int(tlength/dt)

    F=8.
    xf=np.ones(nx)*F

    np.random.seed(5)
    xic=F/2. * np.random.normal(size=nx)
    dx=np.random.normal(size=nx)*xic*1.E-3

    x0=xic+dx
    x1out=[]
    x1out.append(x0.copy())
    for it in range(nt):
        x1=lz_rk4(x0, xf, dt)
        x1out.append(x1.copy())
        x0=x1.copy()
    x1out=np.array(x1out)

    x0=xic.copy()
    x0_tl=dx.copy()
    x0out, xtlout=[], []
    x0out.append(x0.copy())
    xtlout.append(x0_tl.copy())
    for it in range(nt):
        x1, x1_tl=lz_rk4_tlm(x0, x0_tl, xf, xf.copy()*0., dt)
        x0out.append( x1.copy())
        xtlout.append(x1_tl.copy())
        x0_tl=x1_tl.copy()
        x0   =x1.copy()
    x0out, xtlout=np.array(x0out), np.array(xtlout)

    x=x1out-x0out
    y=xtlout
    diff=np.sqrt( ((x-y)**2).mean(axis=1) )

    plt.plot(diff[0:300])
    plt.savefig('test')
    exit()

    idx=240
    x=x1out[idx]-x0out[idx]
    y=xtlout[idx]
    plt.plot(x, y, ls='none', marker='.', mec='none', mfc='k', ms=3.)
    plt.savefig('test')

