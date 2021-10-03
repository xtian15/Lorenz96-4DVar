import numpy as np
import netCDF4 as nc
from lorenz import lz_rk4

if __name__=='__main__':

    nx=40
    tlength=3
    dt=0.0125
    nt=1 #int(tlength/dt)

    F=8.
    xf=np.ones(nx)*F

    np.random.seed(5)
    x0=F/2. * np.random.normal(size=nx)

    xout=[]
    xout.append(x0.copy())
    for it in range(nt):
        x1=lz_rk4(x0, xf, dt)
        xout.append(x1.copy())
        x0=x1.copy()
    xout=np.array(xout)

    osigma=1.
    obsloc=[]
    obsval=[]
    interval=1
    for i in range(interval, nt+1, interval):
        obsloc.append(i)
        obs=xout[i,:]+np.random.normal(scale=osigma, size=nx)
        obsval.append(obs)
    obsloc=np.array(obsloc)
    obsval=np.array(obsval)

    import matplotlib.pyplot as plt
    y=xout[:,0]
    plt.plot(np.arange(y.shape[0]), y, color='b', lw=1.)
    plt.plot(obsloc, obsval[:,0], ls='none', marker='o',
             mfc='none', mec='r', ms=2.)
    plt.savefig('test')

    fname='xnature-obs.nc'
    with nc.Dataset(fname, 'w') as of:
        of.createDimension('nx', nx)
        of.createDimension('nt', xout.shape[0])
        dim=('nt', 'nx')
        of.createVariable('xnature', 'f4', dim)[:]=xout
    
        of.createDimension('nobs', obsloc.shape[0])
        odim=('nobs', 'nx')
        of.createVariable('obsloc', 'i4', ('nobs'))[:]=obsloc
        of.createVariable('obsval', 'f4', odim)[:]=obsval
        
