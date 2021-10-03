import numpy as np
import netCDF4 as nc
from lorenz import lz_rk4, lz_rk4_adj
import matplotlib.pyplot as plt

bsigma, osigma=1., 1.
nx=40
tlength=1.
dt=0.0125
#nt=10  # simple testing
nt=int(tlength/dt) # testing purpose

xf=np.ones(nx)*8.

obsflag=np.zeros(nx, bool)
obsflag[::5]=1.  # obs spatial locations

xnature_name='xnature-obs.nc'
np.random.seed(10)
with nc.Dataset(xnature_name, 'r') as infile:
    x_fg=infile['xnature'][0]+np.random.normal(scale=bsigma, size=nx)

def cost_obs(x, obsloc, obsval, oflag):
    Jo=0.

    x0=x.copy()
    for it in range(nt):
        x0=lz_rk4(x0, xf, dt)
        if it+1 in obsloc:
            oloc=np.argmin(np.absolute(it+1-obsloc))
            Jo+= ( ((x0[oflag]-obsval[oloc][oflag])**2)/osigma**2 ).sum()
    return Jo*0.5

def grad_cost_obs(x, obsloc, obsval, oflag):
    xstate=np.zeros([nt+1, nx], float)

    # ----- forward run -----
    x0=x.copy()
    for it in range(nt):
        xstate[it]=x0.copy()
        x0=lz_rk4(x0, xf, dt)
    xstate[nt]=x0
    # ----- end forward run -----

    # ----- adjoint run -----
    xf_ad=np.zeros(nx)
    x0_ad=np.zeros(nx)
    xout_ad=np.zeros(nx)
    for it in range(nt, 0, -1):
        if it in obsloc:
            oloc=np.argmin(np.absolute(it-obsloc))
            x0_ad[oflag]+= (xstate[it]-obsval[oloc])[oflag] /osigma**2
        xout_ad=xout_ad+x0_ad
        x0_ad[()]=0.
        x0_ad, xout_ad, xf_ad=lz_rk4_adj(xstate[it-1], x0_ad, xout_ad, xf, xf_ad, dt)
    # ----- end adjoint run -----
    return x0_ad

def cost_mod(x):
    Jb=( ((x-x_fg)**2)/bsigma**2 ).sum()
    return Jb*0.5

def grad_cost_mod(x):
    return (x-x_fg)/bsigma**2

def cost(x):
    with nc.Dataset(xnature_name, 'r') as infile:
        obsloc=infile['obsloc'][()]
        obsval=infile['obsval'][()]
    Jo=cost_obs(x, obsloc, obsval, obsflag)
    Jb=cost_mod(x)
    return Jb+Jo

def grad_cost(x):
    with nc.Dataset(xnature_name, 'r') as infile:
        obsloc=infile['obsloc'][()]
        obsval=infile['obsval'][()]
    dx_o=grad_cost_obs(x, obsloc, obsval, obsflag)
    dx_b=grad_cost_mod(x)
    return dx_o+dx_b

def check_grad():

    with nc.Dataset(xnature_name, 'r') as infile:
        obsloc=infile['obsloc'][()]
        obsval=infile['obsval'][()]
    #J0=cost_obs(x_fg, obsloc, obsval, obsflag)
    #dJ=grad_cost_obs(x_fg, obsloc, obsval, obsflag)
    J0=cost(x_fg)
    dJ=grad_cost(x_fg)

    for i in range(1, 15):
        np.random.seed(10)
        pert=np.random.normal(size=nx)*(10.**-i)
        #J1=cost_obs(x_fg+pert, obsloc, obsval, obsflag)
        J1=cost(x_fg+pert)
        print( i, (J1-J0)/(dJ*pert).sum() )

def minimizer():
    from scipy.optimize import minimize
    res=minimize(cost, x_fg, method='L-BFGS-B', jac=grad_cost,
                 options={'gtol': 1.E-9, 'disp': True})
    return res

if __name__=='__main__':
    #check_grad()
    res=minimizer()
    x_an=res.x

    xsol_an=[]
    x0=x_an.copy()
    xsol_an.append(x0.copy())
    for it in range(nt):
        xout=lz_rk4(x0, xf, dt)
        xsol_an.append(xout.copy())
        x0=xout.copy()
    xsol_an=np.array(xsol_an)

    xsol_fg=[]
    x0=x_fg.copy()
    xsol_fg.append(x0.copy())
    for it in range(nt):
        xout=lz_rk4(x0, xf, dt)
        xsol_fg.append(xout.copy())
        x0=xout.copy()
    xsol_fg=np.array(xsol_fg)

    with nc.Dataset(xnature_name, 'r') as infile:
        xtrue=infile['xnature'][()]

    #clevs=np.arange(-15, 15+1.E-6, 1.5)
    z=xsol_an - xtrue
    plt.contourf(z, extend='both', cmap='coolwarm')#, levels=clevs)
    plt.colorbar()
    plt.savefig('test')
    
