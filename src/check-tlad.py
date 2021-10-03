import numpy as np
from lorenz import lorenzn, lorenzn_tlm, lorenzn_adj, \
    lz_rk4, lz_rk4_tlm, lz_rk4_adj

def check_tlm():

    nx=40
    dt=0.0125

    F=8.
    xf=np.ones(nx)*F

    #np.random.seed(5)
    x0=F/2. * np.random.normal(size=nx)

    xout=[]
    xout.append(x0.copy())
    for ip in range(0, 15):
        pt=10.**(-ip)
        xout1=lz_rk4(x0*(1.+pt), xf*(1.+pt), dt)

        xout0, xout_tl=lz_rk4_tlm(x0, x0*pt, xf, xf*pt, dt)

        idx=0
        print( (xout1[idx]-xout0[idx])/xout_tl[idx] )

def check_adj():
    nx=40
    dt=0.0125

    F=8.
    xf=np.ones(nx)*F

    #np.random.seed(5)
    x0=F/2. * np.random.normal(size=nx)
    xic=x0.copy()

    pt=0.1
    x0_tl, xf_tl=x0*pt, xf*pt

    # ----- multiple time steps -----
    # TLM run
    xstate=[]
    for it in range(100):
        xstate.append(x0.copy())
        xout, xout_tl=lz_rk4_tlm(x0, x0_tl, xf, xf_tl, dt)
        x0_tl=xout_tl.copy()
        x0=xout.copy()
    lhs=(x0_tl**2).sum()

    # ADJ run
    xout_ad=xout_tl*0.
    x0_ad=x0_tl.copy()
    xf_ad=xf_tl*0.
    for it in range(99, -1, -1):
        xout_ad=xout_ad+x0_ad
        x0_ad[()]=0.
        x0_ad, xout_ad, xf_ad=lz_rk4_adj(xstate[it], x0_ad, xout_ad, xf, xf_ad, dt)

    rhs=(x0_ad*xic*pt).sum() + (xf_ad*xf_tl).sum()

    print(lhs, rhs, (lhs-rhs)/lhs)


if __name__=='__main__':
    check_adj()
