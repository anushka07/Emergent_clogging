r"""
Riemann solvers for Burgers equation

.. math:: 
    \frac{\partial R^2\phi}{\partial t} + \frac{\partial R^2 us \phi}{\partial z} = 0.

:Author:
    Anushka Herale(2025): Initial version
"""

from __future__ import absolute_import
num_eqn = 1
num_waves = 3

import numpy as np
from numpy import pi, log
from scipy.optimize import leastsq


def Full_model(q_l,q_r,aux_l,aux_r,problem_data):
    r"""
    Riemann solver for mass conservation of particle fraction in 1d
    See :ref:`pyclaw_rp` for more details.
    
    """
    
    
    num_rp = q_l.shape[1]
    m_ml = np.empty( (num_eqn, num_rp) )
    m_mr = np.empty( (num_eqn, num_rp) )
    m_pl = np.empty( (num_eqn, num_rp) )
    m_pr = np.empty( (num_eqn, num_rp) )
    m_rl = np.empty( (num_eqn, num_rp) )
    m_rr = np.empty( (num_eqn, num_rp) )

    # Output arrays
    wave = np.empty( (num_eqn, num_waves, num_rp) )
    s = np.empty( (num_waves, num_rp) )
    amdq = np.zeros( (num_eqn, num_rp) )
    apdq = np.zeros( (num_eqn, num_rp) )

   
    # Basic solve
    
    # dFlux/dz = dFlux/dphim phim_z + dFlux/dR R_z + dFlux/dmu_w muw_z  # 3 waves, 1 equation
    
    q=.5 * (q_r[0,:] + q_l[0,:])
    phim = .5*(aux_l[0,:] + aux_r[0,:])
    R = .5*(aux_l[2,:] + aux_r[2,:])
    
    phim_r = aux_r[0,:]
    phim_l = aux_l[0,:]
    R_r = aux_r[2,:]
    R_l = aux_l[2,:]
    
    D = problem_data['D']
    mu1 = problem_data['mu1']

    def F1(m1,mw):

        mu=m1

        if(mw>mu):
            F1 = 2*(mw*(-3 + 2*(mw - mu)**(1/2)) + 3*mu + 2*(mw - mu)**(1/2)*(3 + 2*mu) - 6*(1 + mu)*log(1 + (mw - mu)**(1/2)))/(3*mw**2) + (mu/mw)**2
        elif(mw==mu):
            F1 = mw**0
        else:
            F1 = 1000

        return(F1)


    
    def func(mw,*arg):
        
        mu1,R,phim,q = arg

        return ( q/((R**2)*phim) - F1(mu1,mw) )
                

    
    
    x_l=4
    x_r=4


    for x in range(num_rp):
        
        g = 1
        arg = (mu1,R[x], phim[x], q_l[0,x] ) #muw definition since q = phi_m*R**2*f1(muw) from fsolve
        X,ss = leastsq(func,x_l, args = arg)
        m_ml[0,x]=X[0]
        
        x_l = m_ml[0,x]
        if(m_ml[0,x]<mu1):
            raise Exception("Wrong output:mw<mu")
        
        arg = (mu1,R[x], phim[x], q_r[0,x])
        g = 2
        if(x==0):
            X,ss = leastsq(func,x_l, args = arg)
            m_mr[0,x]=X[0]
            x_r = m_mr[0,x]
        else:
            X,ss = leastsq(func,x_r, args = arg)
            m_mr[0,x]=X[0]
            x_r = m_mr[0,x]

        # mu_w for flux_phim

        arg = (mu1,R[x], phim_l[x], q[x] )
        X,ss = leastsq(func,x_l, args = arg)
        m_pl[0,x]=X[0]

        arg = (mu1,R[x], phim_r[x], q[x])
        X,ss = leastsq(func,x_r, args = arg)
        m_pr[0,x]=X[0]

        # mu_w for flux_R

        arg = (mu1,R_l[x], phim[x], q[x] )
        X,ss = leastsq(func,x_l, args = arg)
        m_rl[0,x]=X[0]

        arg = (mu1,R_r[x], phim[x], q[x])
        X,ss = leastsq(func,x_r, args = arg)
        m_rr[0,x]=X[0]




    def F2(mu,mw):
        F = (2*mu - mw)/(2*mw**2) - (mw - mu)**(1/2)/(7*mw) - (mu**2 - 1)/(2*mw**3) - mu/(2*mw**2) + 1/(6*mw) + (mu + 1)/(4*mw**2) - mu**3/(6*mw**4) - ((mw - mu)**(1/2)*(- 64*mu**3 + 140*mu**2*mw - 119*mu**2 - 70*mu*mw**2 + 210*mu*mw + 70*mu - 105*mw**2 + 105))/(105*mw**4) + ((8*mu - 7)*(mw - mu)**(1/2))/(35*mw**2) - (mu*(- mu**2 + 2*mu*mw - mw**2 + 1))/(2*mw**4) + (mu**2*(mu - 1))/(4*mw**4) + ((mw - mu)**(1/2)*(32*mu**2 - 70*mu*mw + 7*mu + 35*mw**2 - 35))/(105*mw**3) + (mu**2*(mu - mw)**2)/(4*mw**4)  + (log((-mw + mu + 1)*(mw - mu + 2*(mw - mu)**(1/2) + 1)/(mu - mw + 1))*(mu + 1)*(- mu**2 + 2*mu*mw - mw**2 + 1))/(2*mw**4) 
        F2 = np.where(mw>mu,F,0)
        return(F2)
    
    
    
    def F3(mu1,muw):
        rc = mu1/muw
        F = (rc**2*(rc - 1/2))/2 - rc/6 + (rc**2*(rc - 1)**2)/4 - (5*rc**4)/24 + 1/8
        F3 = np.where(muw>mu1,F,0)
        return(F3)
      
    
    def F4(mu,mw,pm):
        # K(phi) = (1-phi)^3/phi^2
        F = -(16*mu**2*(mw - mu)**(1/2) - 24*mw**2*(mw - mu)**(1/2) + 60*pm**3*(mw - mu)**(1/2) + 15*mu*mw**2 + 30*mu*pm**3 + 45*mw**2*pm - 30*mw*pm**3 - 5*mu**3 - 15*mw**2 - 10*mw**3 - 30*pm**3*log((mw - mu + 2*(mw - mu)**(1/2) + 1)) + 15*mu**2*pm**3 - 45*mw**2*pm**2 + 8*mu*mw*(mw - mu)**(1/2) - 30*mu*pm**3*log((mw - mu + 2*(mw - mu)**(1/2) + 1)) - 24*mu**2*pm*(mw - mu)**(1/2) + 40*mu*pm**3*(mw - mu)**(1/2) + 36*mw**2*pm*(mw - mu)**(1/2) + 20*mw*pm**3*(mw - mu)**(1/2) - 12*mu*mw*pm*(mw - mu)**(1/2))/(15*mw**2*pm**2)
        F4 = np.where(mw>mu,F,(1 - pm)**3/(pm**2))

        return(F4)
      

    def Flux(pm,mu1,muw,D,R):

        f2 = F2(mu1,muw)
        f3 = F3(mu1,muw)
        f4 = F4(mu1,muw,pm)
        
        Flux =(1/pi)*pm*f2/(f3 + D*f4/(R**2))

        return(Flux)
        
        

    
    
    wave[0,0,:] = q_r - q_l            #Wave mu_w 
    wave[0,1,:] = phim_r - phim_l      #Wave phim
    wave[0,2,:] = R_r - R_l            #Wave R
    
    
    ## Calculate fluxes on both sides of the cell edge
    
    #mu_w varying across cell edge and R, phim const = avg
    Flux_muw_R = Flux(phim,mu1,m_mr,D,R)
    Flux_muw_L = Flux(phim,mu1,m_ml,D,R)
    
    #phim varying across cell edge and R, mu_w const = avg
    Flux_phim_R = Flux(phim_r,mu1,m_pr,D,R)
    Flux_phim_L = Flux(phim_l,mu1,m_pl,D,R)
    
    #R varying across cell edge and phim, mu_w const = avg
    Flux_R_R = Flux(phim,mu1,m_rr,D,R_r)
    Flux_R_L = Flux(phim,mu1,m_rl,D,R_l)
    
    
    #Ensure if phi=phim then Flux F = 0 since it is clogged
    Fm_r = np.where(m_mr>mu1,Flux_muw_R,0)
    Fm_l = np.where(m_ml>mu1,Flux_muw_L,0)
    
    Fp_r = np.where(m_pr>mu1,Flux_phim_R,0)
    Fp_l = np.where(m_pl>mu1,Flux_phim_L,0)
    
    Fr_r = np.where(m_rr>mu1,Flux_R_R,0)
    Fr_l = np.where(m_rl>mu1,Flux_R_L,0)

    
    ##Calculate the speeds of each wave - partial dFlux/d*var* 
    x=0

    for x in range(num_rp):
        if(abs(q_r[0,x]-q_l[0,x])>1e-10):
            s[0,x] =(Fm_r[0,x] - Fm_l[0,x])/(q_r[0,x] - q_l[0,x])
        else:
            s[0,x] =0

        if(abs(phim_r[x]-phim_l[x])>1e-10):
            s[1,x] =(Fp_r[0,x] - Fp_l[0,x])/(phim_r[x] - phim_l[x])
        else:
            s[1,x] =0

        if(abs(R_r[x]-R_l[x])>1e-10):
            s[2,x] =(Fr_r[0,x] - Fr_l[0,x])/(R_r[x] - R_l[x])
        else:
            s[2,x] =0



# define the A del Q left and right going fluctuations 
# which means just seperating the +ve and -ve velocities
  
    s_index = np.zeros((2,num_rp))
    for mw in range(num_waves):
      s_index[0,:] = s[mw,:]
      amdq[0,:] += np.min(s_index,axis=0) * wave[0,mw,:]  
      apdq[0,:] += np.max(s_index,axis=0) * wave[0,mw,:]
    

    return wave, s, amdq, apdq
