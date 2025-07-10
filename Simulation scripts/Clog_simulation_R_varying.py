
from __future__ import absolute_import
import numpy as np
from numpy import pi,log
from scipy.special import erf
from clawpack import riemann 
import time
    
start1 = time.time()
# Function to output parameters \phi and \mu_w goes into _p file
def particle_fraction(state):
    """ Compute phi from q (mu_w)  and store in state.p."""
    import numpy as np
    from numpy import log
    from scipy.optimize import leastsq

    q = state.q[0,:]
    R = state.aux[2,:]
    mu1 = state.problem_data['mu1']
    D = state.problem_data['D']
    phim = state.aux[0,:]

    num_rp = state.q.shape[1]
    m = np.empty( (num_rp) )


    def F1(m1,mw):
        
        m = mw
        rc = m1/mw
        mu=m1

        if(mw>mu):
            F1 = 2*(mw*(-3 + 2*(mw - mu)**(1/2)) + 3*mu + 2*(mw - mu)**(1/2)*(3 + 2*mu) - 6*(1 + mu)*log(1 + (mw - mu)**(1/2)))/(3*mw**2) + (mu/mw)**2
        elif(mw==mu):
            F1 = mw**0
        else:
            F1 = 1000

        return(F1)
        
    def F2(mu1,muw):

        mu=mu1
        mw = muw
        F = (2*mu - mw)/(2*mw**2) - (mw - mu)**(1/2)/(7*mw) - (mu**2 - 1)/(2*mw**3) - mu/(2*mw**2) + 1/(6*mw) + (mu + 1)/(4*mw**2) - mu**3/(6*mw**4) - ((mw - mu)**(1/2)*(- 64*mu**3 + 140*mu**2*mw - 119*mu**2 - 70*mu*mw**2 + 210*mu*mw + 70*mu - 105*mw**2 + 105))/(105*mw**4) + ((8*mu - 7)*(mw - mu)**(1/2))/(35*mw**2) - (mu*(- mu**2 + 2*mu*mw - mw**2 + 1))/(2*mw**4) + (mu**2*(mu - 1))/(4*mw**4) + ((mw - mu)**(1/2)*(32*mu**2 - 70*mu*mw + 7*mu + 35*mw**2 - 35))/(105*mw**3) + (mu**2*(mu - mw)**2)/(4*mw**4)  + (log((-mw + mu + 1)*(mw - mu + 2*(mw - mu)**(1/2) + 1)/(mu - mw + 1))*(mu + 1)*(- mu**2 + 2*mu*mw - mw**2 + 1))/(2*mw**4) #+ (log((mw - mu + 2*(mw - mu)**(1/2) + 1)/(mu - mw + 1))*(mu + 1)*(mu - mw + 1)*(mw - mu + 1))/(2*mw**4)  
        
        F2 = np.where(mw>mu,F,0)
        return(F2)

    
    
    def F3(mu1,muw):
        m = muw
        rc = mu1/muw
        F = (rc**2*(rc - 1/2))/2 - rc/6 + (rc**2*(rc - 1)**2)/4 - (5*rc**4)/24 + 1/8
        F3 = np.where(muw>mu1,F,0)
        return(F3)
      
    
    def F4(mu1,muw,pm):
        mw=muw
        mu=mu1
        F = -(16*mu**2*(mw - mu)**(1/2) - 24*mw**2*(mw - mu)**(1/2) + 60*pm**3*(mw - mu)**(1/2) + 15*mu*mw**2 + 30*mu*pm**3 + 45*mw**2*pm - 30*mw*pm**3 - 5*mu**3 - 15*mw**2 - 10*mw**3 - 30*pm**3*log((mw - mu + 2*(mw - mu)**(1/2) + 1)) + 15*mu**2*pm**3 - 45*mw**2*pm**2 + 8*mu*mw*(mw - mu)**(1/2) - 30*mu*pm**3*log((mw - mu + 2*(mw - mu)**(1/2) + 1)) - 24*mu**2*pm*(mw - mu)**(1/2) + 40*mu*pm**3*(mw - mu)**(1/2) + 36*mw**2*pm*(mw - mu)**(1/2) + 20*mw*pm**3*(mw - mu)**(1/2) - 12*mu*mw*pm*(mw - mu)**(1/2))/(15*mw**2*pm**2)
        F4 = np.where(mw>mu,F,(1 - pm)**3/(pm**2))
        return(F4)
      
    
    def func(mw,*arg):
        mu1,R,phim,q = arg
        return ( q/((R**2)*phim) - F1(mu1,mw) )

    x_in = 4

    for x in range(num_rp):
        arg = (mu1,R[x], phim[x], q[x])
        X,ss = leastsq(func,x_in, args = arg)
        m[x]=X[0]
        x_in = m[x]


    phi = q/(R**2)
    # F_2 = F2(mu1,m)
    # F_3 = F3(mu1,m)
    # F_4 = F4(mu1,m,phim)
    # Flux = (phim*F_2)/(F_3 + D*F_4/(R**2))

    state.p[0,:] = phi
    state.p[1,:] = m
    # state.p[2,:] = Flux   # If needed, can output Flux as well
    
    

    
# Function to setup problem with relevant parameters and solver. Specify output folder by setting outdir 
def setup(use_petsc=0,kernel_language='Python',outdir='./_output_Rstepvary',solver_type='classic',disable_output=False):

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if kernel_language == 'Python': 
        riemann_solver = riemann.Mu_I.Full_model    # Assigning the riemann solver to the one we want to use 
    elif kernel_language == 'Fortran':
        riemann_solver = riemann.Mu_I.Full_model


    if solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver1D(riemann_solver)
    else:
        solver = pyclaw.ClawSolver1D(riemann_solver)
        solver.limiters = pyclaw.limiters.tvd.vanleer


    
    solver.kernel_language = kernel_language
        
    solver.bc_lower[0] = pyclaw.BC.extrap
    solver.bc_upper[0] = pyclaw.BC.extrap
    solver.aux_bc_lower[0]=pyclaw.BC.extrap
    solver.aux_bc_upper[0]=pyclaw.BC.extrap

    x = pyclaw.Dimension(-0.0,1.0,2500,name='x')    # Coarser mesh works for R varying case
    domain = pyclaw.Domain(x)
    num_eqn = 1
    num_aux=3
    state = pyclaw.State(domain,num_eqn,num_aux)
    state.mp = 2

  
    # Define variaiton in phi_m and R
    c = 0.0     # Kept phi_m constant
    phim0 = 0.65
    rc = 0.5
    phin =0.35  # Initial value of particle fraction in the whole channel
    
    # For step phim change
    b = c/2 
    a = phim0-b
    d = 10
    
    # For step R change
    b_r = rc/2 
    a_r = 1-b_r
    d_r = 10

    xc = state.grid.x.centers
    state.aux[0,:] = a + b*(erf(-(xc-0.5)*d))                   # phi_m(z)
    state.aux[1,:] = -2*b*d*np.exp(-((xc-0.5)*d)**2)/pi**(0.5)  # phi_mz
    state.aux[2,:] = a_r + b_r*(erf(-(xc-0.5)*d_r))             # R(z)
     
    state.q[0,:] = phin * (state.aux[2,:])**2                    #q = phi* R^2 
    
 
    state.problem_data['D'] = .01
    state.problem_data['mu1'] = 0.3
    
    claw = pyclaw.Controller()
    claw.tfinal = 14.0
    claw.solution = pyclaw.Solution(state,domain)

    claw.solver = solver
    claw.outdir = outdir
    claw.setplot = setplot
    claw.keep_copy = True
    claw.num_output_times=1400
    claw.output_style =1
    claw.write_aux_always = True
    claw.overwrite = True
    

    if disable_output:
        claw.output_format = None
    claw.compute_p = particle_fraction
    claw.write_aux_init = False
    claw.file_prefix_p = 'fort'

    return claw


# Can be used to plot using Iplotclaw - Clawpacks interactive plotting routine. We currently plot output files using matlab routine plotclaw1.m from visclaw package of Clawpack 
def setplot(plotdata):
    """ 
    Plot solution using VisClaw.
    """ 
    plotdata.clearfigures()  # clear any old figures,axes,items data

    # Figure for q[0]
    plotfigure = plotdata.new_plotfigure(name='q[0]', figno=0)


    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = [0, 3]
    plotaxes.title = 'q[0]'
    plotaxes.afteraxes = add_true_solution
  

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = 0
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    plotitem.kwargs = {'linewidth': .5}
    
    return plotdata

    
def add_true_solution(current_data):
   import matplotlib.pyplot as plt
   from pylab import sin, plot, title
   phim = current_data.aux
   x = current_data.x
   t = current_data.t
   q = current_data.q
   aux = current_data.aux
   phim = aux[0,:]
   R = aux[2,:]
   
   
   phi = np.empty(q.shape)
   phi = q[0,:]/(R**2)   
   
   #plot(x,sin(x),'r')
   plot(x,phim,'r')
   plot(x,R,'b')
   plot(x,phi,'g' )
   plt.xlabel(r'$z$')
   plt.ylabel(r'$\overline{\phi}$')
   title("Solution at time t = %10.4e" % t, fontsize=10)
    


if __name__=="__main__":
    #print('Line 95')
    #x0=2
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup,setplot)
