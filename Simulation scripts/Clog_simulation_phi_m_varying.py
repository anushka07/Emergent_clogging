

r"""
=========================

Solve the Mass conservation equation:

.. math:: 
    \frac{\partial R^2\phi}{\partial t} + \frac{\partial R^2 us \phi}{\partial z} = 0.

Here we solve for q= R^2\phi. Since R is constant in this simulation q = \phi.
"""
from __future__ import absolute_import
import numpy as np
from numpy import pi,log
from scipy.special import erf
from clawpack import riemann  

    

    
# Function to setup problem with relevant parameters and solver. Specify output folder by setting outdir 
def setup(use_petsc=0,kernel_language='Python',outdir='./_output_phimstepvary_',solver_type='classic',disable_output=False):

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

    x = pyclaw.Dimension(0.0,1.0,5000,name='x') # Domain definition , z from 0 to 1 
    domain = pyclaw.Domain(x)
    num_eqn = 1
    num_aux=3
    state = pyclaw.State(domain,num_eqn,num_aux)


  
    # Define variaiton in phi_m and R
    c = 0.2 
    phim0 = 0.8
    rc = 0.0    # Kept R const in this case
    phin = 0.35 # Initial value of particle fraction in the whole channel
    
    # For step phi_m change
    b = c/2 
    a = phim0-b
    d = 10

    xc = state.grid.x.centers
    state.aux[0,:] = a + b*(erf(-(xc-0.5)*d))                   # phi_m(z)
    state.aux[1,:] = -2*b*d*np.exp(-((xc-0.5)*d)**2)/pi**(0.5)  # phi_mz
    state.aux[2,:] = 1 - rc*np.exp(-((xc-0.5)*5)**2)            # R(z)
     
    state.q[0,:] = phin * (state.aux[2,:])**2                    #q = phi* R^2 
    
 
    state.problem_data['D'] = .1
    state.problem_data['mu1'] = 0.5
    
    
    claw = pyclaw.Controller()
    claw.tfinal = 10.0
    claw.solution = pyclaw.Solution(state,domain)

    claw.solver = solver
    claw.outdir = outdir
    claw.setplot = setplot
    claw.keep_copy = True
    claw.num_output_times=1000
    claw.output_style =1
    claw.write_aux_always = True
    claw.overwrite = True
    

    if disable_output:
        claw.output_format = None
    claw.write_aux_init = False

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
  
   

   plot(x,phim,'r')
   plot(x,R,'b')
   plot(x,phi,'g' )
   plt.xlabel(r'$z$')
   plt.ylabel(r'$\overline{\phi}$')
   title("Solution at time t = %10.4e" % t, fontsize=10)
   


if __name__=="__main__":

    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup,setplot)
