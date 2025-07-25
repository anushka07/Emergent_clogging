# Emergent_clogging
Repository with all Python scripts and Riemann solvers used to solve the clogging problem. This code should be used within the **clawpack package** [[1]](#1).
The riemann solver file should be uploaded into PATH/clawpack/riemann/ and \_init_.py, setup.py and static.py need to be updated to add this riemann solver. Just add the following lines to the corresponding python scripts as indicated - <br>
<ol>
<li> 
  
`_init_.py` <br>
<ul>
  
  `from . import Mu_I` <br> 
</ul>
</li>

 <li> 
   
`setup.py` <br> 
<ul>

`one_d_riemann = ['Mu_I',`<br>
                  `]`
 </ul>                 
<li> 

  
`static.py` <br>
<ul>

`num_eqn = {'Mu_I' : 1,` <br>
        `'Full_model' : 1,` <br>
        `}`<br>
`num_waves = {'Mu_I' : 3,`<br>
        `'Full_model' : 3,`<br>
        `}`
</ul>
</ol>

Once the riemann solver is added you can run the required python script from the `Simulation scripts` folder of the repository. The script `Clog_simulation_Step_ini_R_varying.py` should be run after running the MATLAB script `Steady_phi_R_vary.m` which generates the required `.mat` file acting as the initial condition of the simulation. The resulting output files can be made into a simulation video by modifying `plotclaw1.m` and `afterframe.m` files from `visclaw`. Below are two simulations that give rise to a "non-clogged" and a "clogged" steady state for input `phi = 0.35` and `phi = 0.45` respectively.



https://github.com/user-attachments/assets/28f74b2f-6bf9-4bf5-99be-4f61a30e9d02




https://github.com/user-attachments/assets/680ce1fe-2159-4a56-84c8-2939e240b54b



## References
<a id="1">[1]</a> 
Ketcheson, David I. and Mandli, Kyle T. and Ahmadia, Aron J. and Alghamdi, Amal and Quezada de Luna, Manuel and Parsani, Matteo and Knepley, Matthew G. and Emmett, Matthew (2012). 
PyClaw: Accessible, Extensible, Scalable Tools for Wave Propagation Problems. 
SIAM Journal on Scientific Computing,34(4), C210--C231.
