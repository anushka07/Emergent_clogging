# Emergent_clogging
Repository with all Python scripts and Riemann solvers used to solve clogging problem. This code is to be used within the **clawpack package** [[1]](#1).
The riemann solver file is supposed to be uploaded into the PATH/clawpack/riemann/ and _init_.py, setup.py and static.py need to be updated to add this riemann solver. Just add the following lines to the corresponding .py files - 
(1) _init_.py
from . import Mu_I

(2) setup.py
one_d_riemann = ['Mu_I',
                  ...
                  
(3) static.py
num_eqn = {
        'Mu_I' : 1,
        'Full_model' : 1,
        ...
num_waves = {
        'Mu_I' : 3,
        'Full_model' : 3,
        ...

## References
<a id="1">[1]</a> 
Ketcheson, David I. and Mandli, Kyle T. and Ahmadia, Aron J. and Alghamdi, Amal and {Quezada de Luna}, Manuel and Parsani, Matteo and Knepley, Matthew G. and Emmett, Matthew (2012). 
PyClaw: Accessible, Extensible, Scalable Tools for Wave Propagation Problems. 
SIAM Journal on Scientific Computing,34(4), C210--C231.
