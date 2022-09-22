# adding the package path to import mypkg classes
import sys
import os
sys.path.insert(1, os.path.abspath('.\..\..'))
from mypkg.Iteration1D import Iteration1D
from mypkg.prini import prini

from scipy.special import erf
import numpy as np

##### Problem 1c #####
# specify function
f = lambda x: 35*erf(x/(2*np.sqrt(0.715392))) - 15
# specify derivative
fprime = lambda x: (35/np.sqrt(0.715392*np.pi))*np.exp(-x**2/2.861568)
# insantiate Iteration1D object
find = Iteration1D(f,'newton')
# set fprime and p0
find.fprime = fprime; find.p0 = 0.01
# set tol and Nmax
find.tol = 1e-13; find.Nmax = 100;
# find the root
x_newton = find.root()
# Print the root
p = prini("real", "The approximated root is:", x_newton[1])
p.print()
# Print the error code
print(f"The error code is: {x_newton[2]}")
# Print number of iterations
print(f"Newton's took {len(x_newton[0])} iterations.")
##### END #####