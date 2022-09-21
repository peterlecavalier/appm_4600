# adding the package path to import mypkg classes
import sys
import os
sys.path.insert(1, os.path.abspath('.\..\..'))
from mypkg.Iteration1D import Iteration1D
from mypkg.prini import prini

from scipy.special import erf
import numpy as np

##### Problem 1b #####
# specify function
f = lambda x: 35*erf(x/(2*np.sqrt(0.715392))) - 15
# insantiate Iteration1D object
find = Iteration1D(f,'bisection')
# set initial interval
find.a = 0; find.b = 1
# set tol and Nmax
find.tol = 1e-13; find.Nmax = 100
# find the root
x_bisection = find.root()
# Print the root
p = prini("real", "The approximated root is:", x_bisection[0])
p.print()
# Print the error code
print(f"The error code is: {x_bisection[1]}")
# Print the number of iterations
print(f"Number of iterations to reach approximation: {x_bisection[2]}")
##### END #####