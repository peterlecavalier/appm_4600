# adding the package path to import mypkg classes
import sys
import os
sys.path.insert(1, os.path.abspath('.\..\..'))
from mypkg.Iteration1D import Iteration1D
from mypkg.prini import prini

import numpy as np

##### Problem 3b #####
# specify function
f = lambda x: x**3 + x - 4
# insantiate Iteration1D object
find = Iteration1D(f,'bisection')
# set initial interval
find.a = 1; find.b = 4
# set tol and Nmax
find.tol = 1e-3; find.Nmax = 100
# find the root
x_bisection = find.root()
# Print the root
p = prini("real", "The approximated root is:", x_bisection[0])
p.print()
# Print the error code
p = prini("real", "The error code is:", x_bisection[1])
p.print()
# Print the number of iterations
p = prini("real", "Number of iterations to reach approximation:", x_bisection[2])
p.print()
##### END #####