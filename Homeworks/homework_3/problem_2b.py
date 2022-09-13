# adding the package path to import mypkg classes
import sys
import os
sys.path.insert(1, os.path.abspath('.\..\..'))
from mypkg.Iteration1D import Iteration1D
from mypkg.prini import prini

import numpy as np

##### Problem 2b #####
# specify function
# This is the expanded version of (x-5)^9
f = lambda x: x**9 - 45*x**8 + 900*x**7 - 10500*x**6 + 78750*x**5 - 393750*x**4 + 1312500*x**3 - 2812500*x**2 + 3515625*x - 1953125
# insantiate Iteration1D object
find = Iteration1D(f,'bisection')
# set initial interval
find.a = 4.82; find.b = 5.2
# set tol and Nmax
find.tol = 1e-16; find.Nmax = 100
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