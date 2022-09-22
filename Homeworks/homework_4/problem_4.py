# adding the package path to import mypkg classes
import sys
import os
sys.path.insert(1, os.path.abspath('.\..\..'))
from mypkg.Iteration1D import Iteration1D
from mypkg.prini import prini

import numpy as np
import matplotlib.pyplot as plt

##### Problem 4 #####

# PART i
# specify function
f = lambda x: np.exp(3*x) - 27*x**6 + 27*x**4*np.exp(x) - 9*x**2*np.exp(2*x)
# specify derivative
fprime = lambda x: 3*np.exp(3*x) - 162*x**5 + 108*x**3*np.exp(x) + 27*np.exp(x)*x**4 - 18*x*np.exp(2*x) - 18*x**2*np.exp(2*x)
# insantiate Iteration1D object
find = Iteration1D(f,'newton')
# set fprime and p0
find.fprime = fprime; find.p0 = 3.5
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

# compute the convergence rate and constant
# also plots and outputs lambda/alpha
[fit, diff1, diff2] = find.compute_order(x_newton[0], x_newton[1], 'problem_4_plot_1.png')

# PART ii
# Now, mu(x) = f(x)/f'(x). p_n+1 = p_n - mu(x)/mu'(x).
# After a lot of algebra, this simplifies to the following:
#mu = lambda x: (np.exp(3*x) - 27*x**6 + 27*x**4*np.exp(x) - 9*x**2*np.exp(2*x)) / (3*np.exp(3*x) - 162*x**5 + 108*x**3*np.exp(x) + 27*np.exp(x)*x**4 - 18*x*np.exp(2*x) - 18*x**2*np.exp(2*x))
mu = lambda x: (np.exp(x) - 3*x**2) / (3*(np.exp(x) - 6*x))
# Derivative (a big mess)
#muprime = lambda x: (4374*x**10 + 729*np.exp(x)*x**10 - 2916*np.exp(x)*x**9 - 4374*np.exp(x)*x**8 - 972*np.exp(2*x)*x**8 + 3888*np.exp(2*x)*x**7 + 972*np.exp(2*x)*x**6 + 486*np.exp(3*x)*x**6 - 1944*np.exp(3*x)*x**5 - 108*np.exp(4*x)*x**4 + 324*np.exp(3*x)*x**4 + 432*np.exp(4*x)*x**3 + 9*np.exp(5*x)*x**2 - 162*np.exp(4*x)*x**2 - 36*np.exp(5*x)*x + 18*np.exp(5*x)) / (3*np.exp(3*x) - 162*x**5 + 108*np.exp(x)*x**3 + 27*np.exp(x)*x**4 - 18*np.exp(2*x)*x - 18*np.exp(2*x)*x**2)**2
muprime = lambda x: ((x**2 - 4*x + 2)*np.exp(x) + 6*x**2) / (np.exp(x) - 6*x)**2
# insantiate Iteration1D object
find = Iteration1D(mu,'newton')
# set fprime and p0
find.fprime = muprime; find.p0 = 3.5
# set tol and Nmax
find.tol = 1e-13; find.Nmax = 100;
# find the root
x_newtons = find.root()
# Print the root
p = prini("real", "The approximated root is:", x_newtons[1])
p.print()
# Print the error code
print(f"The error code is: {x_newtons[2]}")
# Print number of iterations
print(f"Class Modified Newton's took {len(x_newtons[0])} iterations.")

print(x_newtons[0])

# compute the convergence rate and constant
# also plots and outputs lambda/alpha
[fit, diff1, diff2] = find.compute_order(x_newtons[0], x_newtons[1], 'problem_4_plot_2.png')


##### END #####