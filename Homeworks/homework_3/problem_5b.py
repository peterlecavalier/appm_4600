# adding the package path to import mypkg classes
import sys
import os
sys.path.insert(1, os.path.abspath('.\..\..'))
from mypkg.Iteration1D import Iteration1D
from mypkg.my2DPlot import my2DPlot as myplt
from mypkg.prini import prini

import numpy as np

##### Problem 5b #####
# specify function
f = lambda x: -np.sin(2*x) + (5*x)/4 - (3/4)
# insantiate Iteration1D object
find = Iteration1D(f,'fixedpt_mod')
# set tol and Nmax
find.tol = 0.4e-10; find.Nmax = 100

# set initial guess
find.p0 = -0.8984
# find the root
x_fixed_pt = find.root()
# Print the root
if x_fixed_pt[1] == 0:
    p = prini("real", "The approximated root is:", x_fixed_pt[0])
    p.print()
else:
    # Print the error code
    p = prini("real", "There was an error: code", x_fixed_pt[1])
    p.print()
# ROOT IS NOT FOUND

find.p0 = -0.5444
x_fixed_pt = find.root()
if x_fixed_pt[1] == 0:
    p = prini("real", "The approximated root is:", x_fixed_pt[0])
    p.print()
else:
    p = prini("real", "There was an error: code", x_fixed_pt[1])
    p.print()
# ROOT IS FOUND

find.p0 = 1.7321
x_fixed_pt = find.root()
if x_fixed_pt[1] == 0:
    p = prini("real", "The approximated root is:", x_fixed_pt[0])
    p.print()
else:
    p = prini("real", "There was an error: code", x_fixed_pt[1])
    p.print()
# ROOT IS NOT FOUND

find.p0 = 3.1618
x_fixed_pt = find.root()
if x_fixed_pt[1] == 0:
    p = prini("real", "The approximated root is:", x_fixed_pt[0])
    p.print()
else:
    p = prini("real", "There was an error: code", x_fixed_pt[1])
    p.print()
# ROOT IS FOUND

find.p0 = 4.518
x_fixed_pt = find.root()
if x_fixed_pt[1] == 0:
    p = prini("real", "The approximated root is:", x_fixed_pt[0])
    p.print()
else:
    p = prini("real", "There was an error: code", x_fixed_pt[1])
    p.print()
# ROOT IS NOT FOUND

# Plotting the original function:
plt = myplt(lambda x: x - 4*np.sin(2*x) - 3,-1,6, 'Original f(x)')
plt.labels('x','y')
# Plotting |f'(x)|:
plt.addPlot(lambda x: np.abs(-2*np.cos(2*x) + (5/4)), "Absolute f'(x)")
# Adding y = 1 to the plot
plt.addPlot(lambda x: np.zeros(x.shape) + 1)
plt.dotted()
plt.color('black')
# Adding y = 0 to the plot
plt.addPlot(lambda x: 0*x)
plt.dotted()
plt.color('black')
plt.legend()
plt.save('5b_plot.png')
plt.show()
##### END #####