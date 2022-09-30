# adding the package path to import mypkg classes
import sys
import os
sys.path.insert(1, os.path.abspath('.\..\..'))
from mypkg.IterationND import nd_iteration
import numpy as np
import matplotlib.pyplot as plt


##### Problem 1c #####
def evalF(x):
    F = np.array([[3 * np.power(x[0], 2) - np.power(x[1], 2)], [3 * x[0] * np.power(x[1], 2) - np.power(x[0], 3) - 1]])
    return F
    
def evalJ(x):
    J = np.array([[6*x[0], -2*x[1]], [3*x[1]**2 - 3*x[0]**2, 6*x[0]*x[1]]])
    return J
x0 = np.array([1.0, 1.0])
Nmax = 200
tol = 1e-15

nd_iter = nd_iteration(evalF, evalJ)

[all_x, xstar,ier,its] =  nd_iter.Newton(x0,tol,Nmax)
print(f"Approximated root is: {xstar}")
print(f'The error message reads: {ier}')
print(f'Number of iterations is: {its}')

# Compute order of convergence
nd_iter.compute_order(all_x[:-1], xstar, "problem_1c_plot.png")
##### END #####

##### Problem 1d #####
# Since we already have the result from the code above,
# it's a lot easier to just analytically evaluate here.

# This F evaluation should return approximately [0, 0] to verify
# our solution is accurate.
sol = evalF(xstar)
print(f"The analytical evaluation of our solution is: {sol}")
##### END #####