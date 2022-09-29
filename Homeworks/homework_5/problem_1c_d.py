# adding the package path to import mypkg classes
import sys
import os
sys.path.insert(1, os.path.abspath('.\..\..'))
from mypkg.IterationND import nd_iteration
import numpy as np
import matplotlib.pyplot as plt


##### Problem 1c #####
def evalF(x):
    F = np.zeros(2)
    
    F[0] = 3*x[0]**2 - x[1]**2
    F[1] = 3*x[0]*x[1]**2 - x[1]**3 - 1
    return F
    
def evalJ(x):
    J = np.array([[6*x[0], -2*x[1]], [3*x[1]**2 - 3*x[0]**2, 6*x[0]*x[1]]])
    return J
x0 = [1, 1]
Nmax = 100
tol = 1e-32

nd_iter = nd_iteration(evalF, evalJ)

[xstar,ier,its] =  nd_iter.Newton(x0,tol,Nmax)
print(f"Approximated root is: {xstar}")
tmp = print(f'The error message reads: {ier}')
tmp = print(f'Number of iterations is: {its}')
##### END #####

##### Problem 1d #####
# Since we already have the result from the code above,
# it's a lot easier to just analytically evaluate here.

# This F evaluation should return approximately [0, 0] to verify
# our solution is accurate.
sol = evalF(xstar)
print(f"The analytical evaluation of our solution is: {sol}")
##### END #####