# adding the package path to import mypkg classes
import sys
import os
sys.path.insert(1, os.path.abspath('.\..\..'))
from mypkg.IterationND import nd_iteration
import numpy as np

##### Problem 3b #####
# Our initial guess
x0 = np.array([[1], [1], [1]])

# define our function
f = lambda x: x[0]**2 + 4*x[1]**2 + 4*x[2]**2 - 16

Nmax = 100
tol = 1e-32

all_x = np.zeros((Nmax + 1, 3, 1))

xstar = None

for i in range(Nmax):
    all_x[i] = x0

    # Calculate partial derivatives
    fx = 2*x0[0]
    fy = 8*x0[1]
    fz = 8*x0[2]

    # Evaluate f
    f_eval = f(x0)
    
    # Evaluate d
    d = f_eval / (fx**2 + fy**2 + fz**2)

    # Vector to house all dfs
    df = np.array([fx*d,fy*d,fz*d])

    # Calculate the iterate
    x1 = x0 - df

    # Check if it's within the tol
    if np.linalg.norm(x1 - x0) < tol:
        xstar = x1
        all_x[i+1] = x1
        all_x = all_x[:i+2]
        break
    else:
        x0 = x1

if xstar is None:
    print("Iteration never converged.")
else:
    print(f"Iteration converged to the value: \n{xstar}\nin {all_x.shape[0] - 1} iterations.")

    # Calculate convergence order
    nd_iter = nd_iteration(None, None)
    # Here we truncate the last 3 to avoid dividing by 0.
    nd_iter.compute_order(all_x[:-3], xstar, fig_fp='problem_3b_plot.png')

##### END #####