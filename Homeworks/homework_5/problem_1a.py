# adding the package path to import mypkg classes
import sys
import os
sys.path.insert(1, os.path.abspath('.\..\..'))
from mypkg.IterationND import nd_iteration

import numpy as np
import matplotlib.pyplot as plt

##### Problem 1a #####
# How many max iterations
num_iters = 100
# Store all the iterations
all_x = np.zeros((2, num_iters))
# Scalar matrix to multiply by in each step
scalar_mat = [[1/6, 1/18], [0, 1/6]]

# Our tolerance
tol = 1e-15

x = np.transpose(np.array([1, 1]))
all_x[..., 0] = x
for i in range(1, num_iters):
    # Calculate f and g
    f = np.array([[3*x[0]**2 - x[1]**2], [3*x[0]*x[1]**2 - x[0]**3 - 1]])
    # Calculate the new x
    x = np.squeeze(np.expand_dims(x, -1) - np.matmul(scalar_mat, f))
    all_x[..., i] = x
    if np.linalg.norm(all_x[..., i-1] - x) < tol:
        all_x = all_x[..., :i+1]
        break

    

# final_x is the last x value in our convergence
final_x = x

print(f"Approximated root is: {final_x}")
print(f'Number of iterations is: {all_x.shape[1]}')

# Calculate convergence order
nd_iter = nd_iteration(None, None)
# Here we truncate the last 3 to avoid dividing by 0.
nd_iter.compute_order(np.swapaxes(all_x, 0, 1), final_x, 'problem_1a_convergence_plot.png')
##### END #####