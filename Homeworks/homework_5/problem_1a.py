import numpy as np
import matplotlib.pyplot as plt

# How many iterations
num_iters = 60
# Store all the iterations
all_x = np.zeros((2, num_iters))
# Scalar matrix to multiply by in each step
scalar_mat = [[1/6, 1/18], [0, 1/6]]


x = np.transpose(np.array([1, 1]))
all_x[..., 0] = x
for i in range(1, num_iters):
    # Calculate f and g
    f = np.array([[3*x[0]**2 - x[1]**2], [3*x[0]*x[1]**2 - x[0]**3 - 1]])
    # Calculate the new x
    x = np.squeeze(np.expand_dims(x, -1) - np.matmul(scalar_mat, f))
    all_x[..., i] = x

# final_x is the last x value in our convergence
final_x = x

# Plot absolute error
fig, ax = plt.subplots(1, 2, figsize=(10, 5))
ax[0].scatter(np.arange(0, num_iters, 1), np.abs(all_x[0] - final_x[0]))
ax[1].scatter(np.arange(0, num_iters, 1), np.abs(all_x[1] - final_x[1]), c='orange')
ax[0].set_yscale('log')
ax[1].set_yscale('log')

ax[0].set_xlabel('Iteration Number')
ax[1].set_xlabel('Iteration Number')
ax[0].set_ylabel('Absolute Error x')
ax[1].set_ylabel('Absolute Error y')
plt.suptitle('Convergence of the Iteration Scheme')
plt.savefig('problem_1a_plot.png')
plt.show()