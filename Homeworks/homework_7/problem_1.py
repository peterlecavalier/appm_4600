import sys
import os
sys.path.insert(1, os.path.abspath('.\..\..'))
import numpy as np
import matplotlib.pyplot as plt

##### Problem 1 #####
def vandermonde_mat(xs):
    """
    xs is a [n,1] array of x values
    """
    V = np.zeros((len(xs), len(xs)))
    for i in range(len(xs)):
        V[:, i] = np.power(xs, i)
    return V

f = lambda x: 1/(1+(10*x)**2)

# Finer grain points for plotting f(x)
x_fine = np.linspace(-1, 1, 1001)
fx_fine = f(x_fine)

fig, ax = plt.subplots(6, 3, figsize=(16, 10), sharex=True)

for n in range(2, 20):
    # All xi points
    xi = np.linspace(1, n, n)
    h = 2/(len(xi) - 1)
    xi = -1 + (xi - 1)*h
    fxi = f(xi)
    v = vandermonde_mat(xi)
    # Solve for coefficients
    c = np.linalg.solve(v, np.expand_dims(fxi, -1))
    # Evaluate coefficients at fine points
    px = np.polyval(np.flip(c.flatten()), x_fine)

    # Calculate where to put the plot
    row = int(np.floor((n-2)/3))
    col = (n-2) % 3

    # Plot everything
    ax[row][col].plot(xi, fxi, 'o')
    ax[row][col].plot(x_fine, px, label=f'p(x), N={n}')
    ax[row][col].plot(x_fine, fx_fine, label='f(x)')
    ax[row][col].legend()
plt.savefig('problem1_plot.png')
plt.show()
##### END #####