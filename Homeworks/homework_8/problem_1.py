import sys
import os
sys.path.insert(1, os.path.abspath('.\..\..'))
from mypkg.Interp1D import Interp1D

import numpy as np
import matplotlib.pyplot as plt

##### Problem 1 #####
interp = Interp1D()

# Function and derivative to evaluate
f = lambda x: 1.0/(1+np.power(x, 2.0))
fprime = lambda x: (-2*x)/((1+x**2)**2)

# Set up each separate figure
lagrange_fig, lagrange_ax = plt.subplots(4, 2, figsize=(8, 10))
plt.subplots_adjust(hspace=0.5)
hermite_fig, hermite_ax = plt.subplots(4, 2, figsize=(8, 10))
plt.subplots_adjust(hspace=0.6)
natcubic_fig, natcubic_ax = plt.subplots(4, 2, figsize=(8, 10))
plt.subplots_adjust(hspace=0.8)
clampcubic_fig, clampcubic_ax = plt.subplots(4, 2, figsize=(8, 10))
plt.subplots_adjust(hspace=1.0)

axes = [lagrange_ax, hermite_ax, natcubic_ax, clampcubic_ax]

for idx, n in enumerate([5, 10, 15, 20]):
    # n is the degree of the polynomial
    # interpolation nodes
    x = np.linspace(-5, 5, n+1)
    # f evaluated at the xs
    fx = f(x)
    # fprime evaluated at the xs
    fprimex = fprime(x)

    # evaluating the polynomial at many points on the interval
    neval = 1000
    xeval = np.linspace(-5, 5, neval)
    feval = f(xeval)

    # Setup the interpolation object
    interp = Interp1D(xint=x, yint=fx, yprimeint=fprimex, N=n)

    # Evaluate all interpolations
    lagrange_f, hermite_f, natcubic_f, clampcubic_f = [], [], [], []
    for i in xeval:
        lagrange_f.append(interp.lagrange(i))
        hermite_f.append(interp.hermite(i))

    # Plot everything
    lagrange_ax[idx][0].plot(xeval, lagrange_f, 'k-')
    lagrange_ax[idx][0].plot(xeval, feval, 'b--')
    hermite_ax[idx][0].plot(xeval, hermite_f, 'k-')
    hermite_ax[idx][0].plot(xeval, feval, 'b--')
    lagrange_ax[idx][1].semilogy(xeval, np.abs(feval - lagrange_f), 'r--')
    hermite_ax[idx][1].semilogy(xeval, np.abs(feval - hermite_f), 'r--')

    for ax in axes:
        if idx == 0:
            ax[idx][0].set_title(f'Actual plot vs. polynomial - N = {n}')
            ax[idx][1].set_title('Absolute error')
        else:
            ax[idx][0].set_title(f'N = {n}')

# Set figure titles
lagrange_fig.suptitle('Lagrange Interpolation')
hermite_fig.suptitle('Hermite Interpolation')
natcubic_fig.suptitle('Natural Cubic Spline Interpolation')
clampcubic_fig.suptitle('Clamped Cubic Spline Interpolation')

# Save each figure
lagrange_fig.savefig('lagrange.png')
hermite_fig.savefig('hermite.png')
natcubic_fig.savefig('natural_cubic_spline.png')
clampcubic_fig.savefig('clamped_cubic_spline.png')

# Show all figures
plt.show()
##### END #####