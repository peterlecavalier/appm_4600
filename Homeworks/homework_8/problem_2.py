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
    # n is the number of nodes
    # interpolation nodes
    # Instead of equally spacing, here we'll use Chebychev points
    xi = np.linspace(1, n, n)
    x = np.sort(5*np.cos((2*(xi) - 1)*np.pi / (2*(len(xi)))))
    # f evaluated at the xs
    fx = f(x)
    # fprime evaluated at the xs
    fprimex = fprime(x)

    # evaluating the polynomial at many points on the interval
    neval = 1000
    xeval = np.linspace(x[0], x[-1], neval)
    feval = f(xeval)

    # Setup the interpolation object
    interp = Interp1D(xint=x, yint=fx, yprimeint=fprimex, N=n-1)

    # Evaluate all interpolations
    lagrange_f, hermite_f, natcubic_f, clampcubic_f = [], [], [], []
    for i in xeval:
        lagrange_f.append(interp.lagrange(i))
        hermite_f.append(interp.hermite(i))
        natcubic_f.append(interp.new_natural_cubic_spline(i))
        clampcubic_f.append(interp.new_clamped_cubic_spline(i))

    # Plot everything
    lagrange_ax[idx][0].plot(xeval, lagrange_f, 'k-', label='interpolation')
    lagrange_ax[idx][0].plot(xeval, feval, 'b--', label='f(x)')
    lagrange_ax[idx][0].plot(x, fx, 'ro', label='nodes')
    lagrange_ax[idx][1].semilogy(xeval, np.abs(feval - lagrange_f), 'r')
    hermite_ax[idx][0].plot(xeval, hermite_f, 'k-', label='interpolation')
    hermite_ax[idx][0].plot(xeval, feval, 'b--', label='f(x)')
    hermite_ax[idx][0].plot(x, fx, 'ro', label='nodes')
    hermite_ax[idx][1].semilogy(xeval, np.abs(feval - hermite_f), 'r')
    natcubic_ax[idx][0].plot(xeval, natcubic_f, 'k-', label='interpolation')
    natcubic_ax[idx][0].plot(xeval, feval, 'b--', label='f(x)')
    natcubic_ax[idx][0].plot(x, fx, 'ro', label='nodes')
    natcubic_ax[idx][1].semilogy(xeval, np.abs(feval - natcubic_f), 'r')
    clampcubic_ax[idx][0].plot(xeval, clampcubic_f, 'k-', label='interpolation')
    clampcubic_ax[idx][0].plot(xeval, feval, 'b--', label='f(x)')
    clampcubic_ax[idx][0].plot(x, fx, 'ro', label='nodes')
    clampcubic_ax[idx][1].semilogy(xeval, np.abs(feval - clampcubic_f), 'r')

    for ax in axes:
        if idx == 0:
            ax[idx][0].set_title(f'Actual plot vs. polynomial - N = {n}')
            ax[idx][1].set_title('Absolute error')
            ax[idx][0].legend()
        else:
            ax[idx][0].set_title(f'N = {n}')

# Set figure titles
lagrange_fig.suptitle('Lagrange Interpolation (Chebychev nodes)')
hermite_fig.suptitle('Hermite Interpolation (Chebychev nodes)')
natcubic_fig.suptitle('Natural Cubic Spline Interpolation (Chebychev nodes)')
clampcubic_fig.suptitle('Clamped Cubic Spline Interpolation (Chebychev nodes)')

# Save each figure
lagrange_fig.savefig('2_lagrange.png')
hermite_fig.savefig('2_hermite.png')
natcubic_fig.savefig('2_natural_cubic_spline.png')
clampcubic_fig.savefig('2_clamped_cubic_spline.png')

# Show all figures
plt.show()
##### END #####