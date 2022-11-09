import numpy as np
from scipy.special import gamma

##### Problem 2c #####
xs = [2, 4, 6, 8, 10]

print("The following approximates the gamma function on the interval [0,inf] using Gauss-Laguerre quadrature:")

for x in xs:
    # Set up our function for this x
    f = lambda t: t**(x-1)

    # We can use x/2 here, since it will correctly integrate
    # polynomials of degree 2*deg-1 or less.
    pts, ws = np.polynomial.laguerre.laggauss(int(x / 2))

    fs = f(pts)
    
    result = np.sum(np.multiply(ws, fs))
    actual_gamma = gamma(x)
    rel_error = np.abs((result - actual_gamma) / actual_gamma)
    print(f"x={x}:")
    print(f"    ->Actual gamma = {actual_gamma}")
    print(f"    ->Approximated gamma = {result:.4}")
    print(f"    ->Relative error = {rel_error:.4}")
##### END #####