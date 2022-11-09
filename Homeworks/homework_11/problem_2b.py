import sys
import os
sys.path.insert(1, os.path.abspath('.\..\..'))
from mypkg.Quadrature import MyQuad

from scipy.integrate import quad
import numpy as np
from scipy.special import gamma

##### Problem 2b #####
# The interval I chose was [0,42]
xs = [2, 4, 6, 8, 10]

# Defining this just for comparison to the scipy approx
q = MyQuad(a=0, b=42, n=1721)

print("The following approximates the gamma function on the interval [0,42]")
print("using SciPy's built-in quadrature:")

for x in xs:
    f = lambda t: t**(x-1) * np.exp(-t)
    cur_trap_result = q.composite_trap(f=f)
    cur_result, _, info = quad(f, 0, 42, full_output=1)
    actual_gamma = gamma(x)
    rel_error = np.abs((cur_result - actual_gamma) / actual_gamma)
    trap_rel_error = np.abs((cur_trap_result - actual_gamma) / actual_gamma)
    print(f"The relative error for x={x} is {rel_error:.4}, while Composite Trapezoidal is {trap_rel_error:.4}")
    print(f"    -> Number of function evaluations for SciPy quad is {info['neval']} (vs 1722 for composite trap)")
##### END #####