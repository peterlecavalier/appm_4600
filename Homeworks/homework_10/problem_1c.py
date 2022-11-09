import sys
import os
sys.path.insert(1, os.path.abspath('.\..\..'))
from mypkg.Quadrature import MyQuad

from scipy.integrate import quad

##### Problem 1c #####
# Define function
f = lambda s: 1 / (1+s**2)

# Initialize our integration techniques
q = MyQuad(f, -5, 5)

# Run Composite Trapezoidal
comp_trap_result = q.composite_trap(n=1291)
# Run Composite Simpson's
comp_simp_result = q.composite_simp(n=108)
# Run SciPy built-in quadrature
scipy_result, scipy_err_est, info = quad(f, -5, 5, full_output=1)

# Compute errors
comp_trap_err = abs(comp_trap_result - scipy_result)
comp_simp_err = abs(comp_simp_result - scipy_result)

# Print the results of our integration techniques
print(f"Composite Trapezoidal={comp_trap_result}, with error={comp_trap_err}")
print(f"Composite Simpson's={comp_simp_result}, with error={comp_simp_err}")
print(f"SciPy built-in quadrature={scipy_result}")

# Print number of function evaluations for SciPy built-in quadrature
print(f"SciPy quad with 10^-6 tolerance evaluated f(s) {info['neval']} times.")
scipy_result, scipy_err_est, info = quad(f, -5, 5, epsabs=1e-4, full_output=1)
print(f"SciPy quad with 10^-4 tolerance evaluated f(s) {info['neval']} times.")
##### END #####