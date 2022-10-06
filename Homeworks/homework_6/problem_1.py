import sys
import os
sys.path.insert(1, os.path.abspath('.\..\..'))
from mypkg.IterationND import nd_iteration
import numpy as np

# This allows us to catch warnings
import warnings
warnings.filterwarnings("error")

##### Problem 1 #####
def evalF(x):
    F = np.array([[x[0]**2 + x[1]**2 - 4], [np.exp(x[0]) + x[1] - 1]])
    return F
    
def evalJ(x):
    J = np.array([[2*x[0], 2*x[1]], [np.exp(x[0]), 1]])
    return J

Nmax = 200
tol = 1e-32

# Set up our iteration object
nd_iter = nd_iteration(evalF, evalJ)

def run_methods(nd_iter, x0):
    try:
        [_, xstar_newton,ier_newton,its_newton] =  nd_iter.Newton(x0,tol,Nmax)
        print(f"Newton's approximated root is: {xstar_newton}")
        print(f"Newton's number of iterations is: {its_newton}\n")
    except np.linalg.LinAlgError:
        print("ERROR: Newton's Method attempted an inverse on a singular matrix\n")
    except RuntimeWarning:
        print("ERROR: Newton's Method encountered runtime error (likely over/underflow)\n")
    
    try:
        [xstar_lazy, ier_lazy, its_lazy] =  nd_iter.LazyNewton(x0,tol,Nmax)
        print(f"Lazy Newton's approximated root is: {xstar_lazy}")
        print(f"Lazy Newton's number of iterations is: {its_lazy}\n")
    except np.linalg.LinAlgError:
        print("ERROR: Lazy Newton's Method attempted an inverse on a singular matrix\n")
    except RuntimeWarning:
        print("ERROR: Lazy Newton's Method encountered runtime error (likely over/underflow)\n")

    try:
        [xstar_broyden, ier_broyden, its_broyden] =  nd_iter.Broyden(x0,tol,Nmax)
        print(f"Broyden's approximated root is: {xstar_broyden}")
        print(f"Broyden's number of iterations is: {its_broyden}\n")
    except np.linalg.LinAlgError:
        print("ERROR: Broyden's Method attempted an inverse on a singular matrix\n")
    except RuntimeWarning:
        print("ERROR: Broyden's Method encountered runtime error (likely over/underflow)\n")

# First, we do part i
x0 = np.array([1, 1])
print('Part (i) => (1, 1)')
print('------------------')
run_methods(nd_iter, x0)


# Next, part ii
x0 = np.array([1, -1])
print('Part (ii) => (1, -1)')
print('------------------')
run_methods(nd_iter, x0)

# Finally, part iii
x0 = np.array([0, 0])
print('Part (iii) => (0, 0)')
print('------------------')
run_methods(nd_iter, x0)
##### END #####