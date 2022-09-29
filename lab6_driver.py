from mypkg.IterationND import nd_iteration
import numpy as np

##### Part 1 #####
def evalF(x):
    F = np.zeros(2)
    
    F[0] = 4*x[0]**2 + x[1]**2 - 4
    F[1] = x[0] + x[1] - np.sin(x[0] - x[1])
    return F
    
def evalJ(x):
    J = np.array([[8*x[0], 2*x[1]], [1 - np.cos(x[0] - x[1]), 1 + np.cos(x[0] - x[1])]])
    return J
x0 = [1.5, 0]
Nmax = 100
tol = 1e-10

nd_iter = nd_iteration(evalF, evalJ)

[xstar,ier,its] =  nd_iter.LazyNewton(x0,tol,Nmax)
print(f"Approximated root is: {xstar}")
print(f'The error message reads: {ier}')
print(f'Number of iterations is: {its}')

[xstar,ier,its] =  nd_iter.SlackerNewton(x0,tol,Nmax)
print(f"Approximated root is: {xstar}")
print(f'The error message reads: {ier}')
print(f'Number of iterations is: {its}')
##### END #####