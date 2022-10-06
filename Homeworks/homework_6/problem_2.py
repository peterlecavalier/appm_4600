import sys
import os
sys.path.insert(1, os.path.abspath('.\..\..'))
from mypkg.IterationND import nd_iteration
import numpy as np

##### Problem 2 #####
def evalF(x_inp):
    x = np.squeeze(np.copy(x_inp))
    F = np.array([[x[0] + np.cos(x[0]*x[1]*x[2]) - 1],
     [(1 - x[0])**(1/4) + x[1] + 0.05*x[2]**2 - 0.15*x[2] - 1],
     [-x[0]**2 - 0.1*x[1]**2 + 0.01*x[1] + x[2] - 1]])
    return F
    
def evalJ(x_inp):
    x = np.squeeze(np.copy(x_inp))
    J = np.array([[1 - x[1]*x[2]*np.sin(x[0]*x[1]*x[2]), -x[0]*x[2]*np.sin(x[0]*x[1]*x[2]), -x[0]*x[1]*np.sin(x[0]*x[1]*x[2])],
     [-(1/4)*(1-x[0])**(-3/4), 1, 0.1*x[2] - 0.15],
     [-2*x[0], -0.2*x[1] + 0.01, 1]])
    return J

Nmax = 200
tol = 1e-6

# Set up our iteration object
nd_iter = nd_iteration(evalF, evalJ)

# Narrowed down this initial guess using guess/check
x0 = np.array([[0],[-2],[-2]])

# Newton's Method
[allx_newton, xstar_newton,ier_newton,its_newton] =  nd_iter.Newton(np.squeeze(x0),tol,Nmax)
print(f"Newton's approximated root is: {xstar_newton}")
print(f"Newton's number of iterations is: {its_newton}\n")

# Steepest Descent Method
[allx_steep, xstar_steep, ier_steep, its_steep] = nd_iter.SteepestDescent(x0, tol, Nmax)
print(f"Steepest Descent's approximated root is: {np.squeeze(xstar_steep)}")
print(f"Steepest Descent's number of iterations is: {its_steep}")
print("")

# Hybrid Method
init_tol = 5e-2
[allx_new_steep, xstar_new_steep, ier_steep, its_steep] = nd_iter.SteepestDescent(x0, init_tol, Nmax)
new_x0 = np.squeeze(xstar_new_steep)
[allx_new_newton, xstar_new_newton, ier_newton, its_newton] =  nd_iter.Newton(new_x0,tol,Nmax)
print(f"Hybrid's approximated root is: {xstar_new_newton}")
print(f"Hybrid's number of iterations is: {its_newton + its_steep} ({its_steep} for steepest, {its_newton} for Newton)")

# Compute the order of each equation
nd_iter.compute_order(allx_newton, xstar_newton, fig_fp="problem_2_newtons_plot.png", print_info=False)
nd_iter.compute_order(allx_steep, xstar_steep, fig_fp="problem_2_steepest_plot.png", print_info=False)
allx_hybrid = np.concatenate((np.squeeze(allx_new_steep), allx_new_newton), axis=0)
nd_iter.compute_order(allx_hybrid, xstar_new_newton, fig_fp="problem_2_hybrid_plot.png", print_info=False)
##### END #####