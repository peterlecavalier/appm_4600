# adding the package path to import mypkg classes
import sys
import os
sys.path.insert(1, os.path.abspath('.\..\..'))
from mypkg.Iteration1D import Iteration1D
from mypkg.prini import prini

##### Problem 5a #####

# NEWTON'S
# specify function
f = lambda x: x**6 - x - 1
# specify derivative
fprime = lambda x: 6*x**5 - 1
# insantiate Iteration1D object
find = Iteration1D(f,'newton')
# set fprime and p0
find.fprime = fprime; find.p0 = 2
# set tol and Nmax
find.tol = 1e-10; find.Nmax = 100;
# find the root
x_newton = find.root()
# Print the root
p = prini("real", "The approximated root is:", x_newton[1])
p.print()
# Print the error code
print(f"The error code is: {x_newton[2]}")
# Print number of iterations
print(f"Newton's took {len(x_newton[0]) - 1} iterations.")

print("\n Newton's approximations:")

for i in x_newton[0][1:]:
    print(abs(i[0] - x_newton[0][-1])[0])

# SECANT
# specify function
f = lambda x: x**6 - x - 1
# insantiate Iteration1D object
find = Iteration1D(f,'secant')
# set p0 and p1
find.p0 = 2; find.p1 = 1
# set tol and Nmax
find.tol = 1e-15; find.Nmax = 100;
# find the root
x_secant = find.root()
# Print the root
p = prini("real", "The approximated root is:", x_secant[1])
p.print()
# Print the error code
print(f"The error code is: {x_secant[2]}")
# Print number of iterations
print(f"Secant took {len(x_secant[0]) - 2} iterations.")

print("\n Secant error:")

for i in x_secant[0][2:]:
    print(abs(i[0] - x_secant[0][-1])[0])
##### END #####

##### Problem 5b #####
# Kept these in the same file so we don't have to rerun root-finding
# First, Newton's:
find.compute_order(x_newton[0], x_newton[0][-1], 'problem_5b_plot_1.png')

# Second, secant:
find.compute_order(x_secant[0], x_secant[0][-1], 'problem_5b_plot_2.png')
##### END #####