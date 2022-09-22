from mypkg.Iteration1D import Iteration1D
import numpy as np

f = lambda x : np.exp(x**2 + 7*x - 30) - 1
find = Iteration1D(f,'bisection')

# This covers everything we need for bisection, Newton's, and hybrid
find.fprime = lambda x : (2*x + 7)*np.exp(x**2 + 7*x - 30)
find.fprime2 = lambda x : (4*x**2 + 28*x + 51)*np.exp(x**2 + 7*x - 30)
find.a = 2
find.b = 4.5
find.tol = 1e-32
find.Nmax = 100
find.p0 = 4.5



x_bisection = find.root()
print(f"Bisection Sol: {x_bisection}")

find.method = 'newton'

x_newton = find.root()
print(f"Newton's Sol: {len(x_newton[0])}")

find.method = 'new_bisection'

x_merged = find.root()
print(f"Merged Sol: {x_merged}")
