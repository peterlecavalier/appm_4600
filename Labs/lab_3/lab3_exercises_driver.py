from mypkg.Iteration1D import Iteration1D
f = lambda x : x**2 * (x-1)
find = Iteration1D(f,'bisection')
# set initial interval
find.a = 0.5; find.b = 2
# set tol and Nmax
find.tol = 1e-6; find.Nmax = 100
# find the root
x_bisection = find.root()
print(f"1a: {x_bisection}")

find.a = -1; find.b = 0.5
x_bisection = find.root()
print(f"1b: {x_bisection}")

find.a = -1; find.b = 2
x_bisection = find.root()
print(f"1c: {x_bisection}")

find.method = 'fixedpt'
# recast problem
find.f = lambda x: x**(3/2)
# set initial guess
find.p0 = 0
# find the root
x_fixedpt = find.root()
print(f"2, root 1: {x_fixedpt}")

# set initial guess
find.p0 = 1
# find the root
x_fixedpt = find.root()
print(f"2, root 1: {x_fixedpt}")
