from mypkg.Iteration1D import Iteration1D

# specify function
f = lambda x: (x-2)**3
# insantiate Iteration1D object using bisection
# this should set self.f = f, self.method = 'bisection'
# and other all attributes, like self.tol = None
find = Iteration1D(f,'bisection')
# set initial interval
find.a = 1.5; find.b = 2.5
# set tol and Nmax
find.tol = 1e-6; find.Nmax = 100
# find the root
x_bisection = find.root()
print(x_bisection)

# now try another method
find.method = 'fixedpt'
# recast problem
find.f = lambda x: (8-12*x)/(x**2-6*x)
# set initial guess
find.p0 = 1.2
# find the root
x_fixedpt = find.root()
print(x_fixedpt)