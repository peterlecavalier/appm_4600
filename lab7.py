from mypkg.IterationND import nd_iteration
import numpy as np
import matplotlib.pyplot as plt

def eval_lagrange(xeval,xint,yint,N):
    # xeval is one point we want to evaluate at
    # xint is interp nodes
    # yint is f(interp nodes)
    # N is degree of our polynomial
    lj = np.ones(N+1)
    
    for count in range(N+1):
       for jj in range(N+1):
           if (jj != count):
              lj[count] = lj[count]*(xeval - xint[jj])/(xint[count]-xint[jj])

    yeval = 0.
    
    for jj in range(N+1):
       yeval = yeval + yint[jj]*lj[jj]
  
    return(yeval)


f = lambda x: 1.0/(1+np.power(10.0*x, 2.0))

# degree of polynomial
n = 20
# interpolation nodes
x = np.linspace(-1, 1, n+1)
# f evaluated at the xs
fx = f(x)

# vandermonde construction
v = np.zeros((n+1, n+1))
for i in range(n+1):
    v[:, i] = x**i
# a is the poly coefficients
a = np.linalg.solve(v, fx)

# evaluating the polynomial at the points
neval = 1000
xeval = np.linspace(-1, 1, neval)
feval = f(xeval)
veval = np.zeros((neval, n+1))
for i in range(0, n+1):
    veval[:, i] = xeval**i

fig, ax = plt.subplots(1, 2)
ax[0].plot(xeval, feval, 'k-', )
ax[0].plot(xeval, veval@a, 'b--')
ax[1].semilogy(xeval, np.abs(veval@a - feval), 'r--')
ax[0].set_title('Actual plot vs. polynomial')
ax[1].set_title('Absolute error')
fig.suptitle('Monomial interpolation')
plt.savefig('monomial.png')
plt.show()

# Lagrange
lagrange_f = []
for i in xeval:
    lagrange_f.append(eval_lagrange(i, x, fx, n))

fig, ax = plt.subplots(1, 2)
ax[0].plot(xeval, lagrange_f, 'k-', )
ax[0].plot(xeval, veval@a, 'b--')
ax[1].semilogy(xeval, np.abs(veval@a - lagrange_f), 'r--')
ax[0].set_title('Actual plot vs. polynomial')
ax[1].set_title('Absolute error')
fig.suptitle('Lagrange interpolation')
plt.savefig('lagrange.png')
plt.show()
