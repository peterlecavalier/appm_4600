import numpy as np

class Iteration1D:
    def __init__(self, f, method):
        # we assign to or initialize as None all self attributes
        self.f = f
        self.method = method
        # initial interval for bisection
        self.a = None
        self.b = None
        # initial guess for newton/fixedpt
        self.p0 = None
        # tolerance and max iter
        self.tol = None
        self.Nmax = None
        # info message
        self.info = None
        # root
        self.pstar = None
        # iters for newton or fixedpt
        self.p_iters = None

    def root(self):
        if self.method == 'bisection':
            if self.f is not None and self.a is not None and self.b is not None and self.tol is not None and self.Nmax is not None:
                pstar = bisection(self.f, self.a, self.b, self.tol, self.Nmax)
            else:
                return -1
        elif self.method == 'fixedpt':
            if self.f is not None and self.p0 is not None and self.tol is not None and self.Nmax is not None:
                pstar = fixedpt(self.f, self.p0, self.tol, self.Nmax)
            else:
                return -1
        elif self.method == 'fixedpt_mod':
            if self.f is not None and self.p0 is not None and self.tol is not None and self.Nmax is not None:
                pstar = fixedpt_mod(self.f, self.p0, self.tol, self.Nmax)
            else:
                return -1
        elif self.method == 'fixedpt_mod2':
            if self.f is not None and self.p0 is not None and self.tol is not None and self.Nmax is not None:
                pstar = fixedpt_mod2(self.f, self.p0, self.tol, self.Nmax)
            else:
                return -1
        
        return pstar # the root


def bisection(f,a,b,tol,Nmax):
    '''
    Inputs:
      f,a,b       - function and endpoints of initial interval
      tol, Nmax   - bisection stops when interval length < tol
                  - or if Nmax iterations have occured
    Returns:
      astar - approximation of root
      ier   - error message
            - ier = 1 => cannot tell if there is a root in the interval
            - ier = 0 == success
            - ier = 2 => ran out of iterations
            - ier = 3 => other error ==== You can explain
      count - number of iterations it took to find the root (or 0 if error)
    '''

    '''     first verify there is a root we can find in the interval '''
    fa = f(a); fb = f(b);
    if (fa*fb>0):
       ier = 1
       astar = a
       return [astar, ier, 0]

    ''' verify end point is not a root '''
    if (fa == 0):
      astar = a
      ier = 0
      return [astar, ier, 0]

    if (fb ==0):
      astar = b
      ier = 0
      return [astar, ier, 0]

    count = 0
    while (count < Nmax):
      count = count + 1
      c = 0.5*(a+b)
      fc = f(c)

      if (fc ==0):
        astar = c
        ier = 0
        return [astar, ier, count]

      if (fa*fc<0):
         b = c
      elif (fb*fc<0):
        a = c
        fa = fc
      else:
        astar = c
        ier = 3
        return [astar, ier, count]

      if (abs(b-a)<tol):
        astar = a
        ier =0
        return [astar, ier, count]

    astar = a
    ier = 2
    return [astar,ier, count] 

def fixedpt(f,x0,tol,Nmax):

    ''' x0 = initial guess''' 
    ''' Nmax = max number of iterations'''
    ''' tol = stopping tolerance'''

    count = 0
    while (count <Nmax):
       count = count +1
       x1 = f(x0)
       if (abs(x1-x0) <tol):
          xstar = x1
          ier = 0
          return [xstar,ier]
       x0 = x1

    xstar = x1
    ier = 1
    return [xstar, ier]


def fixedpt_mod(f,x0,tol,Nmax):
    """
    Modified fixed pt method to return approximation within
    ABSOLUTE error tolerance instead of relative.
    """

    ''' x0 = initial guess''' 
    ''' Nmax = max number of iterations'''
    ''' tol = stopping tolerance'''

    count = 0
    while (count <Nmax):
       count = count +1
       x1 = f(x0)
       if (abs(x1-x0)/abs(x0) <tol):
          xstar = x1
          ier = 0
          return [xstar,ier]
       x0 = x1

    xstar = x1
    ier = 1
    return [xstar, ier]

def fixedpt_mod2(f,x0,tol,Nmax):
    """
    Modified fixed pt method to track all iteration approximations
    """

    ''' x0 = initial guess''' 
    ''' Nmax = max number of iterations'''
    ''' tol = stopping tolerance'''

    all_iters = np.zeros((Nmax, 1))

    count = 0
    while (count <Nmax):
       count = count +1
       x1 = f(x0)
       x[count - 1] = x1
       if (abs(x1-x0) <tol):
          xstar = x1
          ier = 0
          return [all_iters, xstar, ier]
       x0 = x1

    xstar = x1
    ier = 1
    return [all_iters, xstar, ier]