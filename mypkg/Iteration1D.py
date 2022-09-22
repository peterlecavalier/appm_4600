import numpy as np
import matplotlib.pyplot as plt

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
        self.p1 = None
        # tolerance and max iter
        self.tol = None
        self.Nmax = None
        # info message
        self.info = None
        # root
        self.pstar = None
        # iters for newton or fixedpt
        self.p_iters = None
        # derivative for newton
        self.fprime = None
        self.fprime2 = None

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
      elif self.method == 'newton':
        if self.f is not None and self.fprime is not None and self.p0 is not None and self.tol is not None and self.Nmax is not None:
          pstar = newtons(self.f, self.fprime, self.p0, self.tol, self.Nmax)
        else:
          return -1
      elif self.method == 'secant':
        if self.f is not None and self.p0 is not None and self.p1 is not None and self.tol is not None and self.Nmax is not None:
          pstar = secant(self.f, self.p0, self.p1, self.tol, self.Nmax)
        else:
          return -1
      elif self.method == 'new_bisection':
        if self.f is not None and self.fprime is not None and self.fprime2 is not None and self.a is not None and self.b is not None and self.tol is not None and self.Nmax is not None:
          pstar = new_bisection(self.f,self.fprime,self.fprime2,self.a,self.b,self.tol,self.Nmax)
        else:
          return -1
        
      return pstar # the root
    
    def compute_order(self, x, xstar, fig_fp=None):
      diff1 = np.abs(x[1:]-xstar)
      # p_n-p (from the first index to the second to last)
      diff2 = np.abs(x[0:-1]-xstar)
      # linear fit to log of differences

      # Following avoids dividing by zero:
      while(True):
        if diff1[-1] == 0 or diff2[-1] == 0:
          diff1 = diff1[:-1]
          diff2 = diff2[:-1]
        else:
          break

      fit = np.polyfit(np.log(diff2.flatten()),np.log(diff1.flatten()),1)
      print('the order equation is')
      print('log(|p_{n+1}-p|) = log(lambda) + alpha*log(|p_n-p|) where')
      print('lambda = ' + str(np.exp(fit[1])))
      print('alpha = ' + str(fit[0]))

      # plot the data
      plt.loglog(diff2,diff1,'ro',label='fixedpt data')
      # plot the fit
      plt.loglog(diff2,np.exp(fit[1]+fit[0]*np.log(diff2)),'b-',label='fixedpt fit')
      # label the plot
      plt.xlabel('$|p_{n}-p|$')
      plt.ylabel('$|p_{n+1}-p|$')
      plt.legend()
      plt.title(f'lambda = {np.exp(fit[1])}, alpha {fit[0]}')
      if fig_fp is not None:
        plt.savefig(fig_fp)
      plt.show()
      return [fit,diff1,diff2]



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
    all_iters[0] = x0
    count = 0
    for count in range(1, Nmax):
       count = count +1
       x1 = f(x0)
       all_iters[count] = x1
       if (abs(x1-x0) <tol):
          xstar = x1
          ier = 0
          all_iters = all_iters[:count+1]
          return [all_iters, xstar, ier]
       x0 = x1

    xstar = x1
    ier = 1
    return [all_iters, xstar, ier]

def newtons(f, fprime, p0, tol, Nmax):
  all_iters = np.zeros((Nmax, 1))
  all_iters[0] = p0
  for j in range(1, Nmax):
    p = p0 - f(p0)/fprime(p0)
    all_iters[j] = p
    if abs(p - p0) < tol:
      pstar = p
      ier = 0
      all_iters = all_iters[:j+1]
      return [all_iters, pstar, ier]
    
    p0 = p

  ier = 1
  pstar = p
  return [all_iters, pstar, ier]

def secant(f, x0, x1, tol, Nmax):
  all_iters = np.zeros((Nmax, 1))
  all_iters[0] = x0
  all_iters[1] = x1
  for j in range(2, Nmax):
    x = x1 - f(x1)*(x0 - x1)/(f(x0) - f(x1))
    all_iters[j] = x
    if abs(x - x1) < tol:
      xstar = x
      ier = 0
      all_iters = all_iters[:j+1]
      return [all_iters, xstar, ier]
    
    x0 = x1
    x1 = x

  ier = 1
  xstar = x
  return [all_iters, xstar, ier]

def new_bisection(f,fprime,fprime2,a,b,tol,Nmax):
    '''
    Inputs:
      f,fprime, fprime2, a,b       - function and endpoints of initial interval
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

      # Instead of checking within the tolerance, check if it's within the basin of convergence
      abs_gprime = f(a)*fprime2(a)/(fprime(a)**2)
      abs_gprime = abs(abs_gprime)
      if (abs_gprime < 1):
        # Calling Newton's with newly found starting point
        [all_iters, astar, ier] = newtons(f, fprime, a, tol, Nmax)
        return [astar, ier, count + len(all_iters)]

    astar = a
    ier = 2
    return [astar,ier, count]

