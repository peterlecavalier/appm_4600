import matplotlib.pyplot as plt
import numpy as np

def find_points(xeval, xint, idx):
  """
  xeval is the larger array of points
  xint is the smaller array of intervals
  idx is the interval index (counting from 1)
  """
  indices = np.where(np.logical_and(xeval <= xint[idx+1], xeval >= xint[idx]))[0]
  return indices

def eval_line(x0, fx0, x1, fx1, xeval):
  m = (fx1 - fx0)/(x1 - x0)
  return m*xeval - m*x0 + fx0


def driver():
    
  f = lambda x: 1/(1+(10*x)**2)
  a = -1
  b = 1
  
  ''' create points you want to evaluate at'''
  Neval = 100
  xeval =  np.linspace(a,b,Neval)
  
  ''' number of intervals'''
  Nint = 10
  
  '''evaluate the linear spline'''
  yeval = eval_lin_spline(xeval,Neval,a,b,f,Nint)
  
  ''' evaluate f at the evaluation points'''
  fex = np.zeros(Neval)
  for j in range(Neval):
    fex[j] = f(xeval[j]) 

  fig, ax = plt.subplots(1, 2)
  ax[0].plot(xeval, yeval)
  ax[0].set_title('Linear Spline')

  err = abs(yeval-fex)
  ax[1].plot(xeval, err)
  ax[1].set_title('Error in spline calc')
  plt.savefig('lab8_plot.png')
  plt.show()

    
    
def  eval_lin_spline(xeval,Neval,a,b,f,Nint):
  '''create the intervals for piecewise approximations'''
  xint = np.linspace(a,b,Nint+1)
  xeval = np.linspace(a,b, Neval)
  
  '''create vector to store the evaluation of the linear splines'''
  yeval = np.zeros(Neval)
  
  for jint in range(Nint):
    '''find indices of xeval in interval (xint(jint),xint(jint+1))'''
    '''let ind denote the indices in the intervals'''
    '''let n denote the length of ind'''
    ind = find_points(xeval, xint, jint)
    '''temporarily store your info for creating a line in the interval of 
      interest'''
    a1= xint[jint]
    fa1 = f(a1)
    b1 = xint[jint+1]
    fb1 = f(b1)
    
    for kk in ind:
        '''use your line evaluator to evaluate the lines at each of the points 
        in the interval'''
        '''yeval(ind(kk)) = call your line evaluator at xeval(ind(kk)) with 
        the points (a1,fa1) and (b1,fb1)'''
        yeval[kk] = eval_line(a1, fa1, b1, fb1, xeval[kk])
  
  return yeval

           
if __name__ == '__main__':
      # run the drivers only if this is called from the command line
      driver()               
