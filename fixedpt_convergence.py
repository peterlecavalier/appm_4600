import numpy as np
import matplotlib.pyplot as plt
########################################################
def fixedpt(f,x0,tol,Nmax):

    ''' x0 = initial guess''' 
    ''' Nmax = max number of iterations'''
    ''' tol = stopping tolerance'''

    count = 0
    # make an array of zeros of length Nmax
    x = np.zeros((Nmax,1))
    # save the initial guess
    x[0] = x0
    while (count < Nmax):
       count = count +1
       x1 = f(x0)
       # save the current iterate
       x[count] = x1
       if (abs(x1-x0) <tol):
          xstar = x1
          ier = 0
          # truncate the array to have only count entries
          x = x[0:count]
          return [xstar,x,ier,count]
       x0 = x1

    xstar = x1
    x = x[0:count]
    ier = 1
    return [xstar,x,ier,count]

def compute_order(x,xstar):
# p_{n+1}-p (from the second index to the end)
  diff1 = np.abs(x[1:]-xstar)
  # p_n-p (from the first index to the second to last)
  diff2 = np.abs(x[0:-1]-xstar)
  # linear fit to log of differences
  fit = np.polyfit(np.log(diff2.flatten()),np.log(diff1.flatten()),1)
  print('the order equation is')
  print('log(|p_{n+1}-p|) = log(lambda) + alpha*log(|p_n-p|) where')
  print('lambda = ' + str(np.exp(fit[1])))
  print('alpha = ' + str(fit[0]))
  return [fit,diff1,diff2]

def aitkens(seq, tol, Nmax):
   xstar = seq[-1]
   p_n = seq[:-2] # all but the last two elements
   p_n1 = seq[1:-1]
   p_n2 = seq[2:]
   p_nhat = np.zeros((Nmax, 1))
   for i in range(Nmax):
      # since we already sliced the correct indices,
      # we can just index into i for each vector.
      p_nhat[i] = p_n[i] - (p_n1[i] - p_n[i])**2 / (p_n2[i] - 2*p_n1[i] + p_n[i])
      if (abs(p_nhat[i]-xstar) < tol):
         p_nhat = p_nhat[:i+1]
         return p_nhat
   return p_nhat

########################################################


### Exercise 2.2
g = lambda x: np.power(10/(x+4),1/2)
[xstar,x,ier,count] = fixedpt(g,1.5,1e-10,100)
print(x)

print(f"Fixed point took {count} iterations")

# compute the convergence rate and constant
[fit, diff1, diff2] = compute_order(x, xstar)
# plot the data
plt.loglog(diff2,diff1,'ro',label='fixedpt data')
# plot the fit
plt.loglog(diff2,np.exp(fit[1]+fit[0]*np.log(diff2)),'b-',label='fixedpt fit')
# label the plot
plt.xlabel('$|p_{n}-p|$')
plt.ylabel('$|p_{n+1}-p|$')
plt.legend()
plt.savefig('lab4_fixed_pt_plot.png')
plt.show()

# Exercise 3.2

# sequence created in the pre-lab is x
print(x)
new_x_seq = aitkens(x, 1e-10, 100)
new_xstar = new_x_seq[-1]
new_count = len(new_x_seq)

# compute the convergence rate and constant
[fit, diff1, diff2] = compute_order(new_x_seq, xstar)
# plot the data
plt.loglog(diff2,diff1,'ro',label='aitkens data')
# plot the fit
plt.loglog(diff2,np.exp(fit[1]+fit[0]*np.log(diff2)),'b-',label='aitkens fit')
# label the plot
plt.xlabel('$|p_{n}-p|$')
plt.ylabel('$|p_{n+1}-p|$')
plt.legend()
plt.savefig('lab4_aitkens_plot.png')
plt.show()
