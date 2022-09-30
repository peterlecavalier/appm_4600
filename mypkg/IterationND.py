from numpy.linalg import inv 
from numpy.linalg import norm
import numpy as np 
import matplotlib.pyplot as plt

'''
MODIFIED FROM CODE PROVIDED BY ADRIANNA GILLMAN
'''

class nd_iteration():
    def __init__(self, F, J):
        self.__evalF = F
        self.__evalJ = J
        
        self.x0 = None
        self.tol = None
        self.Nmax = None
    
    def compute_order(self, x, xstar, fig_fp=None):
      diff1 = x[1:] - xstar
      diff1 = np.array([np.linalg.norm(i) for i in diff1])
      # p_n-p (from the first index to the second to last)
      diff2 = x[0:-1] - xstar
      diff2 = np.array([np.linalg.norm(i) for i in diff2])
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
      plt.loglog(diff2,diff1,'ro',label='iteration data')
      # plot the fit
      plt.loglog(diff2,np.exp(fit[1]+fit[0]*np.log(diff2)),'b-',label='fit')
      # label the plot
      plt.xlabel('$|p_{n}-p|$')
      plt.ylabel('$|p_{n+1}-p|$')
      plt.legend()
      plt.title(f'lambda = {np.exp(fit[1])}, alpha {fit[0]}')
      if fig_fp is not None:
        plt.savefig(fig_fp)
      plt.show()
      return [fit,diff1,diff2]

    def Newton(self, x0=None,tol=None,Nmax=None):

        ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
        ''' Outputs: xstar= approx root, ier = error message, its = num its'''

        # Verify that needed vars are inputted somewhere
        if x0 is None:
            if self.x0 is None:
                raise ValueError("Please input x0")
            else:
                x0 = self.x0
        if tol is None:
            if self.tol is None:
                raise ValueError("Please input tol")
            else:
                tol = self.tol
        if Nmax is None:
            if self.Nmax is None:
                raise ValueError("Please input Nmax")
            else:
                Nmax = self.Nmax

        all_x = np.zeros((Nmax + 1, x0.shape[0]))

        for its in range(Nmax):
            all_x[its] = x0
            J = self.__evalJ(x0)
            Jinv = inv(J)
            F = self.__evalF(x0)
            
            x1 = np.squeeze(np.expand_dims(x0, -1) - np.matmul(Jinv, F))
            
            if (norm(x1-x0) < tol):
                all_x[its+1] = x1
                all_x = all_x[:its+2]
                xstar = x1
                ier = 0
                return [all_x, xstar, ier, its]
                
            x0 = x1
        
        all_x[its+1] = x1
        all_x = all_x[:its+2]
        xstar = x1
        ier = 1
        return [all_x, xstar,ier,its]
            
    def LazyNewton(self,x0=None,tol=None,Nmax=None):

        ''' Lazy Newton = use only the inverse of the Jacobian for initial guess'''
        ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
        ''' Outputs: xstar= approx root, ier = error message, its = num its'''

        # Verify that needed vars are inputted somewhere
        if x0 is None:
            if self.x0 is None:
                raise ValueError("Please input x0")
            else:
                x0 = self.x0
        if tol is None:
            if self.tol is None:
                raise ValueError("Please input tol")
            else:
                tol = self.tol
        if Nmax is None:
            if self.Nmax is None:
                raise ValueError("Please input Nmax")
            else:
                Nmax = self.Nmax

        J = self.__evalJ(x0)
        Jinv = inv(J)
        for its in range(Nmax):

            F = self.__evalF(x0)
            x1 = x0 - Jinv.dot(F)
            
            if (norm(x1-x0) < tol):
                xstar = x1
                ier = 0
                return[xstar, ier, its]
                
            x0 = x1
        
        xstar = x1
        ier = 1
        return[xstar,ier,its]

    def Broyden(self,x0=None,tol=None,Nmax=None):
        '''tol = desired accuracy
        Nmax = max number of iterations'''

        '''Sherman-Morrison 
    (A+xy^T)^{-1} = A^{-1}-1/p*(A^{-1}xy^TA^{-1})
        where p = 1+y^TA^{-1}Ax'''

        '''In Newton
        x_k+1 = xk -(G(x_k))^{-1}*F(x_k)'''


        '''In Broyden 
        x = [F(xk)-F(xk-1)-\hat{G}_k-1(xk-xk-1)
        y = x_k-x_k-1/||x_k-x_k-1||^2'''

        ''' implemented as in equation (10.16) on page 650 of text'''
        
        '''initialize with 1 newton step'''
        
        A0 = self.__evalJ(x0)

        v = self.__evalF(x0)
        A = inv(A0)

        s = -A.dot(v)
        xk = x0+s
        for  its in range(Nmax):
            '''(save v from previous step)'''
            w = v
            ''' create new v'''
            v = self.__evalF(xk)
            '''y_k = F(xk)-F(xk-1)'''
            y = v-w;                   
            '''-A_{k-1}^{-1}y_k'''
            z = -A.dot(y)
            ''' p = s_k^tA_{k-1}^{-1}y_k'''
            p = -np.dot(s,z)                 
            u = np.dot(s,A) 
            ''' A = A_k^{-1} via Morrison formula'''
            tmp = s+z
            tmp2 = np.outer(tmp,u)
            A = A+1./p*tmp2
            ''' -A_k^{-1}F(x_k)'''
            s = -A.dot(v)
            xk = xk+s
            if (norm(s)<tol):
                alpha = xk
                ier = 0
                return[alpha,ier,its]
        alpha = xk
        ier = 1
        return[alpha,ier,its]
    

    def SlackerNewton(self,x0=None,tol=None,Nmax=None):

        ''' Lazy Newton = use only the inverse of the Jacobian for initial guess'''
        ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
        ''' Outputs: xstar= approx root, ier = error message, its = num its'''

        # Verify that needed vars are inputted somewhere
        if x0 is None:
            if self.x0 is None:
                raise ValueError("Please input x0")
            else:
                x0 = self.x0
        if tol is None:
            if self.tol is None:
                raise ValueError("Please input tol")
            else:
                tol = self.tol
        if Nmax is None:
            if self.Nmax is None:
                raise ValueError("Please input Nmax")
            else:
                Nmax = self.Nmax

        J = self.__evalJ(x0)
        Jinv = inv(J)
        update_tol = tol * 1e4
        for its in range(Nmax):

            F = self.__evalF(x0)
            # norm of dot product is less than 10e-6, recompute
            curr_jdot = Jinv.dot(F)
            x1 = x0 - curr_jdot
            
            if (norm(x1-x0) < tol):
                xstar = x1
                ier = 0
                return[xstar, ier, its]
            # Here is where we recompute the Jacobian
            # Recompute if the error is greater than the last one
            elif (norm(curr_jdot) < update_tol):
                J = self.__evalJ(x0)
                Jinv = inv(J)

            #last_err = norm(x1 - x0) 
            x0 = x1
        
        xstar = x1
        ier = 1
        return[xstar,ier,its]
