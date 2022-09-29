from numpy.linalg import inv 
from numpy.linalg import norm
import numpy as np 

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
    
    def Newton(self, x0=None,tol=None,Nmax=None):

        ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
        ''' Outputs: xstar= approx root, ier = error message, its = num its'''

        # Verify that needed vars are inputted somewhere
        if x0 is None:
            if self.x0 is None:
                print("ERROR: Please input x0")
                return
            else:
                x0 = self.x0
        if tol is None:
            if self.tol is None:
                print("ERROR: Please input tol")
                return
            else:
                tol = self.tol
        if Nmax is None:
            if self.Nmax is None:
                print("ERROR: Please input Nmax")
                return
            else:
                Nmax = self.Nmax
        

        for its in range(Nmax):
            J = self.__evalJ(x0)
            Jinv = inv(J)
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
            
    def LazyNewton(self,x0=None,tol=None,Nmax=None):

        ''' Lazy Newton = use only the inverse of the Jacobian for initial guess'''
        ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
        ''' Outputs: xstar= approx root, ier = error message, its = num its'''

        # Verify that needed vars are inputted somewhere
        if x0 is None:
            if self.x0 is None:
                print("ERROR: Please input x0")
                return
            else:
                x0 = self.x0
        if tol is None:
            if self.tol is None:
                print("ERROR: Please input tol")
                return
            else:
                tol = self.tol
        if Nmax is None:
            if self.Nmax is None:
                print("ERROR: Please input Nmax")
                return
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
                print("ERROR: Please input x0")
                return
            else:
                x0 = self.x0
        if tol is None:
            if self.tol is None:
                print("ERROR: Please input tol")
                return
            else:
                tol = self.tol
        if Nmax is None:
            if self.Nmax is None:
                print("ERROR: Please input Nmax")
                return
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
