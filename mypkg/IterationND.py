from numpy.linalg import inv 
from numpy.linalg import norm 

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
