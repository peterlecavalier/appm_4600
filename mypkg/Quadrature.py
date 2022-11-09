import numpy as np

class MyQuad():
    def __init__(self,f=None,a=None,b=None, n=None):
        '''
        f: function
        a: lower endpoint
        b: higher endpoint
        n: number of intervals
        '''
        self.f = f
        self.a = a
        self.b = b
        self.n = n

    def composite_trap(self,f=None,a=None,b=None,n=None):
        # Verify that needed vars are inputted somewhere
        if f is None:
            if self.f is None:
                raise ValueError("Please input f")
            else:
                f = self.f
        if a is None:
            if self.a is None:
                raise ValueError("Please input a")
            else:
                a = self.a
        if b is None:
            if self.b is None:
                raise ValueError("Please input b")
            else:
                b = self.b
        if n is None:
            if self.n is None:
                raise ValueError("Please input n")
            else:
                n = self.n
        
        # Validate endpoints
        if b < a:
            raise ValueError("b must be greater than a")
        elif b == a:
            return 0
        
        h = (b - a)/n

        xj = np.linspace(a, b, n+1)
        fxj = f(xj)

        return (h/2)*(fxj[0] + 2*np.sum(fxj[1:n]) + fxj[-1])

    def composite_simp(self,f=None,a=None,b=None,n=None):
        # Verify that needed vars are inputted somewhere
        if f is None:
            if self.f is None:
                raise ValueError("Please input f")
            else:
                f = self.f
        if a is None:
            if self.a is None:
                raise ValueError("Please input a")
            else:
                a = self.a
        if b is None:
            if self.b is None:
                raise ValueError("Please input b")
            else:
                b = self.b
        if n is None:
            if self.n is None:
                raise ValueError("Please input n")
            else:
                n = self.n
        
        # Validate endpoints
        if b < a:
            raise ValueError("b must be greater than a")
        elif b == a:
            return 0
        
        # Validate that n is even
        if n % 2 != 0:
            raise ValueError("n must be even")
                
        h = (b - a)/n

        xj = np.linspace(a, b, n+1)
        fxj = f(xj)

        return (h/3)*(fxj[0] + 2*np.sum(fxj[2:n:2]) + 4* np.sum(fxj[1:n:2]) + fxj[-1])
        
