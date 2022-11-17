import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la
import scipy.linalg as scila
from time import time

def driver():

     ''' create  matrix for testing different ways of solving a square 
     linear system'''

     '''' N = size of system'''
     for N in [100, 500, 1000, 2000, 4000, 5000]:
          '''Right hand side'''
          b = np.random.rand(N,1)
          A = np.random.rand(N,N)

          normal_start = time()
          x = scila.solve(A,b)
          

          # Instead, solve with LU-fact
          fact_start = time()
          LU, P = scila.lu_factor(A)
          fact_end = time()
          xlu = scila.lu_solve((LU, P), b)
          lu_end = time()

          test = np.matmul(A,x)
          r = la.norm(test-b)

          # Now, we can additionally calculate the lu residuals
          test_lu = np.matmul(A, xlu)
          r_lu = la.norm(test_lu - b)

          '''
          print(f"Residual without LU: {r}")
          print(f"Residual with LU: {r_lu}")
          '''
          
          print(f'---------N = {N}---------')
          print(f"Regular solver: {fact_start - normal_start} seconds")
          print(f"LU Factorization: {fact_end - fact_start} seconds")
          print(f"LU Solve: {lu_end - fact_end} seconds")
          print(f"Total LU: {lu_end - fact_start} seconds")



     ''' Create an ill-conditioned rectangular matrix '''
     N = 10
     M = 5
     A = create_rect(N,M)     
     b = np.random.rand(N,1)


     
def create_rect(N,M):
     ''' this subroutine creates an ill-conditioned rectangular matrix'''
     a = np.linspace(1,10,M)
     d = 10**(-a)
     
     D2 = np.zeros((N,M))
     for j in range(0,M):
        D2[j,j] = d[j]
     
     '''' create matrices needed to manufacture the low rank matrix'''
     A = np.random.rand(N,N)
     Q1, R = la.qr(A)
     test = np.matmul(Q1,R)
     A =    np.random.rand(M,M)
     Q2,R = la.qr(A)
     test = np.matmul(Q2,R)
     
     B = np.matmul(Q1,D2)
     B = np.matmul(B,Q2)
     return B     
          
  
if __name__ == '__main__':
      # run the drivers only if this is called from the command line
      driver()       
