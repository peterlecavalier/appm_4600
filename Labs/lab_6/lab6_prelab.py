from mypkg.Iteration1D import Iteration1D
import numpy as np

##### Part 1 #####
h = 0.01 * 2. ** (-np.arange(0, 10))

f = lambda x : np.cos(x)
x = np.pi / 2


f_diff = (f(x + h) - f(x)) / h
c_diff = (f(x + h) - f(x - h)) / (2*h)

print(f"Forward Difference: {f_diff}")
print(f"Centered Difference: {c_diff}")
##### END #####

##### Part 2 #####
it1d = Iteration1D(None, None)
it1d.compute_order(f_diff[:6], -1)
it1d.compute_order(c_diff[:6], -1)
##### END #####