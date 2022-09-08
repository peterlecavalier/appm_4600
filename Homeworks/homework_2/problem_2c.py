import numpy as np
##### Problem 2c #####
a_inv_b_norm = np.sqrt((3*10**5 + 10**-5)**2 + (-3*10**5 + 10**-5)**2)
print(a_inv_b_norm)
# divide by sqrt(2) since that's the norm of x
print(a_inv_b_norm / np.sqrt(2))
##### END #####
