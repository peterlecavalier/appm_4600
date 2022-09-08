import numpy as np
##### Problem 2b #####
A = (1/2.) * np.array([[1, 1],[1 + 10**-10, 1 - 10**-10]])
_, svds, _ = np.linalg.svd(A)
#print(svds)
cond_numb = svds[0] / svds[1]
print(np.format_float_scientific(cond_numb))
##### END #####