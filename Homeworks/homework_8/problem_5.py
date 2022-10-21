import numpy as np

A = [[1, 3], [6, -1], [4, 0], [2, 7]]
c = [[1], [2], [3], [4]]
D = np.diag([1, 2, 5, 3])

x = np.linalg.inv(((D @ A).T @ (D @ A))) @ ((D @ A).T) @ (D @ c)
print(x)