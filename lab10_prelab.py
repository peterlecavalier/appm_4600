import numpy as np
import matplotlib.pyplot as plt

def eval_legendre(n, x):
    vals = np.zeros((n+1))
    vals[0] = 1; vals[1] = x
    for i in range(2, n+1):
        vals[i] = (1/i) * ((2*(i-1)+1)*x*vals[i-1] - (i-1)*vals[i-2])
    return vals


# Testing that the Legendre works
# (should match the first picture here: https://en.wikipedia.org/wiki/Legendre_polynomials)
all_vals = np.zeros((1000, 6))

x = np.linspace(-1, 1, 1000)
for idx, i in enumerate(x):
    cur = eval_legendre(5, i)
    all_vals[idx] = cur

for i in range(6):
    plt.plot(x, all_vals[:, i], label=f'P{i}')

plt.legend()
plt.show()