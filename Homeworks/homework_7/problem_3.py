import numpy as np
import matplotlib.pyplot as plt

##### Problem 3 #####
f = lambda x: 1/(1+(10*x)**2)

def eval_barycentric(xeval, xs, fxs):
    phi = np.prod(np.array([(xeval - i) for i in xs]))
    total_sum = 0
    for idx, xj in enumerate(xs):
        wj = np.prod(np.array([(1/(xj - xi)) for xi in xs if xi != xj]))
        fxj = fxs[idx]
        total_sum += (wj / (xeval - xj))*fxj
    return phi*total_sum

# Finer grain points for plotting f(x)
x_fine = np.linspace(-1, 1, 1001)
fx_fine = f(x_fine)

fig, ax = plt.subplots(6, 3, figsize=(16, 10), sharex=True)

for n in range(2, 20):
    # All xi points
    xi = np.linspace(1, n, n)
    xi = np.cos((2*(xi) - 1)*np.pi / (2*(len(xi))))
    fxi = f(xi)

    xfine_px = []
    px = []
    for idx, i in enumerate(x_fine):
        if i not in xi:
            xfine_px.append(i)
            px.append(eval_barycentric(i, xi, fxi))
    px = np.array(px)

    # Calculate where to put the plot
    row = int(np.floor((n-2)/3))
    col = (n-2) % 3

    # Plot everything
    ax[row][col].plot(xi, fxi, 'o')
    ax[row][col].plot(xfine_px, px, label=f'p(x), N={n}')
    ax[row][col].plot(x_fine, fx_fine, label='f(x)')
    ax[row][col].legend()
plt.suptitle('Barycentric (Chebyshev) Interpolation')
plt.savefig('problem3_plot.png')
plt.show()
##### END #####