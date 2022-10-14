import numpy as np
import matplotlib.pyplot as plt

##### Problem 2 #####
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

fig, ax = plt.subplots(2, 3, figsize=(16, 10), sharex=True)
fig2, ax2 = plt.subplots(2, 3, figsize=(16, 10), sharex=True)

# Making an n-values list to get specific ns
n_vals = [2, 10, 19, 50, 100, 300]
for n_i, n in enumerate(n_vals):
    # All xi points
    xi = np.linspace(1, n, n)
    h = 2/(len(xi) - 1)
    xi = -1 + (xi - 1)*h
    fxi = f(xi)

    xfine_px = []
    px = []
    new_idxs = []
    for idx, i in enumerate(x_fine):
        if i not in xi:
            new_idxs.append(idx)
            xfine_px.append(i)
            px.append(eval_barycentric(i, xi, fxi))
    px = np.array(px)

    # Calculate where to put the plot
    row = int(np.floor(n_i/3))
    col = int(n_i % 3)

    # Plot everything
    ax[row][col].plot(xi, fxi, 'o')
    ax[row][col].plot(xfine_px, px, label=f'p(x), N={n}')
    ax[row][col].plot(x_fine, fx_fine, label='f(x)')
    ax[row][col].legend()

    new_fx_fine = fx_fine[new_idxs]
    err = np.abs(new_fx_fine - px)

    ax2[row][col].plot(xfine_px, err, label=f'p(x) error, N={n}')
    ax2[row][col].legend()
    ax2[row][col].set_yscale('log')
    if n == 1001:
        print(err)

fig.suptitle('Barycentric Interpolation')
fig.savefig('problem2_plot1.png')

fig2.suptitle('Barycentric Interpolation Error')
fig2.savefig('problem2_plot2.png')

plt.show()
##### END #####