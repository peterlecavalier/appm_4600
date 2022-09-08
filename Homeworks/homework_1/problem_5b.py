import matplotlib.pyplot as plt
import numpy as np
##### Problem 5b #####
def plot_prob(x, filename):
    # logarithmically spaced values in the interval
    delta = np.logspace(-16, 0, 17)
    # original and modified eqn
    og_eqn = lambda delt : np.cos(x + delt) - np.cos(x)
    new_eqn = lambda delt : -2*np.sin(x + delt/2)*np.sin(delt/2)

    og_y = og_eqn(delta)
    new_y = new_eqn(delta)

    y_diff = og_y - new_y

    # Plotting difference between old and new equation
    plt.plot(delta, y_diff)

    plt.xscale('log')
    # Formatting
    plt.xlabel('delta')
    plt.ylabel('Difference in y (original - new)')

    plt.title(f'Difference in y, where x = {x:.3f}')

    plt.savefig(filename)

    plt.show()

# Using Euler's number as x
x = np.exp(1)
plot_prob(x, 'number5_b_plot_1.png')

# Using 10^6 as x

x = 10**6
plot_prob(x, 'number5_b_plot_2.png')
##### END #####