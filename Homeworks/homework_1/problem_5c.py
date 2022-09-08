import matplotlib.pyplot as plt
import numpy as np
##### Problem 5c #####
def plot_prob(x, filename):
    # logarithmically spaced values in the interval
    delta = np.logspace(-16, -1, 16)
    # original and modified eqn
    new_eqn = lambda delt : -2*np.sin(x + delt/2)*np.sin(delt/2)

    # I divided the addition of delta by a constant
    # I tried to mess with this and find a suitable value given the x.
    # Seems like dividing the delta by 3 before adding it produces a reasonably small error
    def approx_eqn(delt):
        xi = [x] * len(delt) + (delt / 3) # division happens here
        return -delt*np.sin(x) - ((delt**2)/2)*np.cos(xi)

    new_y = new_eqn(delta)
    approx_y = approx_eqn(delta)

    # Absolute value of error
    y_diff = np.abs(approx_y - new_y)

    # Plotting difference between old and new equation
    plt.plot(delta, y_diff)

    # Log scale for both x and y
    plt.xscale('log')
    plt.yscale('log')
    # Formatting
    plt.xlabel('delta')
    plt.ylabel('Absolute difference in y (original - new)')

    plt.title(f'Difference between modified y calculation and Taylor approx. \n x = {x:.3f}')

    plt.savefig(filename)

    plt.show()

# Using Euler's number as x
x = np.exp(1)
plot_prob(x, 'number5_c_plot_1.png')

# Using 10^6 as x
x = 10**6
plot_prob(x, 'number5_c_plot_2.png')
##### END #####