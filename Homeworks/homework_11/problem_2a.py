import sys
import os
sys.path.insert(1, os.path.abspath('.\..\..'))
from mypkg.Quadrature import MyQuad

import numpy as np
from scipy.special import gamma
import matplotlib.pyplot as plt

##### Problem 2a #####
''' Experimentation to find endpoint value'''
# First, we can just set our number of intervals to a very high value
# Let's do n=100000:
q = MyQuad(a=0, n=100000)
# Here's the actual gamma values:
actual_gammas = [gamma(x) for x in [2, 4, 6, 8, 10]]

# Set up our right endpoints to loop through and test
endpoints = np.arange(1, 150, 1)
# Create an errors array to store the relative error of each endpoint
errors = np.zeros(endpoints.shape, dtype=np.float64)

for idx, i in enumerate(endpoints):
    results = []
    for gamma_idx, x in enumerate([2, 4, 6, 8, 10]):
        # Define function, iterating through all x values for experimenting
        f = lambda t: t**(x-1) * np.exp(-t)
        cur_result = q.composite_trap(f=f, b=i)
        rel_error = np.abs((cur_result - actual_gammas[gamma_idx]) / actual_gammas[gamma_idx])
        results.append(rel_error)
    # rel_error is now the mean of the error for each x value
    rel_error = np.mean(results)

    errors[idx] = rel_error

# Find the minimum error endpoint value
min_x = endpoints[np.argmin(errors)]
min_y = errors[np.argmin(errors)]

plt.figure(figsize=(14, 7))
plt.semilogy(endpoints, errors)
plt.plot(min_x, min_y, 'ro')
plt.annotate(f"     Minimum error {min_y:.2} at b={min_x}", (min_x, min_y))
plt.xlabel('Right endpoint value')
plt.ylabel('Average relative error')
plt.title('Average error for the gamma function approximation over x=[2, 4, 6, 8, 10] using Composite Trapezoidal rule')
plt.savefig('2a_endpoint_plot.png')
plt.show()

''' Experimentation to find number of intervals'''
# Now that we've picked an endpoint of b=42, we can further expirement to find a good cutoff for our number of intervals.

# In a similar fashion as before, approach is to loop through interval values

# Now, let's put in our found b value of 42
q = MyQuad(a=0, b=42)
# Here's the actual gamma values:
actual_gammas = [gamma(x) for x in [2, 4, 6, 8, 10]]

# Set up our interval numbers to loop through and test
intervals = np.arange(1, 10000, 10)
# Create an errors array to store the relative error of each endpoint
errors = np.zeros(intervals.shape, dtype=np.float64)

for idx, i in enumerate(intervals):
    results = []
    for gamma_idx, x in enumerate([2, 4, 6, 8, 10]):
        # Define function, iterating through all x values for experimenting
        f = lambda t: t**(x-1) * np.exp(-t)
        cur_result = q.composite_trap(f=f, n=i)
        rel_error = np.abs((cur_result - actual_gammas[gamma_idx]) / actual_gammas[gamma_idx])
        results.append(rel_error)
    # rel_error is now the mean of the error for each x value
    rel_error = np.mean(results)

    errors[idx] = rel_error

# After initially looking at the plot below, I realized that I wanted to find the point
# that gave the first crossover into 10^-6 order error:
good_x = intervals[np.where(errors < 10**-5)[0][0]]
good_y = errors[np.where(errors < 10**-5)[0][0]]

plt.figure(figsize=(14, 7))
plt.semilogy(intervals, errors)
plt.plot(good_x, good_y, 'ro')
plt.annotate(f"     Chosen point has error {good_y:.2} at n={good_x} intervals", (good_x, good_y))
plt.xlabel('Number of intervals')
plt.ylabel('Average relative error')
plt.title('Average error for the gamma function approximation over x=[2, 4, 6, 8, 10] using Composite Trapezoidal rule')
plt.savefig('2a_intervals_plot.png')
plt.show()

# Thus, we have shown that choosing b=42 as our endpoint with n=1721 intervals will give us
# a suitable balance between efficiency and accuracy.

'''
END EXPERIMENTATION
'''

# Now, we can apply our previous findings to the approximation.

# Set up quadrature with our previous findings
q = MyQuad(a=0, b=42, n=1721)
xs = [2, 4, 6, 8, 10]

print("The following approximates the gamma function on the interval [0,42]")
print("using Composite Trapezoidal rule with n=1721 intervals:")

for x in xs:
    f = lambda t: t**(x-1) * np.exp(-t)
    cur_result = q.composite_trap(f=f)
    actual_gamma = gamma(x)
    rel_error = np.abs((cur_result - actual_gamma) / actual_gamma)
    print(f"The relative error for x={x} is {rel_error:.4}")
##### END #####