# adding the package path to import mypkg classes
import sys
import os
sys.path.insert(1, os.path.abspath('.\..\..'))
from mypkg.my2DPlot import my2DPlotVect as myplt

import numpy as np

##### Problem 4b #####
theta = np.linspace(0, 2*np.pi, 1000)
# First part:
R = 1.2
delta_r = 0.1
f = 15
p = 0

x = R * (1 + delta_r * np.sin(f*theta + p)) * np.cos(theta)
y = R * (1 + delta_r * np.sin(f*theta + p)) * np.sin(theta)

plt = myplt(x, y)
plt.labels('x','y')
plt.save('4b_plot1.png')
plt.show()

# Second part:
for i in range(1, 11):
    R = i
    delta_r = 0.05
    f = 2 + i
    p = np.random.uniform(0, 2)

    x = R * (1 + delta_r * np.sin(f*theta + p)) * np.cos(theta)
    y = R * (1 + delta_r * np.sin(f*theta + p)) * np.sin(theta)
    if i == 1:
        plt = myplt(x, y)
    else:
        plt.addPlot(x, y)
plt.labels('x','y')
plt.save('4b_plot2.png')
plt.show()
##### END #####