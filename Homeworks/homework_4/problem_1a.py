# adding the package path to import mypkg classes
import sys
import os
sys.path.insert(1, os.path.abspath('.\..\..'))
from mypkg.my2DPlot import my2DPlot as myplt

from scipy.special import erf
import numpy as np

##### Problem 1a #####
# Plotting the original function:
plt = myplt(lambda x: 35*erf(x/(2*np.sqrt(0.715392))) - 15,0,1, 'Original f(x)')
plt.labels('x - Depth (meters)','y - Temperature (Celsius)')
plt.save('1a_plot.png')
plt.show()
##### END #####