# adding the package path to import mypkg classes
import sys
import os
sys.path.insert(1, os.path.abspath('.\..\..'))
from mypkg.my2DPlot import my2DPlot as myplt

import numpy as np

##### Problem 5a #####
plt = myplt(lambda x: x - 4*np.sin(2*x) - 3,-2,6)
plt.labels('x','y')
# Adding y = 0 to the plot
plt.addPlot(lambda x: 0*x)
plt.dotted()
plt.color('black')
plt.save('5a_plot.png')
plt.show()
##### END #####