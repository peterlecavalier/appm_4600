import sys
import os
sys.path.insert(1, os.path.abspath('.\..\..'))
from mypkg.Quadrature import MyQuad

import numpy as np

##### Problem 1 #####
# Preset numpy to not raise a warning with 0 division
np.seterr(divide='ignore')
# Define function
# nan_to_num is used here to set evaluation at 0 to 0 (instead of nan)
f = lambda t: t*np.cos(np.nan_to_num(1/t, nan=0))

# Initialize our integration techniques
q = MyQuad(f, 0, 1)

# Run Composite Simpson's with 4 intervals (5 nodes)
comp_simp_result = q.composite_simp(n=4)

# Print the result of our approximation
print(f"Composite Simpson's with 5 nodes gives approximation = {comp_simp_result}")
##### END #####