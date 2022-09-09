# adding the package path to import mypkg classes
import sys
import os
sys.path.insert(1, os.path.abspath('.\..\..'))
import mypkg.prini as prini

import numpy as np

##### Problem 3d #####
t = np.linspace(0, np.pi, 31)
y = np.cos(t)

print(t)
print(y)

S = 0
for i in range(len(t)):
    S += (t[i] * y[i])

p = prini("real", "the sum is:", S)

p.print()
##### END #####