import matplotlib.pyplot as plt
import numpy as np
##### Problem 1 #####
x = np.arange(1.920,2.081,0.001)
# plot i:
p = lambda x : x**9 - 18*x**8 + 144*x**7 - 672*x**6 + 2016*x**5 - 4032*x**4 + 5376*x**3 - 4608*x**2 + 2304*x - 512
y = p(x)
plt.plot(x, y, label='Coefficients')

# plot ii:
p = lambda x : (x - 2)**9
y = p(x)
plt.plot(x, y, label='(x-2)^9')

# Formatting
plt.xlabel(x)
plt.ylabel(y)

plt.legend()

plt.savefig('number1_plot.png')

plt.show()
##### END #####