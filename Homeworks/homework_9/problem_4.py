import numpy as np
import matplotlib.pyplot as plt

maclaurin = lambda x: x - (x**3)/6 + (x**5)/120
fx = lambda x: np.sin(x)

# P_3^3 is equivalent to P_2^4.
p_3_3 = lambda x: (x - (7/60)*x**3) / (1 + (1/20)*x**2)
p_4_2 = lambda x: x / (1 + (1/6)*x**2 + (7/360)*x**4)

x = np.linspace(0, 5, 1000)
y_3_3 = p_3_3(x)
y_4_2 = p_4_2(x)

y_mac = maclaurin(x)
y_sin = fx(x)


fig, ax = plt.subplots(3, 2, figsize=(12,9))
ax[0][0].plot(x, y_3_3, label='P(3,3)')
ax[0][0].plot(x, y_mac, label='6th order Maclaurin')
ax[0][0].plot(x, y_sin, label='sin(x)')
ax[0][1].plot(x, np.abs(y_mac - y_sin), label='Maclaurin Error')
ax[0][1].plot(x, np.abs(y_3_3 - y_sin), label='Pade Error')
ax[0][0].legend()
ax[0][1].legend()
ax[0][0].set_title('Maclaurin Series vs. Pade plots')
ax[0][1].set_title('Absolute Error plots')

ax[1][0].plot(x, y_4_2, label='P(4,2)')
ax[1][0].plot(x, y_mac, label='6th order Maclaurin')
ax[1][0].plot(x, y_sin, label='sin(x)')
ax[1][1].plot(x, np.abs(y_mac - y_sin), label='Maclaurin Error')
ax[1][1].plot(x, np.abs(y_4_2 - y_sin), label='Pade Error')
ax[1][0].legend()
ax[1][1].legend()

ax[2][0].plot(x, y_3_3, label='P(2,4)')
ax[2][0].plot(x, y_mac, label='6th order Maclaurin')
ax[2][0].plot(x, y_sin, label='sin(x)')
ax[2][1].plot(x, np.abs(y_mac - y_sin), label='Maclaurin Error')
ax[2][1].plot(x, np.abs(y_3_3 - y_sin), label='Function Error')

ax[2][0].legend()
ax[2][1].legend()

plt.savefig('problem_4_plot.png')
plt.show()