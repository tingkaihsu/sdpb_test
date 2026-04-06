import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0, 2, 200)

y1 = -2.483

y2 = -1.840

f1 = 1 + x**4 + y1 * ( x**4/12 + x**2)

f2 = 1 + x**4 + y2 * ( x**4/12 + x**2)

plt.plot(x, f1, label='y = -2.483')
plt.plot(x, f2, label='y = -1.840')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title('Plot of f(x) for different values of y')
plt.legend()
plt.grid()
plt.show()
