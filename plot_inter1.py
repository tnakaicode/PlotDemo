import matplotlib.pyplot as plt
from ipywidgets import interact
import numpy as np

x = np.linspace(0, 2 * np.pi, num=5000)
y = []
param = 50

for i in range(1, param + 1):
    y.append([np.sin(i * j) for j in x])


def f(k):
    plt.plot(x, y[k])
    plt.show()


interact(f, k=(0, param - 1))
