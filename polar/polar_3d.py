import numpy as np
import sys
import time
import os

sys.path.append(os.path.join("../"))
from base import plot3d

# Make data
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = 10 * np.outer(np.cos(u), np.sin(v))
y = 10 * np.outer(np.sin(u), np.sin(v))
z = 10 * np.outer(np.ones(np.size(u)), np.cos(v))

obj = plot3d()
obj.axs.plot_surface(x, y, z, color='b')
#obj.axs.plot_wireframe(x, y, z, color='k', rstride=1, cstride=1)
obj.SavePng()
