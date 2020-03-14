import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
import matplotlib.colors

sys.path.append(os.path.join('../'))
from base import plot2d, plot3d

x = np.arange(3)
X, Y = np.meshgrid(x, x)
Z = np.ones_like(X)

V = np.array([[3, 2, 2], [1, 0, 3], [2, 1, 0]])

obj = plot3d()
norm = matplotlib.colors.Normalize(vmin=0, vmax=3)
obj.axs.plot_surface(X, Y, Z, facecolors=plt.cm.jet(norm(V)), shade=False)
m = cm.ScalarMappable(cmap=plt.cm.jet, norm=norm)
#m.set_array([])
obj.fig.colorbar(m)
obj.SavePng_Serial()
