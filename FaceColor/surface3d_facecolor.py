import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
import time
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm

sys.path.append(os.path.join('../'))
from base import plot2d, plot3d

X = np.loadtxt("./surface3d_facecolor_X.xxt")
Y = np.loadtxt("./surface3d_facecolor_Y.xxt")
Z = np.loadtxt("./surface3d_facecolor_Z.xxt")
V = np.loadtxt("./surface3d_facecolor_V.xxt")

obj = plot3d()
obj.axs.view_init(45, 60)

# Normalize in [0, 1] the DataFrame V that defines the color of the surface.
V_normalized = (V - V.min().min())
V_normalized = V_normalized / V_normalized.max().max()

obj.axs.plot_surface(X, Y, Z, facecolors=plt.cm.jet(V_normalized))
obj.axs.set_xlabel('x', fontsize=18)
obj.axs.set_ylabel('y', fontsize=18)
obj.axs.set_zlabel('z', fontsize=18)

m = cm.ScalarMappable(cmap=cm.jet)
m.set_array(V)
plt.colorbar(m)

obj.SavePng_Serial()
