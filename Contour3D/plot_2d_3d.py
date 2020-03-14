import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm

sys.path.append(os.path.join('../'))
from base import plot2d, plot3d

px = np.linspace(-1, 1, 100) * 100
py = np.linspace(-1, 1, 200) * 200
mesh = np.meshgrid(px, py)
func = 2 * mesh[0] + 3 * mesh[1]**2

obj = plot3d()
obj.axs.contourf(*mesh, func, cmap="jet")
obj.SavePng_Serial()

obj.new_2Dfig()
obj.axs.contourf(*mesh, func, cmap="jet")
obj.SavePng_Serial()
