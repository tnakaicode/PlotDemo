'''
===============
Filled contours
===============

contourf differs from contour in that it creates filled contours, ie.
a discrete number of colours are used to shade the domain.

This is like a contourf plot in 2D except that the shaded region corresponding
to the level c is graphed on the plane z=c.
'''

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm

sys.path.append(os.path.join('../'))
from base import plot3d

X, Y, Z = axes3d.get_test_data(0.01)

obj = plot3d()
img = obj.axs.contourf(X, Y, Z, cmap=cm.coolwarm)
obj.axs.clabel(img, fontsize=9, inline=1)
obj.SavePng_Serial()

obj.new_fig()
img = obj.axs.contour(X, Y, Z, cmap=cm.coolwarm)
obj.axs.clabel(img, fontsize=9, inline=1)
obj.SavePng_Serial()
