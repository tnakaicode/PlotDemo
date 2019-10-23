""" From "A SURVEY OF COMPUTATIONAL PHYSICS", Python eBook Version
   by RH Landau, MJ Paez, and CC Bordeianu
   Copyright Princeton University Press, Princeton, 2012; Book  Copyright R Landau, 
   Oregon State Unv, MJ Paez, Univ Antioquia, C Bordeianu, Univ Bucharest, 2012.
   Support by National Science Foundation , Oregon State Univ, Microsoft Corp"""

# Simple3Dplot.py: Simple matplotlib 3D plot; rotate & scale via mouse

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mayavi.mlab import *

from base import plot3d

obj = plot3d()

print("Please be patient, I have packages to import & points to plot")
delta = 0.1
x = np.arange(-3., 3., delta)
y = np.arange(-3., 3., delta)
X, Y = np.meshgrid(x, y)
Z = np.sin(X) * np.cos(Y)

mesh(X, Y, Z, colormap='YlGnBu')
# surface height
#s = surf(X, Y, Z)

plt.show()
