"""
=================================
Triangular 3D filled contour plot
=================================

Filled contour plots of unstructured triangular grids.

The data used is the same as in the second plot of trisurf3d_demo2.
tricontour3d_demo shows the unfilled version of this example.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import sys
import os
import time
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm

sys.path.append(os.path.join('../'))
from base import plot2d, plot3d

# First create the x, y, z coordinates of the points.
n_angles = 48
n_radii = 8
min_radius = 0.25

# Create the mesh in polar coordinates and compute x, y, z.
radii = np.linspace(min_radius, 0.95, n_radii)
angles = np.linspace(0, 2 * np.pi, n_angles, endpoint=False)
angles = np.repeat(angles[..., np.newaxis], n_radii, axis=1)
angles[:, 1::2] += np.pi / n_angles

x = (radii * np.cos(angles)).flatten()
y = (radii * np.sin(angles)).flatten()
z = (np.cos(radii) * np.cos(3 * angles)).flatten()

# Create a custom triangulation.
triang = tri.Triangulation(x, y)

# Mask off unwanted triangles.
triang.set_mask(np.hypot(x[triang.triangles].mean(axis=1),
                         y[triang.triangles].mean(axis=1))
                < min_radius)

obj = plot3d()
obj.axs.tricontourf(triang, z, cmap=plt.cm.CMRmap)
obj.SavePng_Serial()

obj.axs.view_init(elev=45.)
obj.SavePng_Serial()

obj.axs.view_init(elev=-45.)
obj.SavePng_Serial()

obj.new_3Dfig()
obj.axs.plot_trisurf(x, y, z, cmap=plt.cm.CMRmap)
obj.SavePng_Serial()
