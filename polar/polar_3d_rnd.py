import numpy as np
import sys
import time
import os

sys.path.append(os.path.join("../"))
from base import plot3d


def sample_spherical(npoints, ndim=3):
    vec = np.random.randn(ndim, npoints)
    vec /= np.linalg.norm(vec, axis=0)
    return vec


phi = np.linspace(0, np.pi, 20)
theta = np.linspace(0, 2 * np.pi, 40)
x = np.outer(np.sin(theta), np.cos(phi))
y = np.outer(np.sin(theta), np.sin(phi))
z = np.outer(np.cos(theta), np.ones_like(phi))

xi, yi, zi = sample_spherical(100)

obj = plot3d()
obj.axs.plot_wireframe(x, y, z, color='k', rstride=1, cstride=1)
obj.axs.scatter(xi, yi, zi, s=100, c='r', zorder=10)
obj.SavePng()
