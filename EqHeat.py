""" From "COMPUTATIONAL PHYSICS", 3rd Ed, Enlarged Python eTextBook  
    by RH Landau, MJ Paez, and CC Bordeianu
    Copyright Wiley-VCH Verlag GmbH & Co. KGaA, Berlin;  Copyright R Landau,
    Oregon State Unv, MJ Paez, Univ Antioquia, C Bordeianu, Univ Bucharest, 2015.
    Support by National Science Foundation"""

# EqHeat.py: solves heat equation via finite differences, 3-D plot

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from base import plot3d

Nx = 101
Nt = 3000
Dx = 0.03
Dt = 0.9
KAPPA = 210.
SPH = 900.
RHO = 2700.  # Conductivity, specf heat, density
T = np.zeros((Nx, 2), float)
Tpl = np.zeros((Nx, 31), float)

print("Working, wait for figure after count to 10")

for ix in range(1, Nx - 1):
    T[ix, 0] = 100.0
    # Initial T
T[0, 0] = 0.0
T[0, 1] = 0.
# 1st & last T = 0
T[Nx - 1, 0] = 0.
T[Nx - 1, 1] = 0.0
# constant
cons = KAPPA / (SPH * RHO) * Dt / (Dx * Dx)
m = 1

for t in range(1, Nt):
    for ix in range(1, Nx - 1):
        T[ix, 1] = T[ix, 0] + cons * \
            (T[ix + 1, 0] + T[ix - 1, 0] - 2. * T[ix, 0])
    if t % 300 == 0 or t == 1:
        # Every 300 steps
        for ix in range(1, Nx - 1, 2):
            Tpl[ix, m] = T[ix, 1]
        print(m)
        m = m + 1
    for ix in range(1, Nx - 1):
        T[ix, 0] = T[ix, 1]
x = list(range(1, Nx - 1, 2))
y = list(range(1, 30))
X, Y = np.meshgrid(x, y)


def functz(Tpl):
    z = Tpl[X, Y]
    return z


Z = functz(Tpl)

obj = plot3d()
obj.axs.plot_wireframe(X, Y, Z, color='r')
obj.axs.set_xlabel('Position')
obj.axs.set_ylabel('time')
obj.axs.set_zlabel('Temperature')
plt.show()
