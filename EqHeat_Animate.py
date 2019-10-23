""" From "COMPUTATIONAL PHYSICS", 3rd Ed, Enlarged Python eTextBook  
    by RH Landau, MJ Paez, and CC Bordeianu
    Copyright Wiley-VCH Verlag GmbH & Co. KGaA, Berlin;  Copyright R Landau,
    Oregon State Unv, MJ Paez, Univ Antioquia, C Bordeianu, Univ Bucharest, 2015.
    Support by National Science Foundation"""

# EqHeat.py Animated heat equation soltn via fine differences

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d

from base import plot3d

obj = plot3d()
obj.axs.set_xlabel('Position')
obj.axs.set_ylabel('time')
obj.axs.set_zlabel('Temperature')

# Parameters
Nx = 101
Nt = 3000
Dx = 0.01414
Dt = 1.
# Thermal conductivity
KAPPA = 210.
# Specific heat
SPH = 900.
# Density
RHO = 2700.

# Constant combo in algorthim
cons = KAPPA / (SPH * RHO) * Dt / Dx**2

# Temp @ first 2 times
T = np.zeros((Nx, 2), float)
T[:, 0] = 100.0

# Ends of bar at T = 0
T[0, 0] = 0.0
T[0, 1] = 0.
T[-1, 0] = 0.
T[-1, 1] = 0.0

for i in range(0, Nx - 1):
    # Scaled x's
    tempe.x[i] = 2.0 * i - 100.0
    # Scaled y's (Temp)
    tempe.y[i] = 0.8 * T[i, 0] - 20.0

img = []
for k in range(10):
    for ix in range(1, Nx - 1):
        T[ix, 1] = T[ix, 0] + cons * \
            (T[ix + 1, 0] + T[ix - 1, 0] - 2. * T[ix, 0])

    for i in range(0, Nx):
        # Scale 0<x<100 -> -100<x<100
        tempe.x[i] = 2.0 * i - 100.0
        tempe.y[i] = 0.6 * T[i, 1] - 20.0

    for ix in range(1, Nx - 1):
        # Row of 100 positions at t = m
        T[ix, 0] = T[ix, 1]
