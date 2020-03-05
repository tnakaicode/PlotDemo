""" From "COMPUTATIONAL PHYSICS", 3rd Ed, Enlarged Python eTextBook  
    by RH Landau, MJ Paez, and CC Bordeianu
    Copyright Wiley-VCH Verlag GmbH & Co. KGaA, Berlin;  Copyright R Landau,
    Oregon State Unv, MJ Paez, Univ Antioquia, C Bordeianu, Univ Bucharest, 2015.
    Support by National Science Foundation"""

# Beam.py: solves Navier-Stokes equation for the flow around a beam

import numpy as np
import matplotlib.pyplot as plt
import sys
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import image

from base import plot2d

# Grid parameters
Nxmax = 100
Nymax = 20

# Stream
u = np.zeros((Nxmax + 1, Nymax + 1), float)
# Vorticity
w = np.zeros((Nxmax + 1, Nymax + 1), float)

# Initial v
V0 = 1.0
# Relaxation param
omega = 0.1
# Geometry
IL = 10

H = 8
T = 8
h = 1.
# Viscosity
nu = 1.

# Reynold number, normal units
R = V0 * h / nu
print("Working, wait for the figure, count to 30")


def borders():
    # Initialize stream,vorticity, sets BC
    for i in range(0, Nxmax + 1):
        for j in range(0, Nymax + 1):
            w[i, j] = 0.
            u[i, j] = j * V0

    # Fluid surface
    for i in range(0, Nxmax + 1):
        u[i, Nymax] = u[i, Nymax - 1] + V0 * h
        w[i, Nymax - 1] = 0.

    # Inlet
    for j in range(0, Nymax + 1):
        u[1, j] = u[0, j]
        w[0, j] = 0.

    # Centerline
    for i in range(0, Nxmax + 1):
        if i <= IL and i >= IL + T:
            u[i, 0] = 0.
            w[i, 0] = 0.
    # Outlet
    for j in range(1, Nymax):
        w[Nxmax, j] = w[Nxmax - 1, j]
        u[Nxmax, j] = u[Nxmax - 1, j]  # Borders


def beam():
    # BC for the beam
    # Beam sides
    for j in range(0, H + 1):
        # Front side
        w[IL, j] = -2 * u[IL - 1, j] / h**2
        # Back side
        w[IL + T, j] = -2 * u[IL + T + 1, j] / h**2
    for i in range(IL, IL + T + 1):
        # Top
        w[i, H - 1] = -2 * u[i, H] / h**2
    for i in range(IL, IL + T + 1):
        for j in range(0, H + 1):
            # Front
            u[IL, j] = 0.

            # Back
            u[IL + T, j] = 0.

            # Top
            u[i, H] = 0


def relax():
    # Method to relax stream
    # Reset conditions at beam
    beam()

    # Relax stream function
    for i in range(1, Nxmax):
        for j in range(1, Nymax):
            ui = u[i + 1, j] + u[i - 1, j]
            uj = u[i, j + 1] + u[i, j - 1]
            r1 = omega * ((ui + uj + h**2 * w[i, j]) / 4 - u[i, j])
            u[i, j] += r1

    # Relax vorticity
    for i in range(1, Nxmax):
        for j in range(1, Nymax):
            a1 = w[i + 1, j] + w[i - 1, j] + w[i, j + 1] + w[i, j - 1]
            a2 = (u[i, j + 1] - u[i, j - 1]) * (w[i + 1, j] - w[i - 1, j])
            a3 = (u[i + 1, j] - u[i - 1, j]) * (w[i, j + 1] - w[i, j - 1])
            r2 = omega * ((a1 - (R / 4.) * (a2 - a3)) / 4.0 - w[i, j])
            w[i, j] += r2


m = 0
i = 0
borders()
while (i <= 300):
    sys.stdout.write("\r{:d} / {:d}".format(i, 300))
    sys.stdout.flush()
    i += 1
    relax()


for i in range(0, Nxmax + 1):
    for j in range(0, Nymax + 1):
        u[i, j] = u[i, j] / (V0 * h)  # stream in V0h units

obj = plot2d()
img = obj.axs.contourf(u, origin='lower', cmap="jet")
obj.axs.set_title('Stream function - 2D Flow over a beam')
obj.fig.colorbar(img)
obj.fig.savefig(obj.tmpdir + obj.rootname + "-Stream.png")

obj.new_fig()
img = obj.axs.contourf(w, origin='lower', cmap="jet")
obj.axs.set_title('Vorticity - 2D Flow over a beam')
obj.fig.colorbar(img)
obj.fig.savefig(obj.tmpdir + obj.rootname + "-Velocity.png")
plt.show()
