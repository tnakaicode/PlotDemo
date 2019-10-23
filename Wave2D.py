""" From "COMPUTATIONAL PHYSICS", 3rd Ed, Enlarged Python eTextBook  
    by RH Landau, MJ Paez, and CC Bordeianu
    Copyright Wiley-VCH Verlag GmbH & Co. KGaA, Berlin;  Copyright R Landau,
    Oregon State Unv, MJ Paez, Univ Antioquia, C Bordeianu, Univ Bucharest, 2015.
    Support by National Science Foundation"""

# Waves2Danal.py: analytical solution Helmholtz eqn rectangular membrane

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d

from base import plot3d

t = 0
c = np.sqrt(180. / 390.)   # speed, tension N/m2, density kg/m2
s5 = np.sqrt(5)
N = 32


def membrane(t, X, Y):
    return np.cos(c * s5 * t) * np.sin(2 * X) * np.sin(Y)


xs = np.linspace(0, np.pi, 32)
ys = np.linspace(0, np.pi, 32)  # 0->pi
X, Y = np.meshgrid(xs, ys)
Z = membrane(0, X, Y)    # x,y grid; init Z


obj = plot3d()
obj.axs.set_title('Vibrating Membrane')

img = []
for t in np.linspace(0, 10, 40):
    Z = membrane(t, X, Y)
    wframe = obj.axs.plot_wireframe(X, Y, Z)
    img.append([wframe])

ani = animation.ArtistAnimation(obj.fig, img, interval=100)
ani.save('./tmp/Wave2d.gif', writer='ffmpeg', dpi=100)
