""" From "COMPUTATIONAL PHYSICS", 3rd Ed, Enlarged Python eTextBook  
    by RH Landau, MJ Paez, and CC Bordeianu
    Copyright Wiley-VCH Verlag GmbH & Co. KGaA, Berlin;  Copyright R Landau,
    Oregon State Unv, MJ Paez, Univ Antioquia, C Bordeianu, Univ Bucharest, 2015.
    Support by National Science Foundation"""

# FDTD.py  FDTD solution of Maxwell's equations in 1-D

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from base import plot2d

# scene = display(x=0,y=0,width= 800, height= 500, \
#         title= 'E: cyan, H: red. Periodic BC',forward=(-0.6,-0.5,-1))
# Efield = curve(x=list(range(0,xmax)),color=color.cyan,radius=1.5,display=scene)
# Hfield = curve(x=list(range(0,xmax)),color=color.red, radius=1.5,display=scene)
# vplane= curve(pos=[(-xmax,ymax),(xmax,ymax),(xmax,-ymax),(-xmax,-ymax),
#                    (-xmax,ymax)],color=color.cyan)
# zaxis=curve(pos=[(-xmax,0),(xmax,0)],color=color.magenta)
# hplane=curve(pos=[(-xmax,0,zmax),(xmax,0,zmax),(xmax,0,-zmax),(-xmax,0,-zmax),
#                    (-xmax,0,zmax)],color=color.magenta)
# ball1 = sphere(pos = (xmax+30, 0,0), color = color.black, radius = 2)


class FDTD1d (plot2d):

    def __init__(self):
        plot2d.__init__(self)
        self.xmax = 201
        self.ymax = 100
        self.zmax = 100

        self.ts = 2
        self.beta = 0.01
        self.Ex = np.zeros((self.xmax, self.ts), float)
        self.Hy = np.zeros((self.xmax, self.ts), float)
        self.ti = 0

        k = np.arange(self.xmax)
        self.Ex[:self.xmax, 0] = 0.1 * np.sin(2 * np.pi * k / 100.0)
        self.Hy[:self.xmax, 0] = 0.1 * np.sin(2 * np.pi * k / 100.0)

        self.kn, = self.axs.plot(self.Ex[:, 0], self.Ex[:, 1], 'ro')
        self.ln, = self.axs.plot(self.Hy[:, 0], self.Hy[:, 1], 'bo')
        self.dt = 0

        self.ani = FuncAnimation(
            self.fig, self.update, frames=np.linspace(0, 2 * np.pi, 128)
        )
        self.ani.save(self.tmpdir + "Sample.gif", writer='imagemagick')

    def update(self, frame):
        self.dt += 10
        self.Ex[1:self.xmax - 1, 1] = self.Ex[1:self.xmax - 1, 0] + \
            self.beta * (self.Hy[0:self.xmax - 2, 0] - self.Hy[2:self.xmax, 0])
        self.Hy[1:self.xmax - 1, 1] = self.Hy[1:self.xmax - 1, 0] + \
            self.beta * (self.Ex[0:self.xmax - 2, 0] - self.Ex[2:self.xmax, 0])
        self.Ex[0, 1] = self.Ex[0, 0] + self.beta * \
            (self.Hy[self.xmax - 2, 0] - self.Hy[1, 0])
        self.Ex[self.xmax - 1, 1] = self.Ex[self.xmax - 1, 0] + \
            self.beta * (self.Hy[self.xmax - 2, 0] - self.Hy[1, 0])
        self.Hy[0, 1] = self.Hy[0, 0] + self.beta * \
            (self.Ex[self.xmax - 2, 0] - self.Ex[1, 0])
        self.Hy[self.xmax - 1, 1] = self.Hy[self.xmax - 1, 0] + \
            self.beta * (self.Ex[self.xmax - 2, 0] - self.Ex[1, 0])

        # self.xdata.append(frame)
        # self.ydata.append(np.sin(frame))
        self.kn.set_data(self.Ex[:, 0], self.Ex[:, 1])
        self.ln.set_data(self.Hy[:, 0], self.Hy[:, 1])

        self.Ex[:self.xmax, 0] = self.Ex[:self.xmax, 1]
        self.Hy[:self.xmax, 0] = self.Hy[:self.xmax, 1]


if __name__ == "__main__":
    obj = FDTD1d()
