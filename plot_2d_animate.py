import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse
from matplotlib import animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D

from base import SetDir, plot2d


def data_gen(t=0):
    cnt = 0
    while cnt < 1000:
        cnt += 1
        if cnt % 100 == 0:
            print(cnt)
        t += 0.05
        yield t, np.sin(2 * np.pi * t) * np.exp(-t / 10.)


class Decay (plot2d):

    def __init__(self):
        plot2d.__init__(self)
        self.axs.set_aspect("auto")
        self.xdata = []
        self.ydata = []
        self.line, = self.axs.plot(self.xdata, self.ydata, lw=2)

        self.ani = animation.FuncAnimation(
            self.fig, self.run, data_gen, blit=False, interval=0.1, repeat=False, init_func=self.init)

    def run(self, data):
        # update the data
        t, y = data
        self.xdata.append(t)
        self.ydata.append(y)
        xmin, xmax = self.axs.get_xlim()

        if t >= xmax:
            self.axs.set_xlim(xmin, 1.25 * xmax)
            self.axs.figure.canvas.draw()
        self.line.set_data(self.xdata, self.ydata)
        return self.line,

    def init(self):
        self.axs.set_ylim(-1.1, 1.1)
        self.axs.set_xlim(0, 10)
        del self.xdata[:]
        del self.ydata[:]
        self.line.set_data(self.xdata, self.ydata)
        return self.line,

    def SaveGif(self, savename=None):
        if savename == None:
            savename = self.tmpdir + self.rootname + ".gif"
        self.ani.save(savename, writer="pillow")


if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    opt = parser.parse_args()
    print(opt, argvs)

    obj = Decay()
    obj.SaveGif()
    obj.SavePng()
