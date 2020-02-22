import math
import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

from base import LineDrawer

"""
mpl_connect

* 'button_press_event'
* 'button_release_event'
* 'draw_event'
* 'key_press_event'
* 'key_release_event'
* 'motion_notify_event'
* 'pick_event'
* 'resize_event'
* 'scroll_event'
* 'figure_enter_event',
* 'figure_leave_event',
* 'axes_enter_event',
* 'axes_leave_event'
* 'close_event'
"""


class GetCursor(LineDrawer):

    def __init__(self, dirname="./tmp/", txtname="plot_data"):
        LineDrawer.__init__(self, dirname, txtname)
        px = np.linspace(-1, 1, 100) * 500
        self.mesh = np.meshgrid(px, px)

    def run(self):
        self.run_base()
        self.axs.contourf(*self.mesh, self.mesh[0], cmap="jet")

        #fig.canvas.mpl_connect('motion_notify_event', self.onclick)
        self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        self.fig.canvas.mpl_connect('key_press_event', self.onkey)

        animation.FuncAnimation(
            self.fig, self.anim_animate, init_func=self.anim_init, frames=30, interval=100, blit=True)
        #self.anim = anim

    def check_traj(self):
        self.init_fig()
        self.axs.contourf(*self.mesh, self.mesh[0])
        for i in range(self.fg):
            dat_file = self.txtname + "-{:d}.txt".format(i)
            dat = np.loadtxt(dat_file, skiprows=1)
            self.axs.plot(dat[:, 1], dat[:, 2])


if __name__ == '__main__':
    ld = GetCursor()
    ld.run()
    ld.show()
    ld.check_traj()
    ld.show()
