import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button

from base import plot2d
from GaussPlot import gauss_2d, gauss_2d_skew


class ShowSection(plot2d):

    def __init__(self):
        plot2d.__init__(self)
        self.div_axs()
        self.gauss_func()

        self.get_mn([50, 10])
        self.im = self.axs.contourf(*self.mesh, self.func, cmap="jet")
        self.pl_x, = self.ax_x.plot(
            self.mesh[0][self.mx, :], self.func[self.mx, :])
        self.pl_y, = self.ax_y.plot(
            self.func[:, self.my], self.mesh[1][:, self.my])
        self.hx = self.axs.axhline(y=self.sy, xmin=self.xs, xmax=self.xe)
        self.hx.set_color("k")
        self.vy = self.axs.axvline(x=self.sx, ymin=self.ys, ymax=self.ye)
        self.vy.set_color("k")
        plt.colorbar(self.im, ax=self.axs, shrink=0.9)

        self.fig.canvas.mpl_connect('button_press_event', self.onclick)

    def replot_contourf(self, loc=[0, 0]):
        print(*loc)
        self.get_mn(loc)
        self.pl_x.set_xdata(self.mesh[0][self.mx, :])
        self.pl_x.set_ydata(self.func[self.mx, :])
        self.pl_y.set_xdata(self.func[:, self.my])
        self.pl_y.set_ydata(self.mesh[1][:, self.my])
        self.hx.set_ydata(self.sy)
        self.vy.set_xdata(self.sx)
        plt.draw()

    def onclick(self, event):
        if event.button == 1 and event.xdata is not None:
            print(event)
            sx = event.xdata
            sy = event.ydata
            self.replot_contourf([sx, sy])

    def gauss_func(self):
        nx, ny = 500, 500
        lx, ly = 200, 150
        wx, wy = 40, 25
        sx, sy = 50, 10
        kx, ky = 1, 1.0
        deg = 0

        px = np.linspace(-1, 1, nx) * lx
        py = np.linspace(-1, 1, ny) * ly
        self.mesh = np.meshgrid(px, py)
        peak = gauss_2d_skew(self.mesh, [0., sy], [5., 5.], [kx, ky], deg=0)
        data = gauss_2d_skew(self.mesh, [sx, sy], [wx, wy], [kx, ky], deg=30)
        self.func = data + peak

    def get_mn(self, loc=[0, 0]):
        sx, sy = loc
        nx, ny = self.func.shape
        xs, ys = self.mesh[0][0, 0], self.mesh[1][0, 0]
        xe, ye = self.mesh[0][0, -1], self.mesh[1][-1, 0]
        dx = self.mesh[0][0, 1] - self.mesh[0][0, 0]
        dy = self.mesh[1][1, 0] - self.mesh[1][0, 0]
        self.sx, self.sy = sx, sy
        self.xs, self.ys = xs, ys
        self.xe, self.ye = xe, ye
        self.mx = int((sy - ys) / dy)
        self.my = int((sx - xs) / dx)

    def make_button(self, name="Next", pos=[0.7, 0.05, 0.1, 0.075]):
        axs_pos = plt.axes(pos)
        axs_bot = Button(axs_pos, name)
        return axs_bot


if __name__ == "__main__":
    obj = ShowSection()
    plt.show()
