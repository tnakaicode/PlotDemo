import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from scipy import interpolate

from base import plot2d


class GenBoundary (plot2d):

    def __init__(self, num=5, lxy=10):
        plot2d.__init__(self)
        self.xx = []
        self.yy = []
        self.kk = 3
        self.pn = num

        self.axs.plot([-lxy, lxy, lxy, -lxy, -lxy],
                      [-lxy, -lxy, lxy, lxy, -lxy])
        self.spl, = self.axs.plot(self.xx, self.yy)
        self.pts, = self.axs.plot(self.xx, self.yy, 'o')

        self.fig.canvas.mpl_connect('button_press_event', self._OnClick)
        self.fig.canvas.mpl_connect('key_press_event', self._KeyClick)

    def make_button(self, name="Next", pos=[0.7, 0.05, 0.1, 0.075]):
        pos_axs = self.fig.add_axes(pos)
        axs_bot = Button(pos_axs, name)
        return axs_bot

    def _KeyClick(self, event):
        if event.key == 'p':
            self._Spline(event)
        elif event.key == 'c':
            self._Clear(event)

    def _Clear(self, event):
        self.xx = []
        self.yy = []
        self.pn = 2
        self.spl.set_xdata(self.xx)
        self.spl.set_ydata(self.yy)
        self.pts.set_xdata(self.xx)
        self.pts.set_ydata(self.yy)
        plt.draw()

    def _Spline(self, event):
        num = len(self.xx)
        tck, u = interpolate.splprep(
            [self.xx + [self.xx[0]], self.yy + [self.yy[0]]], k=self.kk, s=0
        )
        pu = np.linspace(0, 1, num=self.pn * num)
        spline = interpolate.splev(pu, tck)
        self.spl.set_xdata(spline[0])
        self.spl.set_ydata(spline[1])
        plt.draw()

    def _OnClick(self, event):
        if event.dblclick:
            self._Clear()
        elif event.button == 1 and event.xdata is not None:
            print(event)
            self.xx.append(event.xdata)
            self.yy.append(event.ydata)
            self.pts.set_xdata(self.xx)
            self.pts.set_ydata(self.yy)
            plt.draw()


if __name__ == '__main__':
    obj = GenBoundary()
    obj.Show()

    print(*obj.xx)
    print(*obj.yy)
