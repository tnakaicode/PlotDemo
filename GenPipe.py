import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button

from base import plot2d


class PipePlot(plot2d):

    def __init__(self):
        plot2d.__init__(self)
        px = np.linspace(0, 1, 100)
        self.axs.plot(px, self.logistic_func(px, lxy=[10.0, 15.0]))
        self.axs.plot(px, self.logistic_func(px, lxy=[10.0, 9.0], x0=0.1))
        self.axs.set_aspect('auto')
        self.axs.set_ylim(0, 20)

    def logistic_func(self, px, x0=0.0, lxy=[10, 15.0], grow=1.0):
        val = (lxy[1] - lxy[0]) / (1 + np.exp(-grow * (px - x0))) + lxy[0]
        return val


if __name__ == "__main__":
    obj = PipePlot()
    obj.Show()
