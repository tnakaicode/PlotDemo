import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button

from base import plot2d


class PipePlot(plot2d):

    def __init__(self):
        plot2d.__init__(self)
        px = np.linspace(0, 10, 100)
        self.axs.plot(px, self.logistic_func(
            px, l0=10.0, rx=[2.5, 7.5], lx=[11.0, 9.0], gx=[5.0, 20.0]
        ))
        self.axs.set_aspect('auto')

    def logistic_func(self, px, l0=10.0, lx=[11.0], rx=[0.1], gx=[0.1]):
        val = l0
        for i, v in enumerate(lx):
            idx = i + 1
            if i == 0:
                val += (lx[i] - l0) / (1 + np.exp(-gx[i] * (px - rx[i])))
            else:
                val += (lx[i] - lx[i - 1]) / \
                    (1 + np.exp(-gx[i] * (px - rx[i])))
        return val


if __name__ == "__main__":
    obj = PipePlot()
    obj.Show()
