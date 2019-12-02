import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button

from base import plot2d


class PipePlot(plot2d):

    def __init__(self):
        plot2d.__init__(self)
        px = np.linspace(0, 50, 100)
        self.axs.plot(px, self.logistic_func(
            px,
            l0=5.0,
            rx=[2.5, 15.0, 40.0],
            lx=[6.0, 4.0, 6.0],
            gx=[5.0, 2.0, 2.0]
        ))
        self.axs.set_aspect('auto')
        self.axs.set_ylim(0.0, 10.0)

    def logistic_func(self, px, l0=10.0, lx=[11.0], rx=[0.1], gx=[0.1]):
        val = l0
        for i, v in enumerate(lx):
            idx = i + 1
            if i == 0:
                coef = lx[i] - l0
            else:
                coef = lx[i] - lx[i - 1]
            val += coef / (1 + np.exp(-gx[i] * (px - rx[i])))
        return val


if __name__ == "__main__":
    obj = PipePlot()
    obj.Show()
