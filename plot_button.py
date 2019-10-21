import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button

from base import plot2d


class Index(plot2d):

    def __init__(self):
        plot2d.__init__(self)
        self.div_axs()

        self.ind = 0
        self.freqs = np.arange(2, 20, 3)
        self.pt = np.arange(0.0, 1.0, 0.001)
        self.ps = np.sin(2 * np.pi * self.freqs[0] * self.pt)
        self.pl, = self.axs.plot(self.pt, self.ps)

        # plt.subplots_adjust(bottom=0.2)

    def next(self, event):
        self.ind += 1
        i = self.ind % len(self.freqs)
        ydata = np.sin(2 * np.pi * self.freqs[i] * self.pt)
        self.pl.set_ydata(ydata)
        #self.axs.plot(self.pt, ydata)
        plt.draw()

    def prev(self, event):
        self.ind -= 1
        i = self.ind % len(self.freqs)
        ydata = np.sin(2 * np.pi * self.freqs[i] * self.pt)
        self.pl.set_ydata(ydata)
        plt.draw()

    def make_button(self, name="Next", pos=[0.7, 0.05, 0.1, 0.075]):
        axs_pos = plt.axes(pos)
        axs_bot = Button(axs_pos, name)
        return axs_bot


if __name__ == "__main__":
    obj = Index()

    bnext = obj.make_button("Next", [0.7, 0.05, 0.1, 0.075])
    bnext.on_clicked(obj.next)

    bprev = obj.make_button("Previous", [0.81, 0.05, 0.1, 0.075])
    bprev.on_clicked(obj.prev)
    plt.show()
