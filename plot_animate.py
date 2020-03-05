import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from base import plot2d


class AniFunc (plot2d):

    def __init__(self):
        plot2d.__init__(self)

        self.xdata = []
        self.ydata = []
        self.ln, = self.axs.plot([], [], 'ro')

        ani = FuncAnimation(
            self.fig, self.update,
            frames=np.linspace(0, 2 * np.pi, 128),
            init_func=self.init, blit=True
        )
        ani.save(self.tempname + ".gif", writer="pillow")

    def init(self):
        self.axs.set_xlim(0, 2 * np.pi)
        self.axs.set_ylim(-1, 1)
        return self.ln,

    def update(self, frame):
        self.xdata.append(frame)
        self.ydata.append(np.sin(frame))
        self.ln.set_data(self.xdata, self.ydata)
        return self.ln,


if __name__ == "__main__":
    obj = AniFunc()
    plt.show()
