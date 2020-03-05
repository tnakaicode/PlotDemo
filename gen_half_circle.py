import numpy as np
import matplotlib.pyplot as plt

from base import plot2d


class GenArea (plot2d):

    def __init__(self):
        plot2d.__init__(self)
        self.pnt = []
        self.xy0 = np.array([0, 0])
        self.pnt.append(self.xy0)

        #self.fig, self.axs = plt.subplots()
        # self.axs.set_aspect('equal')
        # self.axs.xaxis.grid()
        # self.axs.yaxis.grid()

    def gen_arc(self, sxy=[0, 0], radi=10, deg=[-30, 60], num=50):
        pd = np.linspace(*deg, num)
        pr = np.deg2rad(pd)
        for i, v in enumerate(pr):
            x = sxy[0] + radi * np.cos(v)
            y = sxy[0] + radi * np.sin(v)
            self.pnt.append(np.array([x, y]))

    def data_ouput(self):
        dat = np.array(obj.pnt)
        np.savetxt(self.tmpdir + "data_area.txt", dat)


if __name__ == '__main__':
    obj = GenArea()
    obj.gen_arc()
    obj.pnt.append(obj.xy0)

    dat = np.array(obj.pnt)
    obj.axs.plot(dat[:, 0], dat[:, 1])
    plt.savefig(obj.tmpdir + "data_area.png")
    plt.show()

    obj.data_ouput()
