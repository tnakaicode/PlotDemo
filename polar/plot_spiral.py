import numpy as np
import matplotlib.pyplot as plt
import json
import sys
import time
import os
import glob
import shutil
import datetime
import argparse

sys.path.append(os.path.join('../'))
from base import plot2d


class SpiralPlot (plot2d):

    def __init__(self):
        plot2d.__init__(self)
        self.ratio = 1.0
        self.xy0 = [0, 0]
        self.xy1 = [1, 0]

        self.plot_xy()
        for i in range(25):
            self.update_vertex()
            self.plot_xy()

    def plot_xy(self):
        px = [self.xy0[0], self.xy1[0]]
        py = [self.xy0[1], self.xy1[1]]
        self.axs.plot(px, py)
        self.axs.plot([0, self.xy1[0]], [0, self.xy1[1]])

    def update_vertex(self):
        x, y = self.xy1
        h = (x**2 + y**2) * self.ratio
        self.xy0 = self.xy1.copy()
        self.xy1 = [x - y / h, y + x / h]
        print(h, self.xy0, self.xy1)


if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--flag", dest="flag", default=1, type=int)
    opt = parser.parse_args()
    print(opt, argvs)

    obj = SpiralPlot()
    obj.SavePng()
