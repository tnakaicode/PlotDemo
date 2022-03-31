"""
po_hermite_modes.py: calculates Gauss-Schell modes

see e.g. https://commons.wikimedia.org/wiki/File:Hermite-gaussian.png

"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.special import hermite

sys.path.append(os.path.join("./"))
from Plot2D import plot_contour_sub
from GaussPlot import gauss_1d


def hermite_gauss(x, w, m=0):
    return (hermite(m)(np.sqrt(2) * x / w) * gauss_1d(x, sx=0, wx=w))


if __name__ == "__main__":
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--lxy", dest="lxy",
                      default=[1.0, 0.5], type=float, nargs=2)
    parser.add_argument("--dxy", dest="dxy",
                      default=[0.01, 0.01], type=float, nargs=2)
    parser.add_argument("--wxy", dest="wxy",
                      default=[0.1, 0.1], type=float, nargs=2)
    parser.add_argument("--nxy", dest="nxy",
                      default=[500, 600], type=int, nargs=2)
    opt = parser.parse_args()
    print(opt, argvs)

    nxy = opt.nxy
    lxy = opt.lxy
    dxy = opt.dxy
    wxy = opt.wxy
    dat = np.loadtxt("./po_hermite.txt", comments="#")

    px = np.linspace(-1, 1, nxy[0]) * lxy[0] + dxy[0]
    py = np.linspace(-1, 1, nxy[1]) * lxy[1] + dxy[1]
    mesh = np.meshgrid(px, py)
    func = np.zeros_like(mesh[0])

    for data in dat:
        print(data)
        m = int(data[0])
        n = int(data[1])
        v = float(data[2])
        func_x = hermite_gauss(mesh[0], wxy[0], m)
        func_y = hermite_gauss(mesh[1], wxy[1], n)
        func += v * (func_x * func_y)

    plot_contour_sub(mesh, func, dirname="./tmp/po_hermite")
