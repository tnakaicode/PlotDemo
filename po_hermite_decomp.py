"""
po_hermite_modes.py: calculates Gauss-Schell modes

see e.g. https://commons.wikimedia.org/wiki/File:Hermite-gaussian.png

"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
from optparse import OptionParser
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.special import hermite
from scipy.integrate import simps

sys.path.append(os.path.join("./"))
from Plot2D import plot_contour_sub
from GaussPlot import gauss_1d


def integrate_simps(mesh, func):
    nx, ny = func.shape
    px, py = mesh[0][int(nx / 2), :], mesh[1][:, int(ny / 2)]
    val = simps(simps(func, px), py)
    return val


def hermite_gauss(x, w, m=0):
    return (hermite(m)(np.sqrt(2) * x / w) * gauss_1d(x, sx=0, wx=w))


def hermite_func(mesh, wxy=[0.1, 0.1], cfgtxt="./po_hermite.txt"):
    dat = np.loadtxt(cfgtxt, comments="#")
    func = np.zeros_like(mesh[0])
    for data in dat:
        m = int(data[0])
        n = int(data[1])
        v = float(data[2])
        func_x = hermite_gauss(mesh[0], wxy[0], m)
        func_y = hermite_gauss(mesh[1], wxy[1], n)
        func += v * (func_x * func_y)
    return func


if __name__ == "__main__":
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    parser.add_option("--lxy", dest="lxy",
                      default=[1.0, 0.5], type="float", nargs=2)
    parser.add_option("--dxy", dest="dxy",
                      default=[0.01, 0.01], type="float", nargs=2)
    parser.add_option("--wxy", dest="wxy",
                      default=[0.1, 0.1], type="float", nargs=2)
    parser.add_option("--nxy", dest="nxy",
                      default=[500, 600], type="int", nargs=2)
    opt, argc = parser.parse_args(argvs)
    print(opt, argc)

    nxy = opt.nxy
    lxy = opt.lxy
    dxy = opt.dxy
    wxy = opt.wxy
    dat = np.loadtxt("./po_hermite.txt", comments="#")

    px = np.linspace(-1, 1, nxy[0]) * lxy[0] + dxy[0]
    py = np.linspace(-1, 1, nxy[1]) * lxy[1] + dxy[1]
    mesh = np.meshgrid(px, py)
    func = hermite_func(mesh, wxy, cfgtxt="./po_hermite.txt")

    for m in range(4):
        for n in range(4):
            outx = hermite_gauss(mesh[0], wxy[0], m)
            outy = hermite_gauss(mesh[1], wxy[1], n)
            g_func = outx * outy
            val = integrate_simps(mesh, func * g_func)
            txt = "{:d}\t{:d}\t{:.5f}".format(m, n, val)
            print(txt)
