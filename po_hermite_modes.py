"""
po_hermite_modes.py: calculates Gauss-Schell modes

see e.g. https://commons.wikimedia.org/wiki/File:Hermite-gaussian.png

"""

import numpy as np
import matplotlib.pyplot as plt
import json
import sys
import time
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.special import hermite

from Plot2D import plot_contour_sub
from GaussPlot import gauss_1d


def hermite_gauss(x, w, m=0):
    return (hermite(m)(np.sqrt(2) * x / w) * gauss_1d(x, sx=0, wx=w))


if __name__ == "__main__":
    size_x = 1.0
    size_y = .50
    n_x = 620
    n_y = 320
    w_x = size_x / 8.
    w_y = size_y / 8.
    m = 1
    n = 0

    X = np.linspace(-size_x / 2, size_x / 2, n_x)
    Y = np.linspace(-size_y / 2, size_y / 2, n_y)
    mesh = np.meshgrid(X, Y)

    for m in range(4):
        for n in range(4):
            outx = hermite_gauss(mesh[0], w_x, m)
            outy = hermite_gauss(mesh[1], w_y, n)
            func = outx * outy
            name = "TE_{:d}_{:d}".format(m, n)

            pngfile = "./tmp/po_hermite_modes_" + name
            plot_contour_sub(mesh, func, title=name, dirname=pngfile)
