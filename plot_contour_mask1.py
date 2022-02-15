from sympy.geometry import Point, Polygon
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import time
import sys
import os

sys.path.append(os.path.join("./"))
from base import plot2d, plot3d

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)


def drawfig(xx, yy, zz, xb, yb):
    # drawing
    plt.figure(figsize=(5, 5), facecolor='w')
    plt.grid(color='#999999', linestyle='solid')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.plot(np.append(xb, xb[0]), np.append(
        yb, yb[0]), '-', color='#000000', lw=1)
    plt.contourf(xx, yy, zz, cmap='jet')
    fnameF = 'fig_test3.png'
    plt.savefig(fnameF, dpi=100, bbox_inches="tight", pad_inches=0.1)


if __name__ == '__main__':
    obj = plot2d()

    start = time.time()

    # boundary
    xb = np.array([1, 3, 3, 2, 1, 1, 3, 2, 6, 5, 6, 7, 7, 5, 7])
    yb = np.array([1, 2, 3, 3, 2, 4, 6, 7, 7, 5, 5, 4, 3, 3, 1])
    bxy = np.zeros((len(xb), 2), dtype=np.float64)
    for i in range(len(xb)):
        bxy[i, 0] = xb[i]
        bxy[i, 1] = yb[i]
    poly = Polygon(*bxy)

    # data of points
    xmin, xmax = 0, 8
    ymin, ymax = 0, 8
    n = 200
    np.random.seed(seed=31)
    x = xmin + (xmax - xmin) * np.random.rand(n)
    y = ymin + (ymax - ymin) * np.random.rand(n)
    z = y

    # makeing mrshgrid and interpolation
    ds = 0.2
    x1 = np.arange(xmin, xmax + ds, ds)
    y1 = np.arange(ymin, ymax + ds, ds)
    xx, yy = np.meshgrid(x1, y1)
    zz = interpolate.griddata((x, y), z, (xx, yy), method='nearest')
    #zz = interpolate.griddata((x,y), z, (xx, yy), method='linear')
    #zz = interpolate.griddata((x,y), z, (xx, yy), method='cubic')

    # judgement (True or False)
    # In case of point beeing just on the boundary, result is 'False'
    xp = xx[0, :]
    yp = yy[:, 0]
    ii = []
    jj = []
    for i in range(len(xp)):
        for j in range(len(yp)):
            s = poly.encloses_point(Point(xp[i], yp[j]))
            if s == False:
                ii = ii + [i]
                jj = jj + [j]
    zz[jj, ii] = np.nan

    print(time.time() - start)

    obj.new_2Dfig()
    obj.axs.plot(
        np.append(xb, xb[0]), np.append(yb, yb[0]),
        '-', color='#000000', lw=1)
    obj.axs.contourf(xx, yy, zz, cmap='jet')
    obj.SavePng()
