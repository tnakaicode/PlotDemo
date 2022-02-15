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


def line1(p1, q1, p2, q2, ds):
    # conversion of point coordinates
    a = p2 - p1
    b = q2 - q1
    if 0 < b:
        dx = -ds / np.sqrt(1 + a**2 / b**2)
        dy = -a / b * dx
    if b < 0:
        dx = ds / np.sqrt(1 + a**2 / b**2)
        dy = -a / b * dx
    if b == 0 and 0 < a:
        dx = 0
        dy = ds
    if b == 0 and a < 0:
        dx = 0
        dy = -ds
    xx1 = p1 + dx
    yy1 = q1 + dy
    xx2 = p2 + dx
    yy2 = q2 + dy
    return xx1, yy1, xx2, yy2


def refb(xb, yb, ds):
    # reference boundary
    xc = np.zeros(len(xb), dtype=np.float64)
    yc = np.zeros(len(yb), dtype=np.float64)
    for i in range(len(xb)):
        if i == 0:
            x1 = xb[-1]
            y1 = yb[-1]
            x2 = xb[i]
            y2 = yb[i]
            x3 = xb[i + 1]
            y3 = yb[i + 1]
        if 1 <= i < len(xb) - 1:
            x1 = xb[i - 1]
            y1 = yb[i - 1]
            x2 = xb[i]
            y2 = yb[i]
            x3 = xb[i + 1]
            y3 = yb[i + 1]
        if i == len(xb) - 1:
            x1 = xb[i - 1]
            y1 = yb[i - 1]
            x2 = xb[i]
            y2 = yb[i]
            x3 = xb[0]
            y3 = yb[0]
        xx1a, yy1a, xx2a, yy2a = line1(x1, y1, x2, y2, ds)
        xx1b, yy1b, xx2b, yy2b = line1(x2, y2, x3, y3, ds)
        if xx2a - xx1a == 0 and xx2b - xx1b != 0:
            aa2 = (yy2b - yy1b) / (xx2b - xx1b)
            bb2 = (xx2b * yy1b - xx1b * yy2b) / (xx2b - xx1b)
            xc[i] = xx2a
            yc[i] = aa2 * xc[i] + bb2
        if xx2a - xx1a != 0 and xx2b - xx1b == 0:
            aa1 = (yy2a - yy1a) / (xx2a - xx1a)
            bb1 = (xx2a * yy1a - xx1a * yy2a) / (xx2a - xx1a)
            xc[i] = xx2b
            yc[i] = aa1 * xc[i] + bb1
        if xx2a - xx1a != 0 and xx2b - xx1b != 0:
            aa1 = (yy2a - yy1a) / (xx2a - xx1a)
            bb1 = (xx2a * yy1a - xx1a * yy2a) / (xx2a - xx1a)
            aa2 = (yy2b - yy1b) / (xx2b - xx1b)
            bb2 = (xx2b * yy1b - xx1b * yy2b) / (xx2b - xx1b)
            xc[i] = (bb2 - bb1) / (aa1 - aa2)
            yc[i] = (aa1 * bb2 - bb1 * aa2) / (aa1 - aa2)
    return xc, yc


def drawfig(xx, yy, zz, xb, yb, xc, yc):
    # drawing
    plt.figure(figsize=(5, 5), facecolor='w')
    plt.grid(color='#999999', linestyle='solid')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.plot(np.append(xb, xb[0]), np.append(
        yb, yb[0]), '-', color='#000000', lw=1)
    plt.plot(np.append(xc, xc[0]), np.append(
        yc, yc[0]), '--', color='#999999', lw=1)
    plt.contourf(xx, yy, zz, cmap='jet')
    fnameF = 'fig_test4.png'
    plt.savefig(fnameF, dpi=100, bbox_inches="tight", pad_inches=0.1)


if __name__ == '__main__':
    obj = plot2d()

    start = time.time()
    ds = 0.2  # interval of grid for interpolation
    # boundary
    xb = np.array([1, 3, 3, 2, 1, 1, 3, 2, 6, 5, 6, 7, 7, 5, 7])
    yb = np.array([1, 2, 3, 3, 2, 4, 6, 7, 7, 5, 5, 4, 3, 3, 1])
    xc, yc = refb(xb, yb, ds * 0.5)
    bxy = np.zeros((len(xc), 2), dtype=np.float64)
    for i in range(len(xc)):
        bxy[i, 0] = xc[i]
        bxy[i, 1] = yc[i]
    poly = Polygon(*bxy)

    # data of points
    xmin, xmax = 0, 8
    ymin, ymax = 0, 8
    n = 200
    np.random.seed(seed=31)
    x = xmin + (xmax - xmin) * np.random.rand(n)
    y = ymin + (ymax - ymin) * np.random.rand(n)
    z = y

    # making meshgrid and interpolation
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
    obj.axs.plot(
        np.append(xc, xc[0]), np.append(yc, yc[0]),
        '--', color='#999999', lw=1)
    obj.axs.contourf(xx, yy, zz, cmap='jet')
    obj.SavePng()
