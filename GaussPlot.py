import numpy as np
import matplotlib.pyplot as plt
import json
import glob
import sys
import time
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
from linecache import getline, clearcache
from scipy.integrate import simps
from scipy.special import erf

from Plot2D import plot_contour_sub, plot_contour_xyz


def gauss_1d(px, sx=0, wx=10):
    py = np.exp(-0.5 * ((px - sx) / wx)**2)
    return py


def cdf_1d(px, sx=0, wx=10, kx=2):
    return (1 + erf(kx * (px - sx) / wx / np.sqrt(wx))) / 2


def gauss_1d_skew(px, sx=0, wx=10, kx=2):
    py = gauss_1d(px, sx, wx)
    py *= cdf_1d(px, sx, wx, kx)
    return py / py.max()


def gauss_2d(mesh, sxy=[0, 0], wxy=[10, 10], deg=0.0):
    x, y = mesh[0] - sxy[0], mesh[1] - sxy[1]
    rot = np.deg2rad(deg)
    px = x * np.cos(rot) - y * np.sin(rot)
    py = y * np.cos(rot) + x * np.sin(rot)
    fx = np.exp(-0.5 * (px / wxy[0])**2)
    fy = np.exp(-0.5 * (py / wxy[1])**2)
    return fx * fy


def gauss_2d_skew(mesh, sxy=[0, 0], wxy=[10, 10], kxy=[2, 2], deg=0.0):
    x, y = mesh[0] - sxy[0], mesh[1] - sxy[1]
    rot = np.deg2rad(deg)
    px = x * np.cos(rot) - y * np.sin(rot)
    py = y * np.cos(rot) + x * np.sin(rot)
    fx = gauss_1d_skew(px, 0, wxy[0], kxy[0])
    fy = gauss_1d_skew(py, 0, wxy[1], kxy[1])
    return fx * fy


if __name__ == "__main__":
    argvs = sys.argv

    nx, ny = 500, 500
    lx, ly = 200, 150
    wx, wy = 40, 25
    sx, sy = 50, 10
    kx, ky = 7.5, 5.0
    deg = 0

    px = np.linspace(-1, 1, nx) * lx
    py = np.linspace(-1, 1, ny) * ly
    mesh = np.meshgrid(px, py)
    func = gauss_2d_skew(mesh, [sx, sy], [wx, wy], [kx, ky], deg=30)
    peak = gauss_2d_skew(mesh, [sx, sy], [1.0, 2.0], [kx, ky], deg=0)

    plot_contour_sub(mesh, func, loc=[sx, sy], pngfile="../tmp/gauss")
    plot_contour_sub(mesh, func + peak,
                     loc=[sx, sy], pngfile="../tmp/gauss_peak")
