import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import imageio
import cv2
import csv
import scipy.misc
from PIL import Image
from PIL.Image import fromarray as toimage
import argparse

sys.path.append(os.path.join("./"))
from base import plot2d

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)


def sxy_to_nxy(mesh, sxy=[0, 0]):
    sx, sy = sxy
    nx, ny = mesh[0].shape
    xs, ys = mesh[0][0, 0], mesh[1][0, 0]
    dx, dy = mesh[0][0, 1] - mesh[0][0, 0], mesh[1][1, 0] - mesh[1][0, 0]
    mx, my = int((sy - ys) / dy), int((sx - xs) / dx)
    return [my, mx]


def rim_mask(mesh, rim_file, idx=0):
    rim = np.loadtxt(rim_file, skiprows=2)
    nxy = []
    for xy in rim:
        nxy.append(sxy_to_nxy(mesh, xy))
        print(xy, nxy[-1])
    (x, y, w, h) = cv2.boundingRect(np.array(nxy))
    print(x, y, w, h)

    pts = np.array(nxy)
    data = np.ones_like(mesh[0])
    mask = np.zeros_like(mesh[0].shape, np.uint8)
    mask = cv2.fillConvexPoly(data, np.array(pts, "int32"),
                              color=(255, 255, 255)) / 255
    if idx == 0:
        return mask
    else:
        mask_out = np.abs(1 - mask)
        return mask_out


if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--pxyz", dest="point",
                      default=[0.0, 0.0, 0.0], type=float, nargs=3)
    opt = parser.parse_args()
    print(opt, argvs)

    px = np.linspace(-1, 1, 400) * 200
    py = np.linspace(-1, 1, 400) * 200
    mesh = np.meshgrid(px, py)
    data = np.ones_like(mesh[0])
    rim_file = "plot_cv2_mask1_poly.rim"
    rim = np.loadtxt(rim_file, skiprows=2)
    mask = rim_mask(mesh, rim_file, idx=1)

    obj = plot2d(aspect="equal")
    obj.contourf_sub(mesh, data,
                     pngname=obj.tempname + "_data.png")
    obj.contourf_sub(mesh, data * mask,
                     pngname=obj.tempname + "_mask.png")

    obj.new_2Dfig(aspect="equal")
    obj.axs.plot(rim[:, 0], rim[:, 1])
    obj.SavePng(obj.tempname + "_rim.png")
