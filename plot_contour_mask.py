import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import imageio
import cv2
import scipy.misc
from PIL.Image import fromarray as toimage
from optparse import OptionParser

sys.path.append(os.path.join("./"))
from base import plot2d

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

if __name__ == '__main__':
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    parser.add_option("--pxyz", dest="point",
                      default=[0.0, 0.0, 0.0], type="float", nargs=3)
    opt, argc = parser.parse_args(argvs)
    print(opt, argc)

    px = np.linspace(-1, 1, 100) * 200 + 50
    py = np.linspace(-1, 1, 200) * 150 - 50
    mesh = np.meshgrid(py, px)
    surf = mesh[0]**2 / 1000 + mesh[1]**2 / 400
    data = np.array(toimage(surf))

    pts = [
        [-100, -50],
        [-250, +275],
        [+150, +200],
        [+50, -100]
    ]
    pxy = [np.array(v) for v in pts]
    pxy.append(pxy[0])
    pxy = np.array(pxy)

    obj = plot2d()
    obj.contourf_div(mesh, surf)
    obj.axs.plot(pxy[:, 0], pxy[:, 1])
    obj.SavePng()

    nxy = []
    for xy in pxy:
        nxy.append(obj.sxy_to_nxy(mesh, xy))
        print(nxy[-1])

    original_frame = surf
    (x, y, w, h) = cv2.boundingRect(np.array(nxy))
    print(x, y, w, h)

    pts = np.array(nxy) - np.array(nxy).min(axis=0)
    mask = np.zeros(original_frame.shape, np.uint8)
    cv2.drawContours(mask, [pts], -1, (255, 255, 255), -1, cv2.LINE_AA)
    mask = cv2.fillConvexPoly(original_frame, np.array(
        pts, "int32"), color=(255, 255, 255))
    result = original_frame * mask
    res_data = np.asanyarray(result)
    print(res_data)
    cv2.imwrite(obj.tempname + "-xxx.png", result)

    obj.new_2Dfig()
    obj.contourf_div(mesh, mask, pngname=obj.tempname)
