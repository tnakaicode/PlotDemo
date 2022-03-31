import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import imageio
import cv2
import scipy.misc
from PIL import Image
from PIL.Image import fromarray as toimage
import argparse

sys.path.append(os.path.join("./"))
from base import plot2d

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--pxyz", dest="point",
                      default=[0.0, 0.0, 0.0], type=float, nargs=3)
    opt = parser.parse_args()
    print(opt, argvs)

    obj = plot2d()

    src_img = cv2.imread("./img/lena.png")

    # parameter
    sg = [1, 1, 1, 200, 300, 300]
    bb = [1, 1, 200, 200]

    # working file
    mask = np.zeros_like(src_img)  # (y, x, c)

    # segmentation data
    sg = np.asarray(sg)
    poly_number = int(len(sg) / 2)
    poly = np.zeros((poly_number, 2))
    for i in range(poly_number):
        poly[i][0] = sg[(i * 2) + 0]  # x
        poly[i][1] = sg[(i * 2) + 1]  # y

    # generate mask
    mask = cv2.fillConvexPoly(
        mask, np.array(poly, 'int32'),
        color=(255, 255, 255))

    # generate src_img and mask_image
    cv2.imwrite(obj.tempname + "-001.png", src_img)
    cv2.imwrite(obj.tempname + "-002.png", mask)

    # masked image
    masked_src_img = np.where(mask == 255, src_img, mask)
    cv2.imwrite(obj.tempname + "-003.png", masked_src_img)

    # cropping img
    bb = np.asarray(bb, 'int32')
    offset_x = bb[0]
    offset_y = bb[1]
    length_x = bb[2]
    length_y = bb[3]
    cv2.imwrite(obj.tempname + "-004.png",
                masked_src_img[offset_y: (offset_y + length_y), offset_x: (offset_x + length_x)])
