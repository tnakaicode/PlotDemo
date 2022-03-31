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

    src = np.array(Image.open('img/lena.png'))
    mask = np.array(Image.open(
        'img/horse.png').resize(src.shape[1::-1], Image.BILINEAR))

    print(mask.dtype, mask.min(), mask.max())
    # uint8 0 255

    mask = mask / 255

    print(mask.dtype, mask.min(), mask.max())
    # float64 0.0 1.0

    dst = src * mask

    Image.fromarray(dst.astype(np.uint8)).save(obj.tempname + '.jpg')
