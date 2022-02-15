import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import sys
import os
import time
import cv2
from optparse import OptionParser
from PIL import Image

sys.path.append(os.path.join("./"))
from base import plot2d, plot3d

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

if __name__ == '__main__':
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    parser.add_option("--pic", dest="pic", default="./pic.bmp")
    parser.add_option("--pxyz", dest="pxyz",
                      default=[0.0, 0.0, 0.0], type="float", nargs=3)
    opt, argc = parser.parse_args(argvs)
    print(opt, argc)

    # pip install opencv-python
    # pip install opencv-contrib-python

    img = cv2.imread(opt.pic, 1)

    # OpenCV: BGR
    # Pillow: RBG

    #img_src = np.asfarray(img, dtype='uint8')
    img_gry = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    #img_gry = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)
    img_dst = cv2.Canny(img_gry, 0, 300)
    #dat = np.asfarray(img)

    obj = plot2d(aspect="equal")
    obj.axs.imshow(img_gry, cmap="jet")

    obj.new_2Dfig()
    obj.axs.imshow(img_dst, cmap="jet")
    obj.SavePng()

    img = Image.open(obj.tempname + ".png")
    img.show()
