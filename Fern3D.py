
""" From "COMPUTATIONAL PHYSICS", 3rd Ed, Enlarged Python eTextBook  
    by RH Landau, MJ Paez, and CC Bordeianu
    Copyright Wiley-VCH Verlag GmbH & Co. KGaA, Berlin;  Copyright R Landau,
    Oregon State Unv, MJ Paez, Univ Antioquia, C Bordeianu, Univ Bucharest, 2015.
    Support by National Science Foundation"""

# Fern3D.py:  Fern in 3D, see Barnsley, "Fractals Everywhere"

import numpy as np
import matplotlib.pyplot as plt
import random
import sys

from base import plot3d, plotocc

if __name__ == "__main__":
    imax = 5*10**4
    x = 0.5
    y = 0.0
    z = -0.2
    xn = 0.0
    yn = 0.0

    graph1 = plotocc()
    for i in range(1, imax):
        sys.stdout.write("\r {:d} / {:d}".format(i, imax))
        sys.stdout.flush()
        r = random.random()
        # random number
        if (r <= 0.1):
            # 10% probability
            xn = 0.0
            yn = 0.18 * y
            zn = 0.0
        elif (r > 0.1 and r <= 0.7):
            # 60% probability
            xn = 0.85 * x
            yn = 0.85 * y + 0.1 * z + 1.6
            zn = -0.1 * y + 0.85 * z
            # print xn,yn,zn
        elif (r > 0.7 and r <= 0.85):
            # 15 % probability
            xn = 0.2 * x - 0.2 * y
            yn = 0.2 * x + 0.2 * y + 0.8
            zn = 0.3 * z
        else:
            # 15% probability
            xn = -0.2 * x + 0.2 * y
            yn = 0.2 * x + 0.2 * y + 0.8
            zn = 0.3 * z
        x = xn
        y = yn
        z = zn

        # linear TF for plot
        xc = 4.0 * x
        yc = 2.0 * y - 7
        zc = z
        graph1.show_pnt([xc, yc, zc])

    print("\n")
    graph1.show()
