import numpy as np
import matplotlib.pyplot as plt
import sys
import time
import os
import argparse

sys.path.append(os.path.join('../'))
from base import plot2d, plotpolar

theta = np.arange(0.0, 4 * 2 * np.pi, 0.01)
r1 = 0.5 * theta
r2 = 0.2 * theta

obj = plotpolar()
obj.plot_polar(theta, r1, arrow=False, label="0.2")
obj.plot_polar(theta, r2, arrow=False, label="0.5")
plt.legend()
obj.SavePng_Serial()
