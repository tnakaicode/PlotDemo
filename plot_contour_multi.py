import numpy as np
import matplotlib.pyplot as plt
import os
import sys

sys.path.append(os.path.join("../"))
from base import SetDir, plot2d

if __name__ == "__main__":
    obj = plot2d()
    obj.new_2Dfig()

    px = np.arctan(np.linspace(-1, 1, 100)) * 100 + 150
    py = np.arctan(np.linspace(-1, 1, 200)) * 100 + 150
    mesh = np.meshgrid(px, py)
    func = mesh[0]**2 / 1000
    obj.axs.contourf(*mesh, func, cmap="jet")

    px = np.linspace(-1, 1, 200) * 100 - 150
    py = np.linspace(-1, 1, 100) * 100 - 150
    mesh = np.meshgrid(px, py)
    func = mesh[1]**2 / 1000
    obj.axs.contourf(*mesh, func, cmap="jet")
    obj.SavePng_Serial()
