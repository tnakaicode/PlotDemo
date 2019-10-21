import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib import animation
from mpl_toolkits.axes_grid1 import make_axes_locatable

from OCC.Display.SimpleGui import init_display
from OCC.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.gp import gp_Ax1, gp_Ax2, gp_Ax3

if __name__ == '__main__':
    display, start_display, add_menu, add_functionto_menu = init_display()

    display.DisplayShape(gp_Pnt())

    display.FitAll()
    start_display()
