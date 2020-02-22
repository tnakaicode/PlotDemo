import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib import animation
from mpl_toolkits.axes_grid1 import make_axes_locatable

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.HLRBRep import HLRBRep_PolyAlgo, HLRBRep_AreaLimit

from base import plotocc

if __name__ == '__main__':
    obj = plotocc()
    obj.show_axs_pln(scale=25)

    axs = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1))
    obj.show_box(axs, lxyz=[100, 50, 100])

    axs = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 1))
    obj.show_box(axs, lxyz=[100, 50, 100])

    axs = gp_Ax3(gp_Pnt(-100, 0, 0), gp_Dir(0, 1, 1))
    obj.show_axs_pln(axs, scale=25)

    axs = gp_Ax3(gp_Pnt(100, 0, 0), gp_Dir(0, 1, 1))
    obj.show_ellipsoid(axs, rxyz=[100, 30, 50])

    obj.show()
