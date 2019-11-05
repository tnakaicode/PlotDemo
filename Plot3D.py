import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib import animation
from mpl_toolkits.axes_grid1 import make_axes_locatable

from OCC.Display.SimpleGui import init_display
from OCC.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.HLRBRep import HLRBRep_PolyAlgo, HLRBRep_AreaLimit

from base import plotocc

if __name__ == '__main__':
    obj = plotocc()
    obj.show_box()

    axs = gp_Ax3(gp_Pnt(-100, 0, 0), gp_Dir(0, 1, 1))
    obj.show_axs_pln(axs, scale=25)

    obj.show()
