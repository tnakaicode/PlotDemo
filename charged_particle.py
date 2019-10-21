import numpy as np
import matplotlib.pyplot as plt

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Lin
from OCCUtils.Construct import make_box, make_line
from OCCUtils.Construct import make_plane, make_polygon
from OCCUtils.Construct import make_wire, make_edge
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import dir_to_vec, vec_to_dir

from base import plot3d, plotocc
from ChargedParticle import Particle


class PlotParticle (Particle, plotocc):

    def __init__(self):
        Particle.__init__(self)
        plotocc.__init__(self)
        self.ini_dat = [0, 0, 0, 2, 1, 1]

    def run(self):
        self.t0 = 0.0
        self.t1 = 5.0
        self.dt = 0.01
        self.pos = []

        while self.solver.successful() and self.solver.t < self.t1:
            print(self.solver.t, *self.solver.y)
            self.solver.integrate(self.solver.t + self.dt)
            self.pos.append(self.solver.y[:3])

        pts = []
        for i, val in enumerate(self.pos[:-1]):
            p0 = gp_Pnt(*val[0:3])
            p1 = gp_Pnt(*self.pos[i + 1][0:3])
            pts.append(p0)
            self.display.DisplayShape(make_edge(p0, p1))


if __name__ == '__main__':
    obj = PlotParticle()

    obj.init_dat = [0, 0, 0, 2, 0, 0]
    obj.init_solver()
    obj.run()

    obj.init_dat = [0, 0, 0, 2, 1, 1]
    obj.init_solver()
    obj.run()

    obj.init_dat = [0, 0, 0, 4, 1, 1]
    obj.init_solver()
    obj.run()
    obj.show()
