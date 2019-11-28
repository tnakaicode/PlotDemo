import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from scipy.integrate import ode, complex_ode
sys.path.append(os.path.join('./'))

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Lin, gp_XYZ
from OCCUtils.Construct import make_box, make_line
from OCCUtils.Construct import make_plane, make_polygon
from OCCUtils.Construct import make_wire, make_edge
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import dir_to_vec, vec_to_dir

from base import plot3d, plotocc
from ChargedParticle import Particle
from GaussPlot import gauss_1d, gauss_1d_skew


class PlotParticle (Particle, plotocc):

    def __init__(self):
        Particle.__init__(self)
        plotocc.__init__(self)
        self.ini_dat = [0, 0, 0, 2, 1, 1]

    def set_solver(self):
        self.solver = ode(self.my_method).set_integrator('dop853')
        self.solver.set_initial_value(self.init_dat, self.t0)

        # vode
        # Real-valued Variable-coefficient Ordinary Differential Equation
        #
        # zvode
        # Complex-valued Variable-coefficient Ordinary Differential Equation
        #
        # lsoda
        # Real-valued Variable-coefficient Ordinary Differential Equation
        #
        # dopri5
        # Runge-Kutta method of order (4)5
        # due to Dormand Prince (with stepsize control and dense output).
        #
        # dop853
        # runge-kutta method of order 8(5,3)
        # due to Dormand Prince (with stepsize control and dense output).

    def my_method(self, t, dat):
        """
        Computes the derivative of the state vector y according to the equation of motion:
        Y is the state vector (x, y, z, u, v, w) === (position, velocity).
        returns dY/dt.
        """
        x, y, z = dat[:3]
        u, v, w = dat[3:]

        alpha = self.q / self.m * self.bz(dat[0:3])
        u1 = alpha * v
        v1 = -alpha * u
        w1 = w
        return np.array([u, v, w, u1, v1, w1])

    def ey(self, xyz):
        x, y, z = xyz
        return 1

    def bz(self, xyz):
        x, y, z = xyz
        radi = np.sqrt(x**2 + y**2)
        return gauss_1d_skew(z, sx=10, wx=10, kx=-2.0) / (0.1 * radi)

    def run(self):
        self.t0 = 0.0
        self.t1 = 10.0
        self.dt = 0.1
        self.pos = []

        while self.solver.successful() and self.solver.t < self.t1:
            txt = "{:.3f} ".format(self.solver.t)
            txt += "{:.3f} {:.3f} {:.3f} ".format(*self.solver.y[:3])
            txt += "{:.3f} {:.3f} {:.3f} ".format(*self.solver.y[3:])
            print(txt)
            self.solver.integrate(self.solver.t + self.dt)
            self.pos.append(self.solver.y[:3])
        self.show_pos()

    def run_until_z(self, pz=100):
        self.dt = 0.05
        self.pos = []

        while self.solver.successful() and self.solver.y[2] <= pz:
            txt = "{:.3f} ".format(self.solver.t)
            txt += "{:.3f} {:.3f} {:.3f} ".format(*self.solver.y[:3])
            txt += "{:.3f} {:.3f} {:.3f} ".format(*self.solver.y[3:])
            print(txt)
            self.solver.integrate(self.solver.t + self.dt)
            self.pos.append(self.solver.y)
        self.show_pos()

    def show_pos(self):
        pts = []
        for i, val in enumerate(self.pos[:-1]):
            p0 = gp_Pnt(*val[0:3])
            p1 = gp_Pnt(*self.pos[i + 1][0:3])
            pts.append(p0)
            self.display.DisplayShape(make_edge(p0, p1))


if __name__ == '__main__':
    obj = PlotParticle()
    obj.show_axs_pln(scale=10)
    obj.show_plane(scale=20.0)
    obj.q = -1

    axs = gp_Ax3(gp_Pnt(0, 0, 100), gp_Dir(0, 0, 1))
    #obj.show_plane(axs, scale=20.0)

    obj.init_dat = [10, 5, 0, 5, 1, 2.0]
    obj.set_solver()
    obj.run_until_z()

    obj.init_dat = [10, 5, 0, 10, 1, 2.0]
    obj.set_solver()
    # obj.run()

    obj.init_dat = [10, 5, -50, 10, 10, 2.0]
    obj.set_solver()
    obj.run_until_z()

    obj.init_dat = [20, 5, -50, 10, 10, 2.0]
    obj.set_solver()
    obj.run_until_z()
    obj.show()
