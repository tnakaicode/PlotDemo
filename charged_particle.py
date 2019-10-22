import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode

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
from GaussPlot import gauss_1d, gauss_1d_skew


class PlotParticle (Particle, plotocc):

    def __init__(self):
        Particle.__init__(self)
        plotocc.__init__(self)
        self.init_dat = [0, 0, 0, 2, 1, 1]

        pz = np.linspace(0, 200, 100)
        fz = self.bz([1 / np.sqrt(2), 1 / np.sqrt(2), pz])
        print(fz.max())

    def set_solver(self):
        self.solver = ode(self.my_method).set_integrator('dopri5')
        self.solver.set_initial_value(self.init_dat, self.t0)

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
        v1 = -alpha * u + self.ey(dat[0:3])
        w1 = w
        return np.array([u, v, w, u1, v1, w1])

    def ey(self, xyz):
        x, y, z = xyz
        return 10

    def bz(self, xyz):
        x, y, z = xyz
        fz = gauss_1d_skew(z, sx=150, wx=50, kx=1.0) / 5.0
        #fz = fz / fz.max()
        radi = np.sqrt(x**2 + y**2)
        return fz / (0.1 * radi)

    def plot_bz(self):
        pz = np.linspace(0, 200, 100)

        plt.figure()
        plt.plot(pz, self.bz([10, 10, pz]))
        plt.plot(pz, self.bz([20, 20, pz]))
        plt.plot(pz, self.bz([30, 30, pz]))
        plt.plot(pz, self.bz([40, 40, pz]))
        plt.show()

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

    def run_until_z(self, pz=500):
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
    obj.q = -10

    obj.init_dat = [10, 5, 0, 5, 1, 2.0]
    obj.set_solver()
    # obj.run()

    obj.init_dat = [10, 5, 0, 5, 1, 2.0]
    obj.set_solver()
    obj.run_until_z()

    obj.init_dat = [10, 5, 0, 10, 1, 2.0]
    obj.set_solver()
    obj.run_until_z()

    obj.init_dat = [20, 5, 0, 10, 10, 2.0]
    obj.set_solver()
    obj.run_until_z()
    obj.show()

    # obj.plot_bz()
