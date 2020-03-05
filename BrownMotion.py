#! /usr/bin/env python
#

import numpy as np
import matplotlib.pyplot as plt
import time
from mpl_toolkits.mplot3d import Axes3D
from sys import exit

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Lin
from OCC.Core.Geom import Geom_Line, Geom_Surface
from OCC.Core.BRep import BRep_Tool_Surface
from OCC.Core.GeomAPI import GeomAPI_IntCS
from OCC.Core.GeomLProp import GeomLProp_SurfaceTool
from OCCUtils.Topology import Topo

from base import plot2d, plot3d, plotocc
from base import gen_ellipsoid, set_trf


class BrownMotion (object):

    def __init__(self, axs=gp_Ax3()):
        self.dm = 3
        self.nx = 200
        self.px = np.zeros([self.dm, self.nx])
        self.pt = np.linspace(0, 1, self.nx)
        self.rxyz = [10., 11., 9.]
        self.r_min = np.min(self.rxyz)
        self.r_max = np.max(self.rxyz)

        self.axis = axs
        self.trsf = set_trf(ax2=self.axis)
        self.ellp = gen_ellipsoid(axs=self.axis, rxyz=self.rxyz)
        self.get_surf()

        lin = Geom_Line(self.axis.Location(), gp_Dir(0, 0, 1))
        self.int_surf = GeomAPI_IntCS(lin, self.ellp_surf)

        self.pnt0 = gp_Pnt(0, 0, self.rxyz[2])

        p0 = self.rand_pnt()
        p1 = self.rand_pnt()
        self.get_pnt_ellipsoid(p0)

    def get_surf(self):
        top_api = Topo(self.ellp)
        print(top_api.number_of_faces())
        for face in top_api.faces():
            ellp_face = face
        self.ellp_surf = BRep_Tool_Surface(ellp_face)

    def rand_point(self):
        dx = np.random.randn(self.dm)
        dx_n = np.sqrt(np.sum(dx ** 2))
        x, y, z = dx / dx_n
        p, q, r = self.rxyz
        return x * p, y * q, z * r

    def rand_point_next(self):
        dx = np.random.uniform(-1, 1, self.dm)
        u = dx[0] * 2 * np.pi
        v = dx[1] * np.pi
        r = dx[2] / 0.5
        dx_n = np.sqrt(np.sum(dx ** 2))
        x, y, z = dx / dx_n
        p, q, r = self.rxyz
        return x * p, y * q, z * r

    def rand_pnt(self):
        dx = np.random.randn(self.dm)
        dx_n = np.sqrt(np.sum(dx ** 2))
        return gp_Pnt(*dx / dx_n)

    def get_pnt_ellipsoid(self, pnt=gp_Pnt(0, 0, 1)):
        vec = gp_Vec(gp_Pnt(), pnt)
        h_lin = Geom_Line(gp_Ax1(gp_Pnt(), gp_Dir(vec)))
        self.int_surf.Perform(h_lin, self.ellp_surf)
        u, v, w = self.int_surf.Parameters(2)
        GeomLProp_SurfaceTool().Value(self.ellp_surf, u, v, gp_Pnt())
        print(pnt)
        print(self.int_surf.Parameters(1))
        print(self.int_surf.Parameters(2))

    def plot_plt3d(self):
        obj = plot3d()
        obj.plot_ball(rxyz=self.rxyz)
        obj.axs.scatter(*self.rand_point())
        obj.axs.scatter(*self.rand_point())
        #obj.plot_ball(rxyz=[1, 2, 1])
        print(*np.random.randn(self.dm))
        print(*np.random.randn(self.dm))
        print(*np.random.randn(self.dm))
        print(*np.random.randn(self.dm))
        plt.show()

    def plot_occ(self):
        obj = plotocc()
        obj.show_pnt([self.rxyz[0], 0, 0])
        obj.show_pnt([0, self.rxyz[1], 0])
        obj.show_pnt([0, 0, self.rxyz[2]])
        obj.display.DisplayShape(self.ellp, transparency=0.5, color="BLUE")
        print(*self.rand_point())
        for t in self.pt:
            obj.show_pnt(self.rand_point())
        #obj.show_ball(scale=self.rxyz[0], trans=0.9)
        #obj.show_ball(scale=self.rxyz[1], trans=0.8)
        #obj.show_ball(scale=self.rxyz[2], trans=0.7)
        obj.show_axs_pln(scale=self.r_min)
        obj.display.View.Dump(obj.tmpdir + "BrownMotion.png")
        obj.show()


if __name__ == '__main__':
    obj = BrownMotion()
    obj.plot_occ()
