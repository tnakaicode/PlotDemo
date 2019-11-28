#! /usr/bin/env python
#

import numpy as np
import matplotlib.pyplot as plt
import time
from mpl_toolkits.mplot3d import Axes3D
from sys import exit

from OCC.Coregp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Coregp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Coregp import gp_Lin
from OCC.CoreGeom import Geom_Line, Geom_Surface
from OCC.CoreBRep import BRep_Tool_Surface
from OCC.CoreGeomAPI import GeomAPI_IntCS
from OCC.CoreGeomLProp import GeomLProp_SurfaceTool
from OCCUtils.Topology import Topo

from base import plot2d, plot3d, plotocc
from base import gen_ellipsoid, set_trf


class BrownMotion (plotocc):

    def __init__(self, axs=gp_Ax3()):
        plotocc.__init__(self)
        self.rxyz = [10., 11., 9.]
        self.r_min = np.min(self.rxyz)
        self.r_max = np.max(self.rxyz)

        self.axis = axs
        self.trsf = set_trf(ax2=self.axis)
        self.ellp = gen_ellipsoid(axs=self.axis, rxyz=self.rxyz)
        self.get_surf()

        lin = Geom_Line(self.axis.Location(), gp_Dir(0, 0, 1))
        self.int_surf = GeomAPI_IntCS(lin.GetHandle(), self.ellp_surf)

        self.dm = 3
        self.get_pnt_ellipsoid(self.rand_pnt())
        self.get_pnt_ellipsoid(self.rand_pnt())
        self.get_pnt_ellipsoid(self.rand_pnt())
        
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
        h_lin = Geom_Line(gp_Ax1(gp_Pnt(), gp_Dir(vec))).GetHandle()
        self.int_surf.Perform(h_lin, self.ellp_surf)
        u, v, w = self.int_surf.Parameters(1)
        GeomLProp_SurfaceTool().Value(self.ellp_surf, u, v, pnt)
        print(pnt)
        print(self.int_surf.Parameters(1))
        print(self.int_surf.Parameters(2))
        self.display.DisplayShape(pnt)

    def plot_occ(self):
        self.show_pnt([self.rxyz[0], 0, 0])
        self.show_pnt([0, self.rxyz[1], 0])
        self.show_pnt([0, 0, self.rxyz[2]])
        self.display.DisplayShape(self.ellp, transparency=0.5, color="BLUE")
        print(*self.rand_point())
        self.show_axs_pln(scale=self.r_min)
        self.show()


if __name__ == '__main__':
    obj = BrownMotion()
    obj.plot_occ()
