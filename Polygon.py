import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

from OCC.BRepBuilderAPI import (BRepBuilderAPI_MakeFace,
                                BRepBuilderAPI_MakePolygon,
                                BRepBuilderAPI_MakeShell,
                                BRepBuilderAPI_MakeSolid)
from OCC.GeomAPI import GeomAPI_Interpolate, GeomAPI_PointsToBSpline
from OCC.GeomLib import GeomLib_Interpolate
from OCC.Geom import Geom_BSplineCurve
from OCC.gp import gp_Ax1, gp_Ax2, gp_Ax3, gp_Dir, gp_Pnt, gp_Vec
from OCC.TColgp import TColgp_Array1OfPnt, TColgp_HArray1OfPnt
from OCC.TColStd import TColStd_Array1OfReal
from OCCUtils.Construct import make_edge, make_face, make_polygon
from scipy.spatial import ConvexHull

from base import plotocc


""" Spatial Transformations
.. autosummary::
   :toctree: generated/

   KDTree      -- class for efficient nearest-neighbor queries
   cKDTree     -- class for efficient nearest-neighbor queries (faster impl.)
   Rectangle
"""

""" Delaunay Triangulation, Convex Hulls and Voronoi Diagrams
.. autosummary::
   :toctree: generated/

   Delaunay    -- compute Delaunay triangulation of input points
   ConvexHull  -- compute a convex hull for input points
   Voronoi     -- compute a Voronoi diagram hull from input points
   SphericalVoronoi -- compute a Voronoi diagram from input points on the surface of a sphere
   HalfspaceIntersection -- compute the intersection points of input halfspaces
"""


class GenPolygon (plotocc):

    def __init__(self):
        plotocc.__init__(self)
        self.gen_area()

        for idx in self.cov.simplices:
            self.tri_plane(idx)

        for i, xyz in enumerate(self.pnt):
            p = gp_Pnt(*xyz)
            self.display.DisplayShape(p)

        p0 = gp_Pnt(*self.pnt[0])
        p1 = gp_Pnt(*self.pnt[1])
        p2 = gp_Pnt(*self.pnt[2])
        p3 = gp_Pnt(*self.pnt[3])

        poly = make_polygon([p0, p1, p2, p3], closed=True)
        #face = make_face(poly)
        # self.display.DisplayShape(poly)
        self.tri_curve()

    def gen_area(self, num=10):
        self.pnt = np.random.rand(num, 3)
        self.cov = ConvexHull(self.pnt[:, 0:3])

    def tri_plane(self, idx=[0, 1, 2]):
        p0 = gp_Pnt(*self.pnt[idx[0]])
        p1 = gp_Pnt(*self.pnt[idx[1]])
        p2 = gp_Pnt(*self.pnt[idx[2]])

        poly = make_polygon([p0, p1, p2], closed=True)
        face = make_face(poly)
        print(idx, poly, face)
        self.display.DisplayShape(face, transparency=0.7, color="BLUE")

    def tri_curve(self, idx=[0, 1, 2]):
        p0 = gp_Pnt(*self.pnt[idx[0]])
        p1 = gp_Pnt(*self.pnt[idx[1]])
        p2 = gp_Pnt(*self.pnt[idx[2]])
        #pts = TColgp_Array1OfPnt(1, 5)
        pts = TColgp_Array1OfPnt(1, 4)
        pts.SetValue(1, p0)
        pts.SetValue(2, p1)
        pts.SetValue(3, p2)
        pts.SetValue(4, p0)
        #pts.SetValue(5, p1)
        wgt = TColStd_Array1OfReal(1, 4)
        wgt.SetValue(1, 0)
        wgt.SetValue(2, 0)
        wgt.SetValue(3, 0)
        wgt.SetValue(4, 0)

        api = GeomLib_Interpolate(0, 4, pts, wgt)
        crv = api.Curve()
        print(crv)
        #self.display.DisplayShape(crv)


if __name__ == '__main__':
    obj = GenPolygon()
    obj.show_axs_pln(scale=0.5)
    obj.show()
