import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from scipy.spatial import ConvexHull

from base import plotocc

from OCC.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCCUtils.Construct import make_edge

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

    def gen_area(self, num=30):
        self.pnt = np.random.rand(num, 3)
        self.cov = ConvexHull(self.pnt[:, 0:3])

    def tri_plane(self, idx=[0, 1, 2]):
        print(idx)
        p0 = gp_Pnt(*self.pnt[idx[0]])
        p1 = gp_Pnt(*self.pnt[idx[1]])
        p2 = gp_Pnt(*self.pnt[idx[2]])

        l0 = make_edge(p0, p1)
        l1 = make_edge(p1, p2)
        l2 = make_edge(p2, p0)
        self.display.DisplayShape(l0)
        self.display.DisplayShape(l1)
        self.display.DisplayShape(l2)


if __name__ == '__main__':
    obj = GenPolygon()
    obj.show_axs_pln(scale=0.5)
    obj.show()
