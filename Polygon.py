import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from scipy.spatial import ConvexHull

from base import plotocc

from OCC.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.gp import gp_Ax1, gp_Ax2, gp_Ax3

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
        for i, xyz in enumerate(self.pnt):
            p = gp_Pnt(*xyz)
            self.display.DisplayShape(p)

    def gen_area(self, num=30):
        self.pnt = np.random.rand(num, 3)
        self.cov = ConvexHull(self.pnt[:, 0:2])


if __name__ == '__main__':
    obj = GenPolygon()
    obj.show()
