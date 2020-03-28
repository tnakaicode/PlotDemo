import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from scipy.spatial import ConvexHull

from base import plotocc

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.GeomAPI import GeomAPI_Interpolate, GeomAPI_PointsToBSpline
from OCC.Core.TColStd import TColStd_Array1OfReal
from OCC.Core.TColgp import TColgp_Array1OfPnt, TColgp_Array2OfPnt
from OCC.Core.TColgp import TColgp_HArray1OfPnt, TColgp_HArray2OfPnt
from OCCUtils.Construct import make_plane, make_polygon
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import dir_to_vec, vec_to_dir


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


class ConvexArea (plotocc):

    def __init__(self, num=30, idx=2):
        plotocc.__init__(self)
        self.pnt = np.random.rand(num, 3)
        self.cov = ConvexHull(self.pnt[:, 0:2])

        for xyz in self.pnt:
            print(xyz)
            self.show_pnt(xyz)

        # print(self.cov)
        print(self.cov.vertices, len(self.cov.vertices))
        h_pts = TColgp_Array1OfPnt(1, len(self.cov.vertices) + 1)
        h_par = TColStd_Array1OfReal(1, len(self.cov.vertices) + 1)
        for i, idx in enumerate(self.cov.vertices):
            x, y, z = self.pnt[idx]
            self.show_pnt([x, y, 0])
            h_pts.SetValue(i + 1, gp_Pnt(x, y, 0))
            h_par.SetValue(i + 1, 1.0)
            print(i + 1, idx, x, y)
        h_pts.SetValue(h_pts.Upper(), h_pts.Value(1))

        api = GeomAPI_PointsToBSpline(h_pts)

        print(self.cov.vertices)
        self.display.DisplayShape(api.Curve())


if __name__ == '__main__':
    obj = ConvexArea(num=50)
    obj.SaveMenu()
    obj.show()
