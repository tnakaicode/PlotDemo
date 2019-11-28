import numpy as np
import math

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeWire
from OCC.Core.TColgp import TColgp_Array2OfPnt
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface
from OCC.Core.GeomAbs import GeomAbs_C2
from OCC.Core.GeomAbs import GeomAbs_Intersection
from OCC.Core.BRepOffset import BRepOffset_MakeOffset, BRepOffset_Skin


def gyroid(x, y, z, t):
    return np.cos(x)*np.sin(y) + np.cos(y)*np.sin(z) + np.cos(z)*np.sin(x) - t


if __name__ == '__main__':
    display, start_display, add_menu, add_function_to_menu = init_display()

    lat = 2.0
    res = 11
    pt = 3.0

    x, y, z = np.pi/2. * np.mgrid[-1:1:res*1j, -1:1:res*1j, -1:1:res*1j] * lat
    vol = gyroid(x, y, z, pt)

    print(vol.shape[-1])
    for iz, vz in enumerate(vol[::, ]):
        print(iz)
        ixy = vz.shape
        pts = TColgp_Array2OfPnt(1, ixy[0], 1, ixy[1])
        for idx, val in np.ndenumerate(vz):
            ix, iy = idx[0], idx[1]
            px, py, pz = x[ix, iy, iz], y[ix, iy, iz], val
            pnt = gp_Pnt(px, py, pz)
            #display.DisplayShape(pnt)
            pts.SetValue(ix+1, iy+1, pnt)

        crve = GeomAPI_PointsToBSplineSurface(
            pts, 3, 8, GeomAbs_C2, 0.001).Surface()
        face = BRepBuilderAPI_MakeFace(crve, 1e-6).Face()
        display.DisplayShape(face)
    
    display.FitAll()
    start_display()
