#!/usr/bin/env python

# Copyright 2009-2015 Jelle Feringa (jelleferinga@gmail.com)
##
# This file is part of pythonOCC.
##
# pythonOCC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
##
# pythonOCC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
##
# You should have received a copy of the GNU Lesser General Public License
# along with pythonOCC.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import math
import time
import sys

from OCC.Core.gp import gp_Pnt2d, gp_Pln
from OCC.Core.Geom import Geom_Plane
from OCC.Core.FairCurve import FairCurve_MinimalVariation, FairCurve_Newton
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
from OCC.Core.FEmTool import FEmTool_ProfileMatrix
from OCC.Display.SimpleGui import init_display


from base import plotocc


class DrawCurve (plotocc):

    def __init__(self):
        plotocc.__init__(self)
        self.pt1 = gp_Pnt2d(0., 0.)
        self.pt2 = gp_Pnt2d(0., 120.)
        self.hit = 100.
        self.pln = Geom_Plane(gp_Pln())
        self.fc_MinVar = FairCurve_MinimalVariation(
            self.pt1, self.pt2, self.hit, 1)
        self.fc_MinVar.SetP1(self.pt1)
        self.fc_MinVar.SetP2(self.pt2)
        self.fc_MinVar.SetConstraintOrder1(2)
        self.fc_MinVar.SetConstraintOrder2(2)
        self.fc_MinVar.SetAngle1(0)
        self.fc_MinVar.SetAngle2(0)
        self.fc_MinVar.SetHeight(self.hit)
        self.fc_MinVar.SetSlope(0)
        self.fc_MinVar.SetFreeSliding(True)
        stat = self.fc_MinVar.Compute()

        for i in range(0, 40):
            # TODO: the parameter slope needs to be visualized
            slope = i / 100.
            self.fc_MinVar.SetAngle1(np.deg2rad(i))
            self.fc_MinVar.SetAngle2(np.deg2rad(-i))
            self.fc_MinVar.SetSlope(slope)
            print(self.fc_MinVar.DumpToString())
            self.fc_MinVar.Compute()
            self.display.EraseAll()
            edge = BRepBuilderAPI_MakeEdge(
                self.fc_MinVar.Curve(), self.pln).Edge()
            self.display.DisplayShape(edge, update=True)
            # time.sleep(0.21)

    def batten_curve(self, slope, angle1, angle2):
        fc = FairCurve_MinimalVariation(self.pt1, self.pt2, self.hit, slope)
        fc.SetConstraintOrder1(2)
        fc.SetConstraintOrder2(2)
        fc.SetAngle1(angle1)
        fc.SetAngle2(angle2)
        fc.SetHeight(self.hit)
        fc.SetSlope(slope)
        fc.SetFreeSliding(True)
        print(fc.DumpToString())
        status = fc.Compute()
        print(error_code(status[0]), error_code(status[1]))
        return fc.Curve()


def error_code(n):
    errors = {0: "FairCurve_OK",
              1: "FairCurve_NotConverged",
              2: "FairCurve_InfiniteSliding",
              3: "FairCurve_NullHeight",
              }
    return errors[n]


if __name__ == "__main__":
    obj = DrawCurve()
    obj.show()
