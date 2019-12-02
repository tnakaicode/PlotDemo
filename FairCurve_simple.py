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

# FairCurve
#
# Constructs the two contact points P1 and P2
# The geometrical characteristics of the batten (elastic beam)
#
# These include the real number values for
#   height of deformation Height,
#   slope value Slope
#   kind of energy PhysicalRatio.
#       Jerk (0)
#       Sagging (1).
#
# Note that the default setting for Physical Ration are initialized
#   FreeSliding = False
#   ConstraintOrder1 = 1
#   ConstraintOrder2 = 1
#   Angle1 = 0
#   Angle2 = 0
#   Curvature1 = 0
#   Curvature2 = 0
#   SlidingFactor = 1
#
# Warning If PhysicalRatio equals 1,
#   you cannot impose constraints on curvature.
#
# Exceptions NegativeValue if Height is less than or equal to 0.
#
# NullValue if the distance between P1 and P2 is less than or equal to
#   the tolerance value for distance in Precision::Confusion: P1.IsEqual(P2, Precision::Confusion()).
#
# The function gp_Pnt2d::IsEqual tests to see if this is the case.
#
# Definition of the geometricals constraints


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
        self.fc_MinVar.SetSlidingFactor(0.1)
        self.fc_MinVar.SetPhysicalRatio(0.5)
        stat = self.fc_MinVar.Compute()

        for i in range(0, 20):
            # TODO: the parameter slope needs to be visualized
            slope = i / 20.
            self.fc_MinVar.SetAngle1(np.deg2rad(i))
            self.fc_MinVar.SetAngle2(np.deg2rad(-i))
            self.fc_MinVar.SetSlope(slope)
            txt = "\r {}".format(self.fc_MinVar.DumpToString())
            # print(self.fc_MinVar.DumpToString())
            self.fc_MinVar.Compute()
            # self.display.EraseAll()
            self.display.DisplayShape(self.fc_MinVar.Curve(), update=True)
            sys.stdout.write(txt)
            sys.stdout.flush()
            print()
            # time.sleep(0.21)
        self.display.DisplayShape(
            self.fc_MinVar.Curve(), update=True, color="BLUE")


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
