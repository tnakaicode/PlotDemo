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
        self.int_surf = GeomAPI_IntCS(lin.GetHandle(), self.ellp_surf)

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
        h_lin = Geom_Line(gp_Ax1(gp_Pnt(), gp_Dir(vec))).GetHandle()
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
        obj.show()


def brownian_displacement_display(k, n, m, d, t, dsq):

    # *****************************************************************************80
    #
    # BROWNIAN_DISPLACEMENT_DISPLAY displays average Brownian motion displacement.
    #
    #  Discussion:
    #
    #    Thanks to Feifei Xu for pointing out a missing factor of 2 in the
    #    displacement calculation.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    10 June 2018
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Parameters:
    #
    #    Input, integer K, the number of repetitions.
    #
    #    Input, integer N, the number of time steps.
    #
    #    Input, integer M, the spatial dimension.
    #
    #    Input, real D, the diffusion coefficient.
    #
    #    Input, real T, the total time.
    #
    #    Input, real DSQ(K,N), the displacements over time for each repetition.
    #

    #  Get the T values.
    tvec = np.linspace(0, t, n)

    #  Select 5 random trajectories for display.
    for s in range(0, 5):
        i = int(k * np.random.rand(1))
        plt.plot(tvec, dsq[i, :], 'b-')

    #  Display the average displacement.
    dsq_ave = np.sum(dsq, 0) / float(k)
    plt.plot(tvec, dsq_ave, 'r-', linewidth=2)

    #  Display the ideal displacment.
    dsq_ideal = 2.0 * m * d * tvec
    plt.plot(tvec, dsq_ideal, 'k-', linewidth=3)

    plt.grid(True)
    plt.xlabel('<--T-->')
    plt.ylabel('<--D^2-->')
    plt.title('Squared displacement (Red), Predicted (Black), Samples (Blue)')

    filename = 'displacement_' + str(m) + '.png'
    plt.savefig(filename)
    plt.show()
    plt.clf()

    print('')
    print('  Plot saved as "%s".' % (filename))

    return


def brownian_displacement_simulation(k=20, n=1001, m=3, d=10.0, t=1.0):

    # *****************************************************************************80
    #
    # BROWNIAN_DISPLACEMENT_SIMULATION simulates Brownian displacement.
    #
    #  Discussion:
    #
    #    Thanks to Feifei Xu for pointing out a missing factor of 2 in the
    #    stepsize calculation, 08 March 2016.
    #
    #    This function computes the square of the distance of the Brownian
    #    particle from the starting point, repeating this calculation
    #    several times.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    10 June 2018
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Parameters:
    #
    #    Input, int  K, the number of repetitions.
    #
    #    Input, int  N, the number of time steps to take, plus 1.
    #
    #    Input, int  M, the spatial dimension.
    #
    #    Input, real D, the diffusion coefficient.  This might be 10.0.
    #    Computationally, this is simply a scale factor between time and space.
    #
    #    Input, real T, the total time.
    #
    #    Output, real DSQ(K,N), the displacements over time for each repetition.
    #    DSQ(:,1) is 0.0, because we include the displacement at the initial time.
    #

    dsq = np.zeros([k, n])

    for i in range(0, k):
        x = brownian_motion_simulation(m, n, d, t)
        dsq[i, 0:n] = np.sum(x[0:m, 0:n] ** 2, 0)

    return dsq


def brownian_motion_display(m, n, x):

    # *****************************************************************************80
    #
    # BROWNIAN_MOTION_DISPLAY displays successive Brownian motion positions.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    26 April 2018
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Parameters:
    #
    #    Input, integer M, the spatial dimension.
    #    M should be 1, 2 or 3.
    #
    #    Input, integer N, the number of time steps.
    #
    #    Input, real X(M,N), the particle positions.
    #

    if (m == 1):
        y = np.linspace(0, n - 1, n) / float(n - 1)
        plt.plot(x[0, :], y[:], 'b', linewidth=2)
        plt.plot(x[0, 0], y[0], 'g.', markersize=35)
        plt.plot(x[0, n - 1], y[n - 1], 'r.', markersize=35)
        plt.grid(True)
        plt.xlabel('<--X-->')
        plt.ylabel('<--Time-->')
        plt.title('Brownian motion simulation in 1D')
        filename = 'motion_' + str(m) + 'd.png'
        plt.savefig(filename)
        plt.show()
        plt.clf()

    elif (m == 2):
        plt.plot(x[0, :], x[1, :], 'b', LineWidth=2)
        plt.plot(x[0, 0], x[1, 0], 'g.', markersize=35)
        plt.plot(x[0, n - 1], x[1, n - 1], 'r.', markersize=35)
        plt.grid(True)
        plt.xlabel('<--X-->')
        plt.ylabel('<--Y-->')
        plt.title('Brownian motion simulation in 2D')
        plt.axis('equal')
        filename = 'motion_' + str(m) + 'd.png'
        plt.savefig(filename)
        plt.show()
        plt.clf()

    elif (m == 3 or 4):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(x[0, :], x[1, :], x[2, :], 'b', linewidth=2)
        ax.scatter(x[0, 0], x[1, 0], x[2, 0], c='g', marker='o', s=100)
        ax.scatter(x[0, n - 1], x[1, n - 1],
                   x[2, n - 1], c='r', marker='o', s=100)
        ax.grid(True)
        ax.set_xlabel('<--X-->')
        ax.set_ylabel('<--Y-->')
        ax.set_zlabel('<--Z-->')
        plt.title('Brownian motion simulation in 3D')
        #plt.axis ( 'equal' )
        filename = 'motion_' + str(m) + 'd.png'
        plt.savefig(filename)
        plt.show()
        plt.clf()

    else:
        print('')
        print('BROWNIAN_MOTION_DISPLAY - Fatal error!')
        print('  Cannot display data except for M = 1, 2, 3.')
        exit('BROWNIAN_MOTION_DISPLAY - Fatal error!')

    print('')
    print('  Plot saved as "%s".' % (filename))
    return


def brownian_motion_simulation(m=3, n=1001, d=10.0, t=1.0):

    # *****************************************************************************80
    #
    # BROWNIAN_MOTION_SIMULATION simulates Brownian motion.
    #
    #  Discussion:
    #
    #    Thanks to Feifei Xu for pointing out a missing factor of 2 in the
    #    stepsize calculation, 08 March 2016.
    #
    #    Thanks to Joerg Peter Pfannmoeller for pointing out a missing factor
    #    of M in the stepsize calculation, 23 April 2018.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    10 June 2018
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Parameters:
    #
    #    Input, int  M, the spatial dimension.
    #
    #    Input, int  N, the number of time steps to take, plus 1.
    #
    #    Input, real D, the diffusion coefficient.
    #
    #    Input, real T, the total time.
    #
    #    Output, real X(M,N), the initial position at time 0.0, and
    #    the N-1 successive locations of the particle.
    #

    #  Set the time step.
    dt = t / float(n - 1)

    #  Compute the individual steps.
    x = np.zeros([m, n])
    if (m >= 4):
        x[0, 0] = 0
        x[1, 0] = 0
        x[2, 0] = 1

    for j in range(1, n):
        #  S is the stepsize
        s = np.sqrt(2.0 * m * d * dt) * np.random.randn(1)

        #  Direction is random.
        if (m == 1):
            dx = np.random.randn(m)
        elif (m == 4):
            dx = np.random.randn(2)
            norm_dx = np.sqrt(np.sum(dx ** 2))
            phi = np.arctan2(x[1, j - 1], x[0, j - 1])
            rho = np.arccos(x[2, j - 1])
            dx[0] = phi + s * dx[0] / norm_dx
            dx[1] = rho + s * dx[1] / norm_dx
        elif (m == 5):
            dx = np.random.randn(3)
            norm_dx = np.sqrt(np.sum(dx ** 2))
            rad = np.sqrt(np.sum(x[0:3, j - 1] ** 2))
            phi = np.arctan2(x[1, j - 1] / rad, x[0, j - 1] / rad)
            rho = np.arccos(x[2, j - 1] / rad)
            dx[0] = phi + s * dx[0] / norm_dx
            dx[1] = rho + s * dx[1] / norm_dx
            dx[2] = 1 + s * dx[2] / norm_dx / 2
        else:
            dx = np.random.randn(m)
            norm_dx = np.sqrt(np.sum(dx ** 2))
            for i in range(0, m):
                dx[i] = s * dx[i] / norm_dx

        #  Each position is the sum of the previous steps.
        if (m == 4):
            x[0, j] = np.cos(dx[0]) * np.sin(dx[1])
            x[1, j] = np.sin(dx[0]) * np.sin(dx[1])
            x[2, j] = np.cos(dx[1])
        elif (m == 5):
            x[0, j] = dx[2] * np.cos(dx[0]) * np.sin(dx[1])
            x[1, j] = dx[2] * np.sin(dx[0]) * np.sin(dx[1])
            x[2, j] = dx[2] * np.cos(dx[1])
        else:
            x[0:m, j] = x[0:m, j - 1] + dx[0:m]
    return x


def brownian_motion_simulation_test():

    # *****************************************************************************80
    #
    # BROWNIAN_MOTION_SIMULATION_TEST tests the BROWNIAN_MOTION_SIMULATION library.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    10 June 2018
    #
    #  Author:
    #
    #    John Burkardt
    #
    print('')
    print('BROWNIAN_MOTION_SIMULATION_TEST')
    print('  Python version')
    print('  Test the BROWNIAN_MOTION_SIMULATION library.')
#
#  Compute the path of a particle undergoing Brownian motion.
#
    m = 1
    n = 1001
    d = 10.0
    t = 1.0
    x = brownian_motion_simulation(m, n, d, t)
    brownian_motion_display(m, n, x)
#
#  Estimate the average displacement of the particle from the origin
#  as a function of time.
#
    k = 40
    n = 1001
    d = 10.0
    t = 1.0

    dsq = brownian_displacement_simulation(k, n, m, d, t)
    brownian_displacement_display(k, n, m, d, t, dsq)
#
#  Terminate.
#
    print('')
    print('BROWNIAN_MOTION_SIMULATION_TEST')
    print('  Normal end of execution.')
    return


def timestamp():

    # *****************************************************************************80
    #
    # TIMESTAMP prints the date as a timestamp.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    06 April 2013
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Parameters:
    #
    #    None
    #

    t = time.time()
    print(time.ctime(t))

    return None


if __name__ == '__main__':
    obj = BrownMotion()
    obj.plot_occ()
    timestamp()

    print('')
    print('BROWNIAN_MOTION_SIMULATION_TEST')
    print('  Python version')
    print('  Test the BROWNIAN_MOTION_SIMULATION library.')

    m = 5
    n = 1001
    d = 10.0
    t = 1.0
    k = 40
    x = brownian_motion_simulation(m, n, d, t)
    dsq = brownian_displacement_simulation(k, n, m, d, t)

    brownian_motion_display(m, n, x)
    brownian_displacement_display(k, n, m, d, t, dsq)

    print('')
    print('BROWNIAN_MOTION_SIMULATION_TEST')
    print('  Normal end of execution.')

    timestamp()
