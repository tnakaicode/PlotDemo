import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import sys
import pickle
from matplotlib import animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Lin
from OCCUtils.Construct import make_box, make_line
from OCCUtils.Construct import make_plane, make_polygon
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import dir_to_vec, vec_to_dir


class plot2d (object):

    def __init__(self):
        self.fig, self.axs = plt.subplots()
        self.axs.set_aspect('equal')
        self.axs.xaxis.grid()
        self.axs.yaxis.grid()

    def div_axs(self):
        self.div = make_axes_locatable(self.axs)
        self.axs.set_aspect('equal')

        self.ax_x = self.div.append_axes(
            "bottom", 1.0, pad=0.5, sharex=self.axs)
        self.ax_x.grid()
        self.ax_y = self.div.append_axes(
            "right", 1.0, pad=0.5, sharey=self.axs)
        self.ax_y.grid()


class plot3d (object):

    def __init__(self):
        self.fig = plt.figure()
        self.axs = self.fig.add_subplot(111, projection='3d')
        #self.axs.set_aspect('equal')

        self.axs.set_xlabel('x')
        self.axs.set_ylabel('y')
        self.axs.set_zlabel('z')

        self.axs.xaxis.grid()
        self.axs.yaxis.grid()
        self.axs.zaxis.grid()


def pnt_trf_vec(pnt=gp_Pnt(), vec=gp_Vec()):
    v = point_to_vector(pnt)
    v += vec
    return vector_to_point(v)


class plotocc (object):

    def __init__(self):
        self.display, self.start_display, self.add_menu, self.add_functionto_menu = init_display()

    def show_box(self):
        self.display.DisplayShape(gp_Pnt())
        self.display.DisplayShape(make_box(100, 100, 100))

    def show_axs_pln(self, axs=gp_Ax3(), scale=100):
        pnt = axs.Location()
        dx = axs.XDirection()
        dy = axs.YDirection()
        dz = axs.Direction()
        vx = dir_to_vec(dx).Scaled(1 * scale)
        vy = dir_to_vec(dy).Scaled(2 * scale)
        vz = dir_to_vec(dz).Scaled(3 * scale)

        pnt_x = pnt_trf_vec(pnt, vx)
        pnt_y = pnt_trf_vec(pnt, vy)
        pnt_z = pnt_trf_vec(pnt, vz)
        self.display.DisplayShape(pnt)
        self.display.DisplayShape(make_line(pnt, pnt_x), color="RED")
        self.display.DisplayShape(make_line(pnt, pnt_y), color="GREEN")
        self.display.DisplayShape(make_line(pnt, pnt_z), color="BLUE")

    def show(self):
        self.display.FitAll()
        self.start_display()


class LineDrawer(object):

    def __init__(self, dirname="./tmp/", txtname="plot_data"):
        self.trajectory = None
        self.xx = []
        self.yy = []
        self.id = 0
        self.fg = 0

        self.dirname = dirname
        self.txtname = dirname + txtname
        self.fp = open(self.txtname + ".txt", "w")

        self.init_fig()

    def run_base(self):
        self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        self.fig.canvas.mpl_connect('key_press_event', self.onkey)
        animation.FuncAnimation(
            self.fig, self.anim_animate, init_func=self.anim_init, frames=30, interval=100, blit=True)

    def init_fig(self):
        self.fig, self.axs = plt.subplots()
        self.axs.set_aspect('equal')
        self.axs.xaxis.grid()
        self.axs.yaxis.grid()
        self.divider = make_axes_locatable(self.axs)

        self.traj_line, = self.axs.plot([], [], 'o', markersize=4, mew=4)
        self.record_line, = self.axs.plot(
            [], [], 'o', markersize=4, mew=4, color='m')
        self.empty, = self.axs.plot([], [])

    def onclick(self, event):
        txt = ""

        # get mouse position and scale appropriately to convert to (x,y)
        if event.xdata is not None:
            self.trajectory = np.array([event.xdata, event.ydata])
            txt += "event.x_f {:.3f} ".format(event.xdata)
            txt += "event.y_f {:.3f} ".format(event.ydata)

        if event.button == 1:
            self.id += 1
            txt += "event {:d} ".format(event.button)
            txt += "event.x_d {:d} ".format(event.x)
            txt += "event.y_d {:d} ".format(event.y)
            txt += "flag {:d} {:d}".format(self.id % 2, self.id)
            self.xx.append(event.xdata)
            self.yy.append(event.ydata)
        elif event.button == 3:
            dat_txt = "data-{:d} num {:d}\n".format(self.fg, self.id)
            for i in range(self.id):
                dat_txt += "{:d} ".format(i)
                dat_txt += "{:.3f} ".format(self.xx[i])
                dat_txt += "{:.3f} ".format(self.yy[i])
                dat_txt += "\n"

            self.fp.write(dat_txt)
            self.fp.write("\n\n")

            self.fq = open(self.txtname + "-{:d}.txt".format(self.fg), "w")
            self.fq.write(dat_txt)
            self.fq.close()

            self.xx = []
            self.yy = []
            self.id = 0
            self.fg += 1

        print(txt)

    def onkey(self, event):
        print("onkey", event, event.xdata)
        # Record
        if event.key == 'r':
            traj = np.array([self.xx, self.yy])
            with open('traj.pickle', 'w') as f:
                pickle.dump(traj, f)
                # f.close()

    def anim_init(self):
        self.traj_line.set_data([], [])
        self.record_line.set_data([], [])
        self.empty.set_data([], [])

        return self.traj_line, self.record_line, self.empty

    def anim_animate(self, i):
        if self.trajectory is not None:
            self.traj_line.set_data(self.trajectory)

        if self.xx is not None:
            self.record_line.set_data(self.xx, self.yy)

        self.empty.set_data([], [])

        return self.traj_line, self.record_line, self.empty

    def show(self):
        try:
            plt.show()
        except AttributeError:
            pass
