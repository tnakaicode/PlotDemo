import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib import animation
from mpl_toolkits.axes_grid1 import make_axes_locatable


def plot_contour_sub(mesh, func, loc=[0, 0], txt="", title="name", dirname="./", level=None):
    """
    plot 2d data
    """
    sx, sy = loc
    nx, ny = func.shape
    xs, ys = mesh[0][0, 0], mesh[1][0, 0]
    xe, ye = mesh[0][0, -1], mesh[1][-1, 0]
    dx, dy = mesh[0][0, 1] - mesh[0][0, 0], mesh[1][1, 0] - mesh[1][0, 0]
    mx, my = int((sy - ys) / dy), int((sx - xs) / dx)
    fig, ax = plt.subplots()
    divider = make_axes_locatable(ax)
    ax.set_aspect('equal')

    ax_x = divider.append_axes("bottom", 1.0, pad=0.5, sharex=ax)
    ax_x.plot(mesh[0][mx, :], func[mx, :])
    ax_x.set_title("y = {:.2f}".format(sy))
    ax_x.xaxis.grid(True, zorder=0)
    ax_x.yaxis.grid(True, zorder=0)

    ax_y = divider.append_axes("right", 1.0, pad=0.5, sharey=ax)
    ax_y.plot(func[:, my], mesh[1][:, my])
    ax_y.set_title("x = {:.2f}".format(sx))
    ax_y.xaxis.grid(True, zorder=0)
    ax_y.yaxis.grid(True, zorder=0)

    tx, ty = 1.1, 0.0
    plt.text(tx, ty, txt, transform=ax_x.transAxes)

    im = ax.contourf(*mesh, func, cmap="jet", levels=level)
    ax.set_title(title)
    ax.xaxis.grid(True, zorder=0)
    ax.yaxis.grid(True, zorder=0)
    plt.colorbar(im, ax=ax, shrink=0.9)
    plt.tight_layout()
    plt.savefig(dirname + ".png")
    plt.close()


def plot_contour_xyz(x, y, z, loc=[0, 0], title="name", pngfile="./name"):
    fig, ax = plt.subplots()
    divider = make_axes_locatable(ax)
    ax.set_aspect('equal')
    ax.tricontour(x, y, z, levels=14, linewidths=0.5, colors='k')
    ax.plot(x, y, 'ko', ms=3)
    cntr = ax.tricontourf(x, y, z, levels=14, cmap="RdBu_r")


if __name__ == '__main__':
    np.random.seed(19680801)
    npts = 400
    x = np.random.uniform(-1, 1, npts) * 100
    y = np.random.uniform(-1, 1, npts) * 100
    z = x * np.exp(-x**2 - y**2)

    plot_contour_xyz(x, y, z)
    plt.show()
