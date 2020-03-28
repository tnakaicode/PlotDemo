import numpy as np
import sys
from matplotlib import pyplot as plt
from IPython import display


from base import plot2d

def main():
    obj = NavierStokes()

    # グラフ描画設定
    xlabel = "$x$"
    ylabel = "$y$"
    zlabel = "$(p-p_\min)/(p_\max-p_\min)$"
    tlabelx = obj.lx / 50
    tlabely = obj.ly / 50
    zmin = obj.p.min()
    zmax = obj.p.max()
    origin = "lower"
    interp = "bicubic"
    cbar_ticks = np.linspace(zmin, zmax, 11)
    extent = (obj.x1d.min(), obj.x1d.max(), obj.y1d.min(), obj.y1d.max())
    timestamp = "%0.4f s"  # timestamp%実数 として使用
    fontsize = 15

    # グラフ描画
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    # velo = ax.quiver(
    #    obj.x[::], obj.y[::], obj.u[::], obj.v[::], units="xy", pivot="middle", color="black", scale=25
    # )
    img = ax.imshow(
        obj.p,
        extent=extent,
        origin=origin,
        interpolation=interp,
        cmap="jet",
        vmin=zmin,
        vmax=zmax,
    )
    cbar = fig.colorbar(img, ticks=cbar_ticks, label=zlabel)
    text = ax.text(tlabelx, tlabely, timestamp %
                   (0.0 * obj.dt), fontsize=15, color="black")
    fig.canvas.draw()
    display.display(fig)
    display.clear_output(wait=True)
    # plt.pause(0.001)

    for n in range(0, 3):
        # print(n)
        obj.computeAuxiallyVelocity()  # 中間速度を計算
        obj.computeDivergenceAuxiallyVelocity()         # 中間速度の発散を計算
        obj.computePressurePoisson()            # 圧力を計算
        obj.computeVelocity()       # 中間速度を修正して速度を計算

        # if (n+1)%10==0:
        text.set_text(timestamp % ((n + 1) * obj.dt))

        #velo.set_UVC(obj.u[::], obj.v[::])
        img.set_data((obj.p[::] - obj.p.min()) /
                     (obj.p.max() - obj.p.min()))  # 圧力を0~1に規格化する

        zmin = 0.0  # p.min()
        zmax = 1.0  # p.max()
        cbar.set_clim(vmin=zmin, vmax=zmax)
        cbar_ticks = np.linspace(zmin, zmax, 11)
        cbar.set_ticks(cbar_ticks)
        cbar.draw_all()

        display.display(fig)
        display.clear_output(wait=True)
        # plt.pause(0.001)
    plt.show()


if __name__ == "__main__":
    obj = NavierStokes()
    obj.run()
