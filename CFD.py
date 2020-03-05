import numpy as np
import sys
from matplotlib import pyplot as plt
from IPython import display


from base import plot2d


class NavierStokes (plot2d):

    def __init__(self):
        plot2d.__init__(self)

        # 流体の密度[kg/m^3]
        self.rho = 1e3
        # 流体の動粘度[m^2/s]
        self.ν = 1e-6
        # フタの移動速度[m/s]
        self.Uwall = 0.01

        # unit [m]
        self.lx = 0.01
        self.ly = 0.01 / 2
        self.nx = 21
        self.ny = 11
        self.x1d = np.linspace(0, self.lx, self.nx)
        self.y1d = np.linspace(0, self.ly, self.ny)
        self.x, self.y = np.meshgrid(self.x1d, self.y1d)
        self.dx = self.lx / (self.nx - 1)
        self.dy = self.ly / (self.ny - 1)

        self.lt = 3 * self.lx / self.Uwall
        self.dt = 0.1 * self.dx / self.Uwall
        self.nt = int(np.ceil(self.lt / self.dt))

        # 加速係数
        self.t = 0
        self.alp = 1.74
        self.eps = 1e-5

        # Veclocity
        self.u = np.zeros_like(self.x)
        self.v = np.zeros_like(self.y)
        self.p = np.zeros_like(self.x)
        self.u_aux = np.zeros_like(self.u)
        self.v_aux = np.zeros_like(self.v)
        self.tetha = np.zeros_like(self.p)
        self.velocityBoundary()

        # figure
        xyz = self.norm_p()
        self.img = self.axs.contourf(self.x, self.y, self.p, cmap="jet")
        self.vec = self.axs.quiver(self.x, self.y, self.u, self.v, scale=1000)
        self.cbr = self.fig.colorbar(self.img)
        self.fig.savefig(self.tmpdir + "CFD_{:d}.png".format(self.t))

    def run(self):
        for n in range(500):
            self.t += 1
            # 中間速度を計算
            self.computeAuxiallyVelocity()
            # 中間速度の発散を計算
            self.computeDivergenceAuxiallyVelocity()
            # 圧力を計算
            self.computePressurePoisson()
            # 中間速度を修正して速度を計算
            self.computeVelocity()

            # self.axs.clear()
            self.cbr.remove()
            for col in self.img.collections:
                self.axs.collections.remove(col)

            # for col in self.fig.collections:
            #    self.fig.collections.remove(col)

            self.img = self.axs.contourf(self.x, self.y, self.p, cmap="jet")
            self.cbr = self.fig.colorbar(self.img)
            self.vec.set_UVC(self.u * 1000, self.v * 1000)

            cbr_ticks = np.linspace(self.p.min(), self.p.max(), 11)
            # self.cbr.update_bruteforce(self.img)
            # self.cbr.set_ticks(cbr_ticks)
            # plt.draw()
            # self.cbr.draw_all()
            self.fig.savefig(self.tmpdir + "CFD_{:d}.png".format(self.t))
            self.fig.savefig(self.tmpdir + "CFD.png")

    def norm_p(self):
        p_min = self.p.min()
        p_max = self.p.max()
        return (self.p - p_min) / (p_max - self.p)

    def velocityBoundary(self):
        """
        速度に境界条件を反映する関数．
        """

        # Left Wall
        self.u[:, 0] = 0.0
        self.v[:, 0] = 0.0
        # Right Wall
        self.u[:, -1] = 0.0
        self.v[:, -1] = 0.0
        # Bot Wall
        self.u[0, :] = 0.0
        self.v[0, :] = 0.0
        # Top Wall
        self.u[-1, :] = self.Uwall
        self.v[-1, :] = self.Uwall / 2

    def computeAuxiallyVelocity(self):
        """
        中間速度を計算する関数．
        """
        # val = u[1:-1, 1:-1] - (
        #     dt * (
        #         u[1:-1, 1:-1] * (u[1:-1, 2:] - u[1:-1, :-2]) / (dx * 2.0) +
        #         v[1:-1, 1:-1] * (u[2:, 1:-1] - u[:-2, 1:-1]) / (dy * 2.0)
        #     ) +
        #     dt * ν * (
        #         (u[1:-1, 2:] - 2.0 * u[1:-1, 1:-1] + u[1:-1, :-2]) / (dx**2) +
        #         (u[2:, 1:-1] - 2.0 * u[1:-1, 1:-1] + u[:-2, 1:-1]) / (dy**2)
        #     )
        # )
        #
        # val = v[1:-1, 1:-1] - (
        #     dt * (
        #         u[1:-1, 1:-1] * (v[1:-1, 2:] - v[1:-1, :-2]) / (dx * 2.0) +
        #         v[1:-1, 1:-1] * (v[2:, 1:-1] - v[:-2, 1:-1]) / (dy * 2.0)
        #     ) +
        #     dt * ν * (
        #         (v[1:-1, 2:] - 2.0 * v[1:-1, 1:-1] + v[1:-1, :-2]) / (dx**2) +
        #         (v[2:, 1:-1] - 2.0 * v[1:-1, 1:-1] + v[:-2, 1:-1]) / (dy**2)
        #     )
        # )

        u_11 = self.u[1:-1, 1:-1]
        u_11_20 = self.u[1:-1, 2:]
        u_11_21 = self.u[1:-1, :-2]
        u_20_11 = self.u[2:, 1:-1]
        u_21_11 = self.u[:-2, 1:-1]

        v_11 = self.v[1:-1, 1:-1]
        v_11_20 = self.v[1:-1, 2:]
        v_11_21 = self.v[1:-1, :-2]
        v_20_11 = self.v[2:, 1:-1]
        v_21_11 = self.v[:-2, 1:-1]

        u_1120 = u_11 * (u_11_20 - u_11_21) / (2 * self.dx)
        u_1121 = v_11 * (u_20_11 - u_21_11) / (2 * self.dy)
        u_2011 = (u_11_20 - 2.0 * u_11 + u_11_21) / self.dx**2
        u_2111 = (u_20_11 - 2.0 * u_11 + u_21_11) / self.dy**2

        v_1120 = u_11 * (v_11_20 - v_11_21) / (2 * self.dx)
        v_1121 = v_11 * (v_20_11 - v_21_11) / (2 * self.dy)
        v_2011 = (v_11_20 - 2.0 * v_11 + v_11_21) / (self.dx**2)
        v_2111 = (v_20_11 - 2.0 * v_11 + v_21_11) / (self.dx**2)

        self.u_aux[1:-1, 1:-1] = u_11
        self.u_aux[1:-1, 1:-1] += -self.dt * (u_1120 + u_1121)
        self.u_aux[1:-1, 1:-1] += +self.dt * self.ν * (u_2011 + u_2111)

        self.v_aux[1:-1, 1:-1] = v_11
        self.v_aux[1:-1, 1:-1] += -self.dt * (v_1120 + v_1121)
        self.v_aux[1:-1, 1:-1] += +self.dt * self.ν * (v_2011 + v_2111)

        self.velocityBoundary()

    def computeDivergenceAuxiallyVelocity(self):
        """
        中間速度の発散を計算する関数．
        """
        u_aux_11_20 = self.u_aux[1:-1, 2:]
        u_aux_11_21 = self.u_aux[1:-1, :-2]
        u_aux_112 = (u_aux_11_20 - u_aux_11_21) / (2 * self.dx)

        v_aux_20_11 = self.v_aux[1:-1, 2:]
        v_aux_21_11 = self.v_aux[1:-1, :-2]
        v_aux_211 = (v_aux_20_11 - v_aux_21_11) / (2 * self.dy)
        self.tetha[1:-1, 1:-1] = u_aux_112 + v_aux_211

    def computePressurePoisson(self):
        """
        圧力Poisson方程式を解いて圧力を計算する関数．
        """

        dp = np.zeros_like(self.p)
        err = 1.0
        ite = 0
        while err > self.eps:
            # forを使った逐次計算
            # for j in range(1,Ny-1):
            #    for i in range(1,Nx-1):
            #        dp[j,i] =( dy**2*(p[j+1,i]+p[j-1,i]) + dx**2*(p[j,i+1]+p[j,i-1]) - (tetha[j,i]*ρ/dt)*(dx**2*dy**2)\
            #                  )/(2.*(dx**2+dy**2)) - p[j,i]
            #        p[j,i] += α*dp[j,i]

            # Black,odd-numbered row
            dp[1:-1:2, 2:-1:2] = (
                self.dy**2 * (self.p[1:-1:2, 1:-2:2] + self.p[1:-1:2, 3::2]) +
                self.dx**2 * (self.p[:-2:2, 2:-1:2] + self.p[2::2, 2:-1:2]) -
                (self.dx * self.dy)**2 * self.rho /
                self.dt * self.tetha[1:-1:2, 2:-1:2]
            ) / (2.0 * (self.dx**2 + self.dy**2)) - self.p[1:-1:2, 2:-1:2]
            self.p[1:-1:2, 2:-1:2] += self.alp * dp[1:-1:2, 2:-1:2]

            # Black,even-numbered row
            dp[2:-1:2, 1:-1:2] = (
                self.dy**2 * (self.p[2:-1:2, :-2:2] + self.p[2:-1:2, 2::2]) +
                self.dx**2 * (self.p[1:-2:2, 1:-1:2] + self.p[3::2, 1:-1:2]) -
                (self.dx * self.dy)**2 * self.rho /
                self.dt * self.tetha[2:-1:2, 1:-1:2]
            ) / (2.0 * (self.dx**2 + self.dy**2)) - self.p[2:-1:2, 1:-1:2]
            self.p[2:-1:2, 1:-1:2] += self.alp * dp[2:-1:2, 1:-1:2]

            # Red,odd-numbered row
            dp[1:-1:2, 1:-1:2] = (
                self.dy**2 * (self.p[1:-1:2, :-2:2] + self.p[1:-1:2, 2::2]) +
                self.dx**2 * (self.p[:-2:2, 1:-1:2] + self.p[2::2, 1:-1:2]) -
                (self.dx * self.dy)**2 * self.rho /
                self.dt * self.tetha[1:-1:2, 1:-1:2]
            ) / (2.0 * (self.dx**2 + self.dy**2)) - self.p[1:-1:2, 1:-1:2]
            self.p[1:-1:2, 1:-1:2] += self.alp * dp[1:-1:2, 1:-1:2]

            # Red,even-numbered row
            dp[2:-1:2, 2:-1:2] = (
                self.dy**2 * (self.p[2:-1:2, 1:-2:2] + self.p[2:-1:2, 3::2]) +
                self.dx**2 * (self.p[1:-2:2, 2:-1:2] + self.p[3::2, 2:-1:2]) -
                (self.dx * self.dy)**2 * self.rho /
                self.dt * self.tetha[2:-1:2, 2:-1:2]
            ) / (2.0 * (self.dx**2 + self.dy**2)) - self.p[2:-1:2, 2:-1:2]
            self.p[2:-1:2, 2:-1:2] += self.alp * dp[2:-1:2, 2:-1:2]

            # ノイマン境界条件
            self.p[0, :] = self.p[1, :]
            self.p[:, 0] = self.p[:, 1]
            self.p[:, -1] = self.p[:, -2]
            self.p[-1, :] = self.p[-2, :]

            err_d = np.sum(np.abs(self.p))
            if err_d < 1e-20:
                err_d = 1.0
                # 全てのpが0だと分母が0になるので，合計小さいときは1にする

            err = np.sum(np.abs(dp[:])) / err_d
            sys.stdout.write(
                "\r {:d} {:.10f} / {:.10f}".format(self.t, err, self.eps))
            sys.stdout.flush()
        print("\n")

    def computeVelocity(self):
        """
        時刻n+1の速度を計算する関数．
        """

        self.u[1:-1, 1:-1] = self.u_aux[1:-1, 1:-1] - self.dt / self.rho * \
            (self.p[1:-1, 2:] - self.p[1:-1, :-2]) / (self.dx * 2.0)
        self.v[1:-1, 1:-1] = self.v_aux[1:-1, 1:-1] - self.dt / self.rho * \
            (self.p[2:, 1:-1] - self.p[:-2, 1:-1]) / (self.dy * 2.0)
        self.velocityBoundary()


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
