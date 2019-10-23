import numpy as np
from matplotlib import pyplot as plt
from IPython import display

from base import plot2d


class NavierStokes (plot2d):

    def __init__(self):
        # 流体の密度[kg/m^3]
        self.rho = 1e3
        # 流体の動粘度[m^2/s]
        self.ν = 1e-6
        # フタの移動速度[m/s]
        self.Uwall = 0.01

        # unit [m]
        self.lx = 0.01
        self.ly = 0.01
        self.nx = 21
        self.ny = 21
        x1d = np.linspace(0, self.lx, self.nx)
        y1d = np.linspace(0, self.ly, self.ny)
        self.x, self.y = np.meshgrid(x1d, y1d)
        self.dx = self.lx / (self.nx - 1)
        self.dy = self.ly / (self.ny - 1)

        self.lt = 3 * self.lx / self.Uwall
        self.dt = 0.1 * self.dx / self.Uwall
        self.nt = int(np.ceil(self.lt / self.dt))

        # 加速係数
        alh = 1.74
        eps = 1e-6

        # Veclocity
        self.u = np.zeros_like(self.x)
        self.v = np.zeros_like(self.y)

        self.p = np.empty_like(self.x)
        self.u_aux = np.empty_like(self.u)
        self.v_aux = np.empty_like(self.v)
        self.tetha = np.zeros_like(self.p)

        self.u[:, :] = 0.0
        self.v[:, :] = 0.0
        self.p[:, :] = 0
        self.velocityBoundary()

    def velocityBoundary(self):
        """
        速度に境界条件を反映する関数．
        入出力：x方向速度（numpy.array）
        入出力：y方向速度（numpy.array）
        入力：フタの移動速度(float型)
        戻り値：なし
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
        self.v[-1, :] = 0.0





    def computeAuxiallyVelocity(self):
        """
        中間速度を計算する関数．
        入出力：x方向中間速度（numpy.array）
        入出力：y方向中間速度（numpy.array）
        入力：x方向速度（numpy.array）
        入力：y方向速度（numpy.array）
        入力：離散点のx方向間隔(float型)
        入力：離散点のy方向間隔(float型)
        入力：計算時間間隔(float型)
        入力：動粘度(float型)
        入力：フタの移動速度(float型)
        戻り値：なし
        """
        u_aux[1:-1, 1:-1] += -dt * (u[1:-1, 1:-1] * (u[1:-1, 2:] - u[1:-1, :-2]) / (dx * 2.0) + v[1:-1, 1:-1] * (u[2:, 1:-1] - u[:-2, 1:-1]) / (dy * 2.0)) 
        u_aux[1:-1, 1:-1] += +dt * ν * ((u[1:-1, 2:] - 2.0 * u[1:-1, 1:-1] + u[1:-1, :-2]) / (dx**2) + (u[2:, 1:-1] - 2.0 * u[1:-1, 1:-1] + u[:-2, 1:-1]) / (dy**2))
        
        v_aux[1:-1, 1:-1] += -dt * (u[1:-1, 1:-1] * (v[1:-1, 2:] - v[1:-1, :-2]) / (dx * 2.0) + v[1:-1, 1:-1] * (v[2:, 1:-1] - v[:-2, 1:-1]) / (dy * 2.0)) 
        v_aux[1:-1, 1:-1] += +dt * ν * ((v[1:-1, 2:] - 2.0 * v[1:-1, 1:-1] + v[1:-1, :-2]) / (dx**2) + (v[2:, 1:-1] - 2.0 * v[1:-1, 1:-1] + v[:-2, 1:-1]) / (dy**2))
                             
        self.velocityBoundary()


def computeDivergenceAuxiallyVelocity(θ, u_aux, v_aux, dx, dy):
    """
    中間速度の発散を計算する関数．
    入出力：中間速度の発散（numpy.array）
    入力：x方向中間速度（numpy.array）
    入力：y方向中間速度（numpy.array）
    入力：離散点のx方向間隔(float型)
    入力：離散点のy方向間隔(float型)
    戻り値：なし
    """
    θ[1:-1, 1:-1] = (u_aux[1:-1, 2:] - u_aux[1:-1, :-2]) / (dx * 2.0)\
        + (v_aux[2:, 1:-1] - v_aux[:-2, 1:-1]) / (dy * 2.0)


def computePressurePoisson(p, θ, α, Err_tol, dx, dy, ρ, dt):
    """
    圧力Poisson方程式を解いて圧力を計算する関数．
    入出力：圧力（numpy.array）
    入力：中間速度の発散（numpy.array）
    入力：加速係数(float型)
    入力：許容誤差(float型)
    入力：離散点のx方向間隔(float型)
    入力：離散点のy方向間隔(float型)
    入力：密度(float型)
    入力：計算時間間隔(float型)
    戻り値：なし
    """

    dp = np.zeros_like(p)
    ε = 1.0
    ite = 0
    while ε > Err_tol:
        # forを使った逐次計算
        # for j in range(1,Ny-1):
        #    for i in range(1,Nx-1):
        #        dp[j,i] =( dy**2*(p[j+1,i]+p[j-1,i])\
        #                  +dx**2*(p[j,i+1]+p[j,i-1])\
        #                  -(θ[j,i]*ρ/dt)*(dx**2*dy**2)\
        #                  )/(2.*(dx**2+dy**2))\
        #                  - p[j,i]
        #        p[j,i] += α*dp[j,i]

        # Black,odd-numbered row
        dp[1:-1:2, 2:-1:2] = (dy**2 * (p[1:-1:2, 1:-2:2] + p[1:-1:2, 3::2])
                              + dx**2 * (p[:-2:2, 2:-1:2] + p[2::2, 2:-1:2])
                              - dx**2 * dy**2 * ρ / dt * θ[1:-1:2, 2:-1:2]
                              ) / (2.0 * (dx**2 + dy**2)) - p[1:-1:2, 2:-1:2]
        p[1:-1:2, 2:-1:2] += α * dp[1:-1:2, 2:-1:2]
        # Black,even-numbered row
        dp[2:-1:2, 1:-1:2] = (dy**2 * (p[2:-1:2, :-2:2] + p[2:-1:2, 2::2])
                              + dx**2 * (p[1:-2:2, 1:-1:2] + p[3::2, 1:-1:2])
                              - dx**2 * dy**2 * ρ / dt * θ[2:-1:2, 1:-1:2]
                              ) / (2.0 * (dx**2 + dy**2)) - p[2:-1:2, 1:-1:2]
        p[2:-1:2, 1:-1:2] += α * dp[2:-1:2, 1:-1:2]
        # Red,odd-numbered row
        dp[1:-1:2, 1:-1:2] = (dy**2 * (p[1:-1:2, :-2:2] + p[1:-1:2, 2::2])
                              + dx**2 * (p[:-2:2, 1:-1:2] + p[2::2, 1:-1:2])
                              - dx**2 * dy**2 * ρ / dt * θ[1:-1:2, 1:-1:2]
                              ) / (2.0 * (dx**2 + dy**2)) - p[1:-1:2, 1:-1:2]
        p[1:-1:2, 1:-1:2] += α * dp[1:-1:2, 1:-1:2]
        # Red,even-numbered row
        dp[2:-1:2, 2:-1:2] = (dy**2 * (p[2:-1:2, 1:-2:2] + p[2:-1:2, 3::2])
                              + dx**2 * (p[1:-2:2, 2:-1:2] + p[3::2, 2:-1:2])
                              - dx**2 * dy**2 * ρ / dt * θ[2:-1:2, 2:-1:2]
                              ) / (2.0 * (dx**2 + dy**2)) - p[2:-1:2, 2:-1:2]
        p[2:-1:2, 2:-1:2] += α * dp[2:-1:2, 2:-1:2]

        p[0, :] = p[1, :]   # ノイマン境界条件
        p[:, 0] = p[:, 1]   # ノイマン境界条件
        p[:, -1] = p[:, -2]  # ノイマン境界条件
        p[-1, :] = p[-2, :]  # ノイマン境界条件

        ε_d = np.sum(np.abs(p[:]))
        if ε_d < 1e-20:
            ε_d = 1.0  # 全てのpが0だと分母が0になるので，合計小さいときは1にする
        ε = np.sum(np.abs(dp[:])) / ε_d
    # print(ε)


def computeVelocity(u, v, u_aux, v_aux, p, dx, dy, ρ, dt, Uwall):
    """
    時刻n+1の速度を計算する関数．
    入出力：x方向速度（numpy.array）
    入出力：y方向速度（numpy.array）
    入力：x方向中間速度（numpy.array）
    入力：y方向中間速度（numpy.array）
    入力：圧力（numpy.array）
    入力：離散点のx方向間隔(float型)
    入力：離散点のy方向間隔(float型)
    入力：計算時間間隔(float型)
    入力：動粘度(float型)
    入力：フタの移動速度(float型)
    戻り値：なし
    """

    u[1:-1, 1:-1] = u_aux[1:-1, 1:-1] - dt / ρ * \
        (p[1:-1, 2:] - p[1:-1, :-2]) / (dx * 2.0)
    v[1:-1, 1:-1] = v_aux[1:-1, 1:-1] - dt / ρ * \
        (p[2:, 1:-1] - p[:-2, 1:-1]) / (dy * 2.0)
    velocityBoundary(u, v, Uwall)


def main():

    ρ = 1e3       # 流体の密度[kg/m^3]
    ν = 1e-6      # 流体の動粘度[m^2/s]
    Uwall = 0.01  # フタの移動速度[m/s]

    Lx = 0.01         # x方向計算領域[m]
    Ly = Lx           # x方向計算領域[m]
    Nx = 21           # x方向の離散点の個数
    Ny = Nx           # y方向の離散点の個数
    dx = Lx / (Nx - 1)  # 離散点のx方向間隔
    dy = Ly / (Ny - 1)  # 離散点のy方向間隔

    Lt = 3. * Lx / Uwall          # 計算終了時間[s]
    dt = 0.1 * dx / Uwall         # 計算時間間隔
    Nt = int(np.ceil(Lt / dt))  # 計算回数

    α = 1.74        # 加速係数
    Err_tol = 1e-6  # 許容誤差

    x1d = np.linspace(0, Lx, Nx)  # x座標
    y1d = np.linspace(0, Ly, Ny)  # y座標
    x, y = np.meshgrid(x1d, y1d)  # 2次元グリッド

    u = np.zeros((Ny, Nx))    # 速度のx方向成分
    v = np.empty_like(u)      # 速度のy方向成分
    p = np.empty_like(u)      # 圧力
    u_aux = np.empty_like(u)  # 中間速度のx方向成分
    v_aux = np.empty_like(v)  # 中間速度のy方向成分
    θ = np.zeros_like(p)      # 中間速度の発散

    # 初期化（0ステップめ）
    initialize(u, v, p, Uwall)

    # グラフ描画設定
    xlabel = "$x$"
    ylabel = "$y$"
    zlabel = "$(p-p_\min)/(p_\max-p_\min)$"
    tlabelx = Lx / 50
    tlabely = Ly / 50
    zmin = p.min()
    zmax = p.max()
    origin = "lower"
    interp = "bicubic"
    cbar_ticks = np.linspace(zmin, zmax, 11)
    extent = (x1d.min(), x1d.max(), y1d.min(), y1d.max())
    timestamp = "%0.4f s"  # timestamp%実数 として使用
    fontsize = 15

    # グラフ描画
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    velo = ax.quiver(
        x[::], y[::], u[::], v[::], units="xy", pivot="middle", color="black", scale=25
    )
    img = ax.imshow(
        p,
        extent=extent,
        origin=origin,
        interpolation=interp,
        cmap="jet",
        vmin=zmin,
        vmax=zmax,
    )
    cbar = fig.colorbar(img, ticks=cbar_ticks, label=zlabel)
    text = ax.text(tlabelx, tlabely, timestamp %
                   (0.0 * dt), fontsize=15, color="black")
    fig.canvas.draw()
    display.display(fig)
    display.clear_output(wait=True)
    # plt.pause(0.001)

    for n in range(0, Nt):
        # print(n)
        computeAuxiallyVelocity(u_aux, v_aux, u, v, dx,
                                dy, dt, ν, Uwall)  # 中間速度を計算
        computeDivergenceAuxiallyVelocity(
            θ, u_aux, v_aux, dx, dy)         # 中間速度の発散を計算
        computePressurePoisson(p, θ, α, Err_tol, dx, dy,
                               ρ, dt)            # 圧力を計算
        computeVelocity(u, v, u_aux, v_aux, p, dx, dy, ρ,
                        dt, Uwall)       # 中間速度を修正して速度を計算

        # if (n+1)%10==0:
        text.set_text(timestamp % ((n + 1) * dt))

        velo.set_UVC(u[::], v[::])
        img.set_data((p[::] - p.min()) / (p.max() - p.min()))  # 圧力を0~1に規格化する

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
    main()
