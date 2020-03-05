import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

from base import plot2d

if __name__ == '__main__':

    n = 512
    dt = 0.01
    f = 1
    t = np.linspace(1, n, n) * dt - dt
    y = np.sin(2 * np.pi * f * t) + 0.5 * np.random.randn(t.size)

    y1 = signal.savgol_filter(y, int(n / 4) + 1, 5)
    y2 = signal.savgol_filter(y, int(n / 4) + 1, 5, mode="nearest")
    y3 = signal.savgol_filter(y, int(n / 4) + 1, 5, mode="mirror")
    y4 = signal.savgol_filter(y, int(n / 4) + 1, 5, mode="constant")
    y5 = signal.savgol_filter(y, int(n / 4) + 1, 5, mode="wrap")

    obj = plot2d()
    plt.plot(t, y)
    plt.plot(t, y1, "m", linewidth=2, label="interp")
    plt.plot(t, y2, "r", linewidth=2, label="nearest")
    plt.plot(t, y3, "c", linewidth=2, label="mirror")
    plt.plot(t, y4, "y", linewidth=2, label="constant")
    plt.plot(t, y5, "k", linewidth=2, label="wrap")
    plt.legend(loc="upper right")
    plt.xlabel("Time [s]")
    plt.ylabel("Amplitude")
    plt.savefig(obj.tmpdir + "Savitzky-Golay.png")
    plt.show()
