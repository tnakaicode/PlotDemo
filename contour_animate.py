import matplotlib.pyplot as plt
import matplotlib.animation
import numpy as np
import time

x = np.linspace(-3 * np.pi, 3 * np.pi)
X, Y = np.meshgrid(x, x)


def f(x, y, alpha, beta):
    return (np.sin(X + alpha) +
            np.sin(Y * (1 + np.sin(beta) * .4) + alpha))**2


alpha = np.linspace(0, 2 * np.pi, num=34)
levels = 10
cmap = plt.cm.magma


fig, axs = plt.subplots()
props = dict(boxstyle='round', facecolor='wheat')
timelabel = axs.text(0.85, 0.9, "", transform=axs.transAxes,
                     ha="right", bbox=props)
t = np.ones(10) * time.time()
p = [axs.contourf(X, Y, f(X, Y, 0, 0), levels, cmap="jet")]


def update(i):
    fig.clear()
    axs = fig.add_subplot(111)
    for tp in p[0].collections:
        tp.remove()
    p[0] = axs.contourf(X, Y, f(X, Y, alpha[i], alpha[i]), levels, cmap="jet")
    t[1:] = t[0:-1]
    t[0] = time.time()
    axs.text(0.85, 0.9, "{:.3f} fps".format(-1. / np.diff(t).mean()),
             transform=axs.transAxes, ha="right", bbox=props)
    fig.colorbar(p[0])
    return p[0].collections + [timelabel]


ani = matplotlib.animation.FuncAnimation(fig, update, frames=len(alpha),
                                         interval=10, blit=True, repeat=True)
ani.save("./img/contour_animate.gif")
plt.show()
