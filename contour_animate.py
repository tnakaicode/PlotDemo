import matplotlib.pyplot as plt
import matplotlib.animation
import numpy as np
import time

x = np.linspace(0, 3 * np.pi)
X, Y = np.meshgrid(x, x)


def f(x, y, alpha, beta):
    return (np.sin(X + alpha) +
            np.sin(Y * (1 + np.sin(beta) * .4) + alpha))**2


alpha = np.linspace(0, 2 * np.pi, num=34)
levels = 10
cmap = plt.cm.magma


fig, ax = plt.subplots()
props = dict(boxstyle='round', facecolor='wheat')
timelabel = ax.text(0.9, 0.9, "", transform=ax.transAxes,
                    ha="right", bbox=props)
t = np.ones(10) * time.time()
p = [ax.contourf(X, Y, f(X, Y, 0, 0), levels, cmap="jet")]


def update(i):
    for tp in p[0].collections:
        tp.remove()
    p[0] = ax.contourf(X, Y, f(X, Y, alpha[i], alpha[i]), levels, cmap="jet")
    t[1:] = t[0:-1]
    t[0] = time.time()
    timelabel.set_text("{:.3f} fps".format(-1. / np.diff(t).mean()))
    # fig.colorbar(p[0])
    return p[0].collections + [timelabel]


ani = matplotlib.animation.FuncAnimation(fig, update, frames=len(alpha),
                                         interval=10, blit=True, repeat=True)
ani.save("./contour_animate.gif")
plt.show()
