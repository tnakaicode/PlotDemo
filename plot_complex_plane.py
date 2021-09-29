# http://flothesof.github.io/branch-cuts-with-square-roots.html
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import hsv_to_rgb


def func_vals(f, re, im, N):
    # evaluates the complex function at the nodes of the grid
    # re and im are  tuples, re=(a, b) and im=(c, d), defining the rectangular region
    # N is the number of discrete points per unit interval

    l = re[1] - re[0]
    h = im[1] - im[0]
    resL = N * l  # horizontal resolution
    resH = N * h  # vertical resolution
    x = np.linspace(re[0], re[1], int(resL))
    y = np.linspace(im[0], im[1], int(resH))
    x, y = np.meshgrid(x, y)
    z = x + 1j * y
    return f(z)


def Hcomplex(z):
    # computes the hue corresponding to the complex number z
    H = np.angle(z) / (2 * np.pi) + 1
    return np.mod(H, 1)


def domaincol_c(w, s):
    # Classical domain coloring
    # w is the complex array of values f(z)
    # s is the constant saturation

    # detects the values w=a+ib, with a or b or both =infinity
    indi = np.where(np.isinf(w))
    indn = np.where(np.isnan(w))  # detects nans

    H = Hcomplex(w)
    S = s * np.ones_like(H)
    modul = np.absolute(w)
    V = (1.0 - 1.0 / (1 + modul**2))**0.2
    # the points mapped to infinity are colored with white; hsv_to_rgb(0, 0, 1)=(1, 1, 1)=white
    H[indi] = 0.0
    S[indi] = 0.0
    V[indi] = 1.0
    # hsv_to_rgb(0, 0, 0.5)=(0.5, 0.5, 0.5)=gray
    H[indn] = 0
    S[indn] = 0
    V[indn] = 0.5
    HSV = np.dstack((H, S, V))
    RGB = hsv_to_rgb(HSV)
    return RGB


def plot_domain(color_func, f, re=[-1, 1], im=[-1, 1], title='',
                s=0.9, N=200, daxis=None):
    w = func_vals(f, re, im, N)
    domc = color_func(w, s)
    plt.xlabel("$\Re(z)$")
    plt.ylabel("$\Im(z)$")
    plt.title(title)
    if(daxis):
        plt.imshow(domc, origin="lower", extent=[re[0], re[1], im[0], im[1]])

    else:
        plt.imshow(domc, origin="lower")
        plt.axis('off')


def numpy_sqrt(z):
    "Complex square root function."
    return np.sqrt(z)


def argand_plot(func):
    "Plots a function in the Argand (complex) plane."
    X, Y = np.meshgrid(np.linspace(-2, 2),
                       np.linspace(-2, 2))
    Z = X + 1j * Y
    ax1 = plt.gcf().add_subplot(211, projection='3d')
    ax2 = plt.gcf().add_subplot(212, projection='3d')
    ax1.plot_surface(X, Y, np.real(func(Z)), rstride=1,
                     cstride=1, cmap='viridis')
    ax1.set_xlabel('real axis')
    ax1.set_ylabel('imaginary axis')
    ax1.set_title('real part')
    ax2.plot_surface(X, Y, np.imag(func(Z)), rstride=1,
                     cstride=1, cmap='magma')
    ax2.set_xlabel('real axis')
    ax2.set_ylabel('imaginary axis')
    ax2.set_title('imaginary part')


def domain_coloring_plot(func):
    "Domain coloring of function func."
    plot_domain(domaincol_c, func, re=[-2, 2], im=[-2, 2], daxis=True)

#plt.rcParams['figure.figsize'] = (6, 7)
# plt.figure()
# argand_plot(numpy_sqrt)


plt.figure()
domain_coloring_plot(numpy_sqrt)

theta = np.linspace(0, 2 * np.pi, num=26, endpoint=False)
unit_circle = [np.exp(1j * _) for _ in theta]


def plot_along_curve(func=numpy_sqrt, param=theta, curve=unit_circle):
    "Plots curve and real/imag values of function func along given curve."
    plt.subplot(121)
    plt.plot(np.real(curve), np.imag(curve), 'o')
    x = np.real(curve)
    y = np.imag(curve)
    plt.quiver(x[:-1], y[:-1], x[1:] - x[:-1], y[1:] -
               y[:-1], scale_units='xy', angles='xy', scale=1)
    plt.xlim(-1.25, 1.25)
    plt.ylim(-1.25, 1.25)
    domain_coloring_plot(func)
    plt.subplot(122)
    plt.plot(param, np.imag(func(curve)), label='imaginary part')
    plt.plot(param, np.real(func(curve)), label='real part')
    plt.legend(loc='lower left')
    plt.xlabel('angle $\\theta$ along the circle (rad)')


plt.figure(figsize=(10, 5))
plot_along_curve()


def square_root(z, theta):
    "Square root with different branch cut defined by alpha parameter."
    argument = np.angle(z)  # between -pi and +pi
    modulus = np.abs(z)
    argument = np.mod(argument + theta, 2 * np.pi) - theta
    return np.sqrt(modulus) * np.exp(1j * argument / 2)


def normal_sqrt(z): return square_root(z, np.pi)


def real_pos_sqrt(z): return square_root(z, 2 * np.pi)


def imag_neg_sqrt(z): return square_root(z, np.pi / 2)


def angle_sqrt(z): return square_root(z, -np.pi / 4)


plt.figure()
argand_plot(normal_sqrt)

plt.figure()
domain_coloring_plot(normal_sqrt)


plt.figure()
argand_plot(real_pos_sqrt)

plt.figure()
domain_coloring_plot(real_pos_sqrt)


plt.figure()
argand_plot(imag_neg_sqrt)

plt.figure()
domain_coloring_plot(imag_neg_sqrt)

plt.figure(figsize=(10, 5))
plot_along_curve(func=imag_neg_sqrt)


plt.figure()
argand_plot(angle_sqrt)

plt.figure()
domain_coloring_plot(angle_sqrt)
plt.figure(figsize=(10, 5))
plot_along_curve(func=angle_sqrt)
plt.show()
