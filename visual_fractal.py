from mpmath import cplot, sin, exp


def f(z):
    return 1j * z * exp(-z / sin(1j * z)) + 1


def g(z):
    return f(f(z)) + 2


cplot(g, [-1, 1], [-3, -2], points=100000)
