# https://flothesof.github.io/charged-particle-trajectories-E-and-B-fields.html
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D


class Particle (object):

    def __init__(self):
        self.solver = ode(self.newton).set_integrator('dopri5')
        self.solver1 = ode(self.newton1).set_integrator('dopri5')
        self.solver2 = ode(self.newton2).set_integrator('dopri5')

    def e_of_x(self, x):
        return 10 * np.sign(np.sin(2 * np.pi * x / 25))

    def b_xyz(self, x=1, y=1, z=1):
        return 2 * x + np.sin(2 * np.pi * y / 25) + np.sin(2 * np.pi * z / 25)

    def newton(self, t, Y, q, m):
        """
        Computes the derivative of the state vector y according to the equation of motion:
        Y is the state vector (x, y, z, u, v, w) === (position, velocity).
        returns dY/dt.
        """
        x, y, z = Y[0], Y[1], Y[2]
        u, v, w = Y[3], Y[4], Y[5]

        alpha = q / m * self.b_xyz(x, y, z)
        return np.array([u, v, w, 0.5, alpha * w + self.e_of_x(x), -alpha * v])

    def newton1(self, t, Y, q, m, B):
        """
        Computes the derivative of the state vector y according to the equation of motion:
        Y is the state vector (x, y, z, u, v, w) === (position, velocity).
        returns dY/dt.
        """
        x, y, z = Y[0], Y[1], Y[2]
        u, v, w = Y[3], Y[4], Y[5]

        alpha = q / m * B
        return np.array([u, v, w, 0, alpha * w, -alpha * v])

    def newton2(self, t, Y, q, m, B, E):
        """
        Computes the derivative of the state vector y according to the equation of motion:
        Y is the state vector (x, y, z, u, v, w) === (position, velocity).
        returns dY/dt.
        """
        x, y, z = Y[0], Y[1], Y[2]
        u, v, w = Y[3], Y[4], Y[5]

        alpha = q / m
        return np.array([u, v, w, 0, alpha * B * w + E, -alpha * B * v])

    def compute_trajectory(self, m=1, q=1):
        r = ode(self.newton1).set_integrator('dopri5')
        r.set_initial_value(initial_conditions,
                            t0).set_f_params(m, q, 1.0, 10.)
        positions = []
        t1 = 200
        dt = 0.05
        while r.successful() and r.t < t1:
            r.set_f_params(m, q, 1.0, self.e_of_x(r.y[0]))
            r.integrate(r.t + dt)
            positions.append(r.y[:3])

        return np.array(positions)


if __name__ == '__main__':
    obj = Particle()

    t0 = 0
    x0 = np.array([0, 0, 0])
    v0 = np.array([1, 1, 0])
    initial_conditions = np.concatenate((x0, v0))

    obj.solver.set_initial_value(initial_conditions, t0)
    obj.solver.set_f_params(-1.0, 1.0)

    obj.solver1.set_initial_value(initial_conditions, t0)
    obj.solver1.set_f_params(1.0, 1.0, 1.0)

    obj.solver2.set_initial_value(initial_conditions, t0)
    obj.solver2.set_f_params(1.0, 1.0, 1.0, 10.0)

    pos = []
    t1 = 10
    dt = 0.05
    while obj.solver.successful() and obj.solver.t < t1:
        print(obj.solver.t, *obj.solver.y)
        obj.solver.integrate(obj.solver.t + dt)
        pos.append(obj.solver.y[:3])

    pos = np.array(pos)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot3D(pos[:, 0], pos[:, 1], pos[:, 2])
    plt.xlabel('x')
    plt.ylabel('y')
    ax.set_zlabel('z')
    plt.show()
