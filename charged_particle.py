import numpy as np
import matplotlib.pyplot as plt

from base import plot3d, plotocc
from ChargedParticle import Particle

if __name__ == '__main__':
    obj = Particle()
    obj.base_solver()

    pos = []
    t1 = 10
    dt = 0.5
    while obj.solver.successful() and obj.solver.t < t1:
        print(obj.solver.t, *obj.solver.y)
        obj.solver.integrate(obj.solver.t + dt)
        pos.append(obj.solver.y[:3])

    pos = np.array(pos)
