#!/usr/bin/env python3

# ------------------------------------------------------------------------------
__author__ = "Nicolas Richart"
__copyright__ = "Copyright (C) 2016-2018, EPFL (Ecole Polytechnique Fédérale" \
                " de Lausanne) Laboratory (LSMS - Laboratoire de Simulation" \
                " en Mécanique des Solides)"
__credits__ = ["Nicolas Richart"]
__license__ = "L-GPLv3"
__maintainer__ = "Nicolas Richart"
__email__ = "nicolas.richart@epfl.ch"
# ------------------------------------------------------------------------------

import numpy as np
import python_fe as pfe
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# sim_u = np.genfromtxt("disp.csv", delimiter=",", names=True)

L = 1.
Ne = 200  # int(sim_u['node'][-1])
F = np.zeros(Ne + 1)  # sim_u['force']
F[-1] = - 9.81
blocked = np.zeros(Ne + 1)  # sim_u['blocked']
blocked[0] = 1
u = np.zeros(Ne + 1)  # sim_u['disp']

trusses = pfe.TrussFE(Ne=Ne,
                      F={f: [np.where(F == f)] for f in np.unique(F)},
                      blocked=(np.where(blocked == 1),
                               u[blocked == 1]))

solver = pfe.DynamicSolver(trusses, delta_t=0.001)

# for s in range(200):
#     solver.solveStep()


fig, ax = plt.subplots()

x = np.arange(Ne+1) * L / Ne        # x-array
line, = ax.plot(x, trusses.u)


def animate(i):
    solver.solveStep()
    line.set_ydata(trusses.u)  # update the data
    plt.ylim(np.min(trusses.u), np.max(trusses.u))
    return line,


ani = animation.FuncAnimation(fig, animate, np.arange(1, 200),
                              interval=25, blit=True)

plt.show()
