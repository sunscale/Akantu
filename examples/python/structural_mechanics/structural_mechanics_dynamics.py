#!/usr/bin/env python
# coding: utf-8

import akantu as aka
import numpy
import numpy as np
try:
    import matplotlib.pyplot as plt
    has_matplotlib = True
except ImportError:
    has_matplotlib = False

# ### Creating the Mesh
# Create a mesh for the two dimensional case
el_type = aka._bernoulli_beam_2
beam = aka.Mesh(2)

# We now create the connectivity array for the beam.
beam.addConnectivityType(el_type)

# We need a `MeshAccessor` in order to change the size of the mesh entities.
beamAcc = aka.MeshAccessor(beam)

# Now we create the array to store the nodes and the connectivities and give them their size.
nb_elem = 40
L = 2
beamAcc.resizeConnectivity(nb_elem, el_type)
beamAcc.resizeNodes(nb_elem + 1)

# #### Setting the Nodes
Nodes = beam.getNodes()
length = L / nb_elem
Nodes[:, :] = 0.
Nodes[:, 0] = np.arange(nb_elem+1) * length

# #### Setting the Connections
Conn = beam.getConnectivity(el_type)

for e in range(nb_elem):
    Conn[e, :] = [e, e + 1]

# #### Ready
# We have to make the mesh ready.
beamAcc.makeReady()

# ### Creating the Model
model = aka.StructuralMechanicsModel(beam)

if el_type == aka._bernoulli_beam_3:
    normal = beam.getDataReal("extra_normal", el_type)

    for e in range(nb_elem):
        normal[e, :] = [0, 0, 1]

# #### Setting up the Modell
# ##### Creating and Inserting the Materials
mat1 = aka.StructuralMaterial()
mat1.E = 1e9
mat1.rho = 10.
mat1.I = 1.
mat1.Iz = 1.
mat1.Iy = 1.
mat1.A = 1.
mat1.GJ = 1.
model.addMaterial(mat1, 'mat1')

# ##### Initializing the Model
model.initFull(aka.AnalysisMethod._implicit_dynamic)

# ##### Assigning the Materials
materials = model.getElementMaterial(el_type)
materials[:, :] = 0

# ##### Setting Boundaries
# Neumann
F = 1e4
no_print = int(nb_elem / 2)

#  Apply a force of `10` at the last (right most) node.
forces = model.getExternalForce()
forces[:, :] = 0
forces[no_print, 1] = F

# Dirichlets
# Block all dofs of the first node, since it is fixed.
#  All other nodes have no restrictions
boundary = model.getBlockedDOFs()
boundary[:, :] = False

boundary[0, 0] = True
boundary[0, 1] = True

if el_type == aka._bernoulli_beam_3:
    boundary[0, 2] = True

boundary[nb_elem, 1] = True

# ### Solving the System
# Set up the system
deltaT = 1e-6
model.setTimeStep(deltaT)
solver = model.getNonLinearSolver()
solver.set("max_iterations", 100)
solver.set("threshold", 1e-8)
solver.set("convergence_type", aka.SolveConvergenceCriteria.solution)

model.assembleMatrix("M")
M_ = model.getDOFManager().getMatrix("M")
M = aka.AkantuSparseMatrix(M_)

model.assembleMatrix("K")
K_ = model.getDOFManager().getMatrix("K")
K = aka.AkantuSparseMatrix(K_)

C_ = model.getDOFManager().getMatrix("C")
C_.add(M_, 0.00001)
C_.add(K_, 0.00001)

def analytical_solution(time, L, rho, E, A, I, F):
    omega = np.pi**2 / L**2 * np.sqrt(E * I / rho)
    sum = 0.
    N = 110
    for n in range(1, N, 2):
        sum += (1. - np.cos(n * n * omega * time)) / n**4

    return 2. * F * L**3 / np.pi**4 / E / I * sum

# Perform N time steps.
#  At each step records the displacement of all three nodes in x direction.
N = 900

mat1 = model.getMaterial('mat1')

disp = model.getDisplacement()
velo = model.getVelocity()
disp[:, :] = 0.

displs = np.zeros(N)

ekin = np.zeros(N)
epot = np.zeros(N)
ework = np.zeros(N)
_ework = 0.

for i in range(1, N):
    model.solveStep()
    displs[i] = disp[no_print, 1]

    _ework += F * velo[no_print, 1] * deltaT

    ekin[i] = model.getEnergy("kinetic")
    epot[i] = model.getEnergy("potential")
    ework[i] = _ework


def sol(x):
    return analytical_solution(x, L, mat1.rho, mat1.E,
                               mat1.A, mat1.I, F)


if has_matplotlib:
    times = np.arange(N) * deltaT
    plt.plot(times, sol(times))
    plt.plot(times, displs)
    plt.plot(times, displs - sol(times))

    # What I do not fully understand is why the middle node first go backwards
    # until it goes forward. I could imagine that there is some vibration,
    # because everything is in rest.
    np.max(displs - sol(times))
    plt.plot(times, ekin+epot)
    plt.plot(times, ework)
