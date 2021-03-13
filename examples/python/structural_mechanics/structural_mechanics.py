#!/usr/bin/env python
# coding: utf-8

# # Test of Structural Mechanics
# In this example a beam, consisting of two elements, three nodes, is created.
# The left most node is fixed and a force is applied at the right most node.
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
beam = aka.Mesh(2)

# We now create the connectivity array for the beam.
beam.addConnectivity(aka._bernoulli_beam_2)

# We need a `MeshAccessor` in order to change the size of the mesh entities.
beamAcc = aka.MeshAccessor(beam)

# Now we create the array to store the nodes and the connectivities and give them their size.
beamAcc.resizeConnectivity(2, aka._bernoulli_beam_2)
beamAcc.resizeNodes(3)

Nodes = beam.getNodes()
Nodes[0, :] = [0., 0.]
Nodes[1, :] = [1., 0.]
Nodes[2, :] = [2., 0.]

# #### Setting the Connections
Conn = beam.getConnectivity(aka._bernoulli_beam_2)
Conn[0, :] = [0, 1]
Conn[1, :] = [1, 2]

# #### Ready
# We have to make the mesh ready.
beamAcc.makeReady()


# ### Creating the Model
model = aka.StructuralMechanicsModel(beam)

# #### Setting up the Modell
# ##### Creating and Inserting the Materials
mat1 = aka.StructuralMaterial()
mat1.E = 1e9
mat1.rho = 1.
mat1.I = 1.
mat1.Iz = 1.
mat1.Iy = 1.
mat1.A = 1.
mat1.GJ = 1.
model.addMaterial(mat1)

mat2 = aka.StructuralMaterial()
mat2.E = 1e9
mat2.rho = 1.
mat2.I = 1.
mat2.Iz = 1.
mat2.Iy = 1.
mat2.A = 1.
mat2.GJ = 1.
model.addMaterial(mat2)

# ##### Initializing the Model
model.initFull(aka._implicit_dynamic)

# ##### Assigning the Materials
materials = model.getElementMaterial(aka._bernoulli_beam_2)
materials[0][0] = 0
materials[1][0] = 1

# ##### Setting Boundaries

# Neumann
#  Apply a force of `10` at the last (right most) node.
forces = model.getExternalForce()
forces[:] = 0
forces[2, 0] = 100.

# Dirichlets
# Block all dofs of the first node, since it is fixed.
#  All other nodes have no restrictions
boundary = model.getBlockedDOFs()
boundary[0, :] = True
boundary[1, :] = False
boundary[2, :] = False

# ### Solving the System

# Set up the system
deltaT = 1e-10
model.setTimeStep(deltaT)
solver = model.getNonLinearSolver()
solver.set("max_iterations", 100)
solver.set("threshold", 1e-8)
solver.set("convergence_type", aka.SolveConvergenceCriteria.solution)

# Perform N time steps.
#  At each step records the displacement of all three nodes in x direction.
N = 1000000

disp1 = np.zeros(N)
disp2 = np.zeros(N)
disp0 = np.zeros(N)
times = np.zeros(N)

for i in range(N):
    model.solveStep()
    disp = model.getDisplacement()
    disp0[i] = disp[0, 0]
    disp1[i] = disp[1, 0]
    disp2[i] = disp[2, 0]
    times[i] = deltaT * i

disps = [disp0, disp1, disp2]
maxMin = [-1.0, 1.0]

for d in disps:
    maxMin[0] = max(np.max(d), maxMin[0])
    maxMin[1] = min(np.min(d), maxMin[1])

if has_matplotlib:
    plt.plot(disp1, times, color='g', label = "middle node")
    plt.plot(disp2, times, color='b', label = "right node")

    plt.title("Displacement in $x$ of the nodes")
    plt.ylabel("Time [S]")
    plt.xlabel("displacement [m]")

    plt.xlim((maxMin[1] * 1.3, maxMin[0] * 1.1))

    plt.legend()

    plt.show()
