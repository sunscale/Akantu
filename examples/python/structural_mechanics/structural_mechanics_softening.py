#!/usr/bin/env python
# coding: utf-8

# # Test of Structural Mechanics
# In this test there is a beam consisting of three parts, all have the same materials.
# The left most node is fixed.
# On the right most node a force is applied in x direction.
#
# After a certain time, the material of the middle _element_ is waekened, lower Young's modulus.
# In each step the modulus is lowered by a coinstant factor.

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
beam.addConnectivityType(aka._bernoulli_beam_2)

# We need a `MeshAccessor` in order to change the size of the mesh entities.
beamAcc = aka.MeshAccessor(beam)

# Now we create the array to store the nodes and the connectivities and give them their size.
beamAcc.resizeConnectivity(3, aka._bernoulli_beam_2)
beamAcc.resizeNodes(4)

# #### Setting the Nodes
Nodes = beam.getNodes()
Nodes[0, :] = [0., 0.]
Nodes[1, :] = [1., 0.]
Nodes[2, :] = [2., 0.]
Nodes[3, :] = [3., 0.]


# #### Setting the Connections
Conn = beam.getConnectivity(aka._bernoulli_beam_2)
Conn[0, :] = [0, 1]
Conn[1, :] = [1, 2]
Conn[2, :] = [2, 3]

#### Ready
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
mat1ID = model.addMaterial(mat1, 'mat1')

mat2 = aka.StructuralMaterial()
mat2.E = 1e9
mat2.rho = 1.
mat2.I = 1.
mat2.Iz = 1.
mat2.Iy = 1.
mat2.A = 1.
mat2.GJ = 1.
mat2ID = model.addMaterial(mat2, 'mat2')

mat3 = aka.StructuralMaterial()
mat3.E = mat2.E / 100000
mat3.rho = 1.
mat3.I = 1.
mat3.Iz = 1.
mat3.Iy = 1.
mat3.A = mat2.A / 100
mat3.GJ = 1.
mat3ID = model.addMaterial(mat3, 'mat3')

# ##### Initializing the Model
model.initFull(aka.AnalysisMethod._implicit_dynamic)

# ##### Assigning the Materials
materials = model.getElementMaterial(aka._bernoulli_beam_2)

materials[0][0] = mat1ID
materials[1][0] = mat2ID
materials[2][0] = mat1ID


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
boundary[3, :] = False

# ### Solving the System
# Set up the system
deltaT = 1e-9
model.setTimeStep(deltaT)
solver = model.getNonLinearSolver()
solver.set("max_iterations", 100)
solver.set("threshold", 1e-8)
solver.set("convergence_type", aka.SolveConvergenceCriteria.solution)

# Perform N time steps.
#  At each step records the displacement of all three nodes in x direction.
N = 10000 * 60

disp0 = np.zeros(N)
disp1 = np.zeros(N)
disp2 = np.zeros(N)
disp3 = np.zeros(N)
times = np.zeros(N)
switchT = None
switchEnd = None

softDuration = 1000
SoftStart = (N // 2) - softDuration // 2
SoftEnd = SoftStart + softDuration
if(softDuration > 0):
    softFactor = (model.getMaterial('mat3').E
                  / model.getMaterial('mat2').E) ** (1.0 / softDuration)

mat2 = model.getMaterial('mat2')

for i in range(N):
    times[i] = deltaT * i

    if((SoftStart <= i <= SoftEnd) and (softDuration > 0)):
        if switchT is None:
            switchT = times[i]
        elif(i == SoftEnd):
            switchEnd = times[i]
        #
        mat2.E *= softFactor
    #

    model.solveStep()
    disp = model.getDisplacement()
    disp0[i] = disp[0, 0]
    disp1[i] = disp[1, 0]
    disp2[i] = disp[2, 0]
    disp3[i] = disp[3, 0]

disps = [disp0, disp1, disp2, disp3]
maxMin = [-1.0, 1.0]

for d in disps:
    maxMin[0] = max(np.max(d), maxMin[0])
    maxMin[1] = min(np.min(d), maxMin[1])

if has_matplotlib:
    #plt.plot(disp0, times, color='k', label = "left node (fix)")
    plt.plot(disp1, times, color='g', label = "middle, left node")
    plt.plot(disp2, times, color='g', linestyle = '--', label = "middle, right node")
    plt.plot(disp3, times, color='b', label = "right node")

    if(softDuration > 0):
        plt.plot((maxMin[1], maxMin[0]), (switchT, switchT),)
        plt.plot((maxMin[1], maxMin[0]), (switchEnd, switchEnd), )

    plt.title("Displacement in $x$ of the nodes")
    plt.ylabel("Time [S]")
    plt.xlabel("displacement [m]")
    plt.xlim((maxMin[1] * 1.3, maxMin[0] * 1.1))
    plt.legend()
    plt.show()

# If the softening is disabled, then the displacement looks wierd.
# Because the displacement first increases and then decreases.
# In this case `softDuration > 0` holds.
#
# However if the softening is enabled, it looks rather good. The left middle
# node will start to vibrate, because it is not pulled in the other direction.
