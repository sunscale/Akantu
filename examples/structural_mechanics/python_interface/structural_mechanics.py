#!/usr/bin/env python
# coding: utf-8

# # Test of Structural Mechanics
# We will now test the python interface of teh structural mechanics part.
# For that we will use the test `test/test_model/test_structural_mechanics_model/test_structural_mechanics_model_bernoulli_beam_2.cc`, which we will simply reproduce.
import akantu as aka

# Creating the Mesh

# Create a mesh for the two dimensional case
beam = aka.Mesh(2)

# read in the mesh description
beam.read("_bernoulli_beam_2.msh", aka.MeshIOType._miot_gmsh_struct)

# Creating the Model
model = aka.StructuralMechanicsModel(beam)

# Setting up the Model

# Creating and Inserting the Materials
mat1 = aka.StructuralMaterial()
mat1.E = 3e10
mat1.I = 0.0025
mat1.A = 0.01
model.addMaterial(mat1)

mat2 = aka.StructuralMaterial()
mat2.E = 3e10
mat2.I = 0.00128
mat2.A = 0.01
model.addMaterial(mat2)

# Initializing the Model
# model.initFull(analysis_method = aka.AnalysisMethod._static)
model.initFull()

# Assigning the Materials
materials = model.getElementMaterialMap(aka.ElementType._bernoulli_beam_2)

print(hex(materials.ctypes.data))

# Once we have written to the `materials` variable, everything becomes unstable.
# And the kernel will die.
materials[0][0] = 0
materials[1][0] = 1

print(materials)


# Setting Boundaries
M = 3600.
q = -6000.
L = 10.

forces = model.getExternalForce()
print(forces)

# Neumann
forces[2, 2] = -M
forces[0, 1] = q * L / 2
forces[0, 2] = q * L * L / 12
forces[1, 1] = q * L / 2
forces[1, 2] = -q * L * L / 12
print(forces)

# Dirichlets
boundary = model.getBlockedDOFs()
boundary[0, :] = True
boundary[1, :] = False
boundary[2, :] = False
boundary[2, 1] = True

print(model.getExternalForce())

model.solveStep()

disp = model.getDisplacement()
d1 = disp[1, 2]
d2 = disp[2, 2]
d3 = disp[1, 0]

print(d1, 5.6 / 4800)
