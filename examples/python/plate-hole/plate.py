#!/usr/bin/env python3

from __future__ import print_function

import akantu
import numpy as np

################################################################
# Dirichlet Boudary condition functor: fix the displacement
################################################################


class FixedValue:

    def __init__(self, value, axis):
        self.value = value
        self.axis = axis

    def operator(self, node, flags, disp, coord):
        # sets the displacement to the desired value in the desired axis
        disp[self.axis] = self.value
        # sets the blocked dofs vector to true in the desired axis
        flags[self.axis] = True

################################################################
# Neumann Boudary condition functor: impose a traction
################################################################


class FromTraction:

    def __init__(self, traction):
        self.traction = traction

    def operator(self, quad_point, force, coord, normals):
        # sets the force to the desired value in the desired axis
        force[:] = self.traction

################################################################


def solve(material_file, mesh_file, traction):
    akantu.parseInput(material_file)
    spatial_dimension = 2

    ################################################################
    # Initialization
    ################################################################
    mesh = akantu.Mesh(spatial_dimension)
    mesh.read(mesh_file)

    model = akantu.SolidMechanicsModel(mesh)
    model.initFull(akantu.SolidMechanicsModelOptions(akantu._static))

    model.setBaseName("plate")
    model.addDumpFieldVector("displacement")
    model.addDumpFieldVector("external_force")
    model.addDumpField("strain")
    model.addDumpField("stress")
    model.addDumpField("blocked_dofs")

    ################################################################
    # Boundary conditions
    ################################################################

    model.applyDirichletBC(FixedValue(0.0, akantu._x), "XBlocked")
    model.applyDirichletBC(FixedValue(0.0, akantu._y), "YBlocked")

    trac = np.zeros(spatial_dimension)
    trac[1] = traction

    print("Solve for traction ", traction)

    model.getExternalForce()[:] = 0
    model.applyNeumannBC(FromTraction(trac), "Traction")

    solver = model.getNonLinearSolver()
    solver.set("max_iterations", 2)
    solver.set("threshold", 1e-10)
    solver.set("convergence_type", akantu._scc_residual)

    model.solveStep()

    model.dump()

################################################################
# main
################################################################


def main():

    import os
    mesh_file = 'plate.msh'
    # if mesh was not created the calls gmsh to generate it
    if not os.path.isfile(mesh_file):
        import subprocess
        ret = subprocess.call(
            'gmsh -2 plate.geo {0}'.format(mesh_file), shell=True)
        if not ret == 0:
            raise Exception(
                'execution of GMSH failed: do you have it installed ?')

    material_file = 'material.dat'

    traction = 1.
    solve(material_file, mesh_file, traction)


################################################################
if __name__ == "__main__":
    main()
