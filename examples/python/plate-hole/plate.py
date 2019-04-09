#!/usr/bin/env python3

from __future__ import print_function

import akantu
import numpy as np

################################################################
# Dirichlet Boudary condition functor: fix the displacement
################################################################


class FixedValue(akantu.DirichletFunctor):

    def __init__(self, value, axis):
        super().__init__(axis)
        self.value = value
        self.axis = int(axis)

    def __call__(self, node, flags, disp, coord):
        # sets the displacement to the desired value in the desired axis
        disp[self.axis] = self.value
        # sets the blocked dofs vector to true in the desired axis
        flags[self.axis] = True

################################################################
# Neumann Boudary condition functor: impose a traction
################################################################


class FromTraction(akantu.NeumannFunctor):

    def __init__(self, traction):
        super().__init__()
        self.traction = traction

    def __call__(self, quad_point, force, coord, normals):
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

    model.applyBC(FixedValue(0.0, akantu._x), "XBlocked")
    model.applyBC(FixedValue(0.0, akantu._y), "YBlocked")

    trac = np.zeros(spatial_dimension)
    trac[1] = traction

    print("Solve for traction ", traction)

    model.getExternalForce()[:] = 0
    model.applyBC(FromTraction(trac), "Traction")

    solver = model.getNonLinearSolver()
    solver.set("max_iterations", int(2))
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
            'gmsh -format msh2 -2 plate.geo {0}'.format(mesh_file), shell=True)
        if not ret == 0:
            raise Exception(
                'execution of GMSH failed: do you have it installed ?')

    material_file = 'material.dat'

    traction = 1.
    solve(material_file, mesh_file, traction)


################################################################
if __name__ == "__main__":
    main()
