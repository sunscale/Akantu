#!/usr/bin/env python3
try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    prank = comm.Get_rank()
except ImportError:
    prank = 0

import akantu as aka
import numpy as np


# -----------------------------------------------------------------------------
def solve(material_file, mesh_file, traction):
    aka.parseInput(material_file)
    spatial_dimension = 2

    # -------------------------------------------------------------------------
    # Initialization
    # -------------------------------------------------------------------------
    mesh = aka.Mesh(spatial_dimension)
    if prank == 0:
        mesh.read(mesh_file)

    mesh.distribute()

    model = aka.SolidMechanicsModel(mesh)
    model.initFull(_analysis_method=aka._static)

    model.setBaseName("plate")
    model.addDumpFieldVector("displacement")
    model.addDumpFieldVector("external_force")
    model.addDumpField("strain")
    model.addDumpField("stress")
    model.addDumpField("blocked_dofs")

    # -------------------------------------------------------------------------
    # Boundary conditions
    # -------------------------------------------------------------------------
    model.applyBC(aka.FixedValue(0.0, aka._x), "XBlocked")
    model.applyBC(aka.FixedValue(0.0, aka._y), "YBlocked")

    model.getExternalForce()[:] = 0
    trac = np.zeros(spatial_dimension)
    trac[1] = traction
    model.applyBC(aka.FromTraction(trac), "Traction")

    solver = model.getNonLinearSolver()
    solver.set("max_iterations", 2)
    solver.set("threshold", 1e-10)
    solver.set("convergence_type", aka.SolveConvergenceCriteria._residual)

    print("Solve for traction ", traction)
    model.solveStep()

    model.dump()


# -----------------------------------------------------------------------------
# main
# -----------------------------------------------------------------------------
def main():
    mesh_file = 'plate.msh'
    material_file = 'material.dat'
    traction = 1.

    solve(material_file, mesh_file, traction)


# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()
