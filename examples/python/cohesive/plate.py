#!/usr/bin/env python3
import akantu as aka
import numpy as np


def solve(material_file, mesh_file, traction):
    aka.parseInput(material_file)
    spatial_dimension = 2

    # -------------------------------------------------------------------------
    # Initialization
    # -------------------------------------------------------------------------
    mesh = aka.Mesh(spatial_dimension)
    mesh.read(mesh_file)

    model = aka.SolidMechanicsModelCohesive(mesh)
    model.initFull(_analysis_method=aka._static,
                   _is_extrinsic=True)

    model.initNewSolver(aka._explicit_lumped_mass)

    model.setBaseName('plate')
    model.addDumpFieldVector('displacement')
    model.addDumpFieldVector('external_force')
    model.addDumpField('strain')
    model.addDumpField('stress')
    model.addDumpField('blocked_dofs')

    model.setBaseNameToDumper('cohesive elements', 'cohesive')
    model.addDumpFieldVectorToDumper('cohesive elements', 'displacement')
    model.addDumpFieldToDumper('cohesive elements', 'damage')
    model.addDumpFieldVectorToDumper('cohesive elements', 'traction')
    model.addDumpFieldVectorToDumper('cohesive elements', 'opening')

    # -------------------------------------------------------------------------
    # Boundary conditions
    # -------------------------------------------------------------------------
    model.applyBC(aka.FixedValue(0.0, aka._x), 'XBlocked')
    model.applyBC(aka.FixedValue(0.0, aka._y), 'YBlocked')

    trac = np.zeros(spatial_dimension)
    trac[int(aka._y)] = traction

    print('Solve for traction ', traction)

    model.getExternalForce()[:] = 0
    model.applyBC(aka.FromTraction(trac), 'Traction')

    solver = model.getNonLinearSolver('static')
    solver.set('max_iterations', 100)
    solver.set('threshold', 1e-10)
    solver.set("convergence_type", aka.SolveConvergenceCriteria.residual)

    model.solveStep('static')
    model.dump()
    model.dump('cohesive elements')

    model.setTimeStep(model.getStableTimeStep()*0.1)

    maxsteps = 100

    for i in range(0, maxsteps):
        print('{0}/{1}'.format(i, maxsteps))
        model.checkCohesiveStress()
        model.solveStep('explicit_lumped')
        if i % 10 == 0:
            model.dump()
            model.dump('cohesive elements')


# -----------------------------------------------------------------------------
# main
# -----------------------------------------------------------------------------
def main():
    mesh_file = 'plate.msh'
    material_file = 'material.dat'

    traction = .095
    solve(material_file, mesh_file, traction)


# -----------------------------------------------------------------------------
if __name__ == '__main__':
    main()
