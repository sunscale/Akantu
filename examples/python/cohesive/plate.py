#!/usr/bin/env python3

from __future__ import print_function

import akantu
import numpy as np

################################################################


def solve(material_file, mesh_file, traction):
    akantu.parseInput(material_file)
    spatial_dimension = 2

    ################################################################
    # Initialization
    ################################################################
    mesh = akantu.Mesh(spatial_dimension)
    mesh.read(mesh_file)

    model = akantu.SolidMechanicsModelCohesive(mesh)
    model.initFull(akantu.SolidMechanicsModelCohesiveOptions(akantu._static,
                                                             True))

    model.initNewSolver(akantu._explicit_lumped_mass)

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

    ################################################################
    # Boundary conditions
    ################################################################

    model.applyBC(akantu.FixedValue(0.0, akantu._x), 'XBlocked')
    model.applyBC(akantu.FixedValue(0.0, akantu._y), 'YBlocked')

    trac = np.zeros(spatial_dimension)
    trac[int(akantu._y)] = traction

    print('Solve for traction ', traction)

    model.getExternalForce()[:] = 0
    model.applyBC(akantu.FromTraction(trac), 'Traction')

    solver = model.getNonLinearSolver('static')
    solver.set('max_iterations', 100)
    solver.set('threshold', 1e-10)
    solver.set('convergence_type', akantu._scc_residual)

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
        # if i < 200:
        #     model.getVelocity()[:] *= .9

    akantu.finalize()

################################################################
# main
################################################################


def main():

    import os
    mesh_file = 'plate.msh'

    if not os.path.isfile(mesh_file):
        import subprocess
        ret = subprocess.call(
            'gmsh -format msh2 -2 plate.geo {0}'.format(mesh_file),
            shell=True)
        if not ret == 0:
            raise Exception(
                'execution of GMSH failed: do you have it installed ?')

    material_file = 'material.dat'

    traction = .095
    solve(material_file, mesh_file, traction)


################################################################
if __name__ == '__main__':
    main()
