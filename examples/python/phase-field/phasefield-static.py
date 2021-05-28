#!/usr/bin/env python
# coding: utf-8

import numpy as np
import py11_akantu as aka


import subprocess

geometry_file = """
element_size = 0.1;
fine_element_size = element_size;

Point(1) = {0.5, 0.5, 0, element_size};
Point(2) = {-0.5, 0.5, 0, element_size};
Point(3) = {-0.5, -0.5, 0, element_size};
Point(4) = {0.5, -0.5, 0, element_size};
Point(5) = {-0.5, 0.001, 0, element_size};
Point(6) = {0., 0.0, 0, fine_element_size};
Point(7) = {0.5, 0.0, 0, fine_element_size};
Point(8) = {-0.5, -0.001, 0, element_size};

Line(1) = {3, 4};
Line(2) = {4, 7};
Line(3) = {7, 1};
Line(4) = {1, 2};
Line(5) = {2, 5};
Line(6) = {5, 6};
Line(7) = {6, 8};
Line(8) = {8, 3};

Line Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};

Plane Surface(1) = {1};

Physical Surface("plate") = {1};

Physical Line("bottom") = {1};
Physical Line("right") = {2, 3};
Physical Line("top") = {4};
Physical Line("left") = {5,8};

"""

with open('plate.geo', 'w') as f:
    f.write(geometry_file)

ret = subprocess.run("gmsh -2 -order 1 -o plate.msh plate.geo", shell=True)
if ret.returncode:
    print("Beware, gmsh could not run: mesh is not regenerated")
else:
    print("Mesh generated")

material_file = """
material phasefield [
name = plate
	 rho = 1.
	 E = 210.0
	 nu = 0.3
         eta = 0.0 
	 Plane_Stress = false
]

phasefield exponential [
      name = plate
      l0 = 0.0075
      gc = 2.7e-3
      E  = 210.0
      nu = 0.3
]
"""

with open('material.dat', 'w') as f:
    f.write(material_file)

aka.parseInput("material.dat")

dim = 2
mesh = aka.Mesh(dim)
mesh.read("plate.msh")

model = aka.CouplerSolidPhaseField(mesh)

solid = model.getSolidMechanicsModel()
phase = model.getPhaseFieldModel()

solid.initFull(_analysis_method=aka._static)
solver = solid.getNonLinearSolver('static')
solver.set('max_iterations', 100)
solver.set('threshold', 1e-8)
solver.set("convergence_type", aka.SolveConvergenceCriteria.solution)


solid.getNewSolver("linear_static", aka.TimeStepSolverType.static,
                   aka.NonLinearSolverType.linear)
solid.setIntegrationScheme("linear_static", "displacement",
                           aka.IntegrationSchemeType.pseudo_time)


phase.initFull(_analysis_method=aka._static)
phase.getNewSolver("nonlinear_static", aka.TimeStepSolverType.static,
                   aka.NonLinearSolverType.newton_raphson)
phase.setIntegrationScheme("nonlinear_static", "damage",
                           aka.IntegrationSchemeType.pseudo_time)
solver = phase.getNonLinearSolver('nonlinear_static')
solver.set('max_iterations', 100)
solver.set('threshold', 1e-4)
solver.set("convergence_type", aka.SolveConvergenceCriteria.solution)


solid.applyBC(aka.FixedValue(0, aka._y), "bottom")
solid.applyBC(aka.FixedValue(0, aka._x), "left")

# Initialization for bulk vizualisation
solid.setBaseName('phasefield-static')
solid.addDumpFieldVector('displacement')
solid.addDumpFieldVector('external_force')
solid.addDumpField('strain')
solid.addDumpField('stress')
solid.addDumpField('damage')
solid.addDumpField('blocked_dofs')


nb_dofs = solid.getMesh().getNbNodes() * dim

increment = solid.getIncrement()
displacement = solid.getDisplacement()
displacement = displacement.reshape(nb_dofs)

blocked_dofs = solid.getBlockedDOFs()
blocked_dofs = blocked_dofs.reshape(nb_dofs)

damage = phase.getDamage()

tolerance = 1e-8

steps = 1500
increment = 1e-5

for n in range(steps):
    print("Computing iteration " + str(n + 1) + "/" + str(steps))

    solid.applyBC(aka.IncrementValue(increment, aka._y), 'top')

    mask = blocked_dofs == False

    iiter = 0
    error_disp = 1
    error_dam = 1

    displacement_prev = displacement[mask].copy()

    damage_prev = damage.copy()
    damage_prev = damage_prev

    # solve using staggered scheme
    while (error_disp > tolerance or error_dam > tolerance):
        model.solve("linear_static", "")

        displacement_new = displacement[mask]
        damage_new = damage

        delta_disp = displacement_new - displacement_prev
        delta_dam = damage_new - damage_prev

        error_disp = np.linalg.norm(delta_disp)
        error_dam = np.linalg.norm(delta_dam)

        iiter += 1

        displacement_prev = displacement_new.copy()
        damage_prev = damage_new.copy()

        print(error_dam, error_disp)
        if iiter > 500:
            raise Exception('Convergence not reached')

    if n % 50 == 0:
        solid.dump()

solid.dump()
