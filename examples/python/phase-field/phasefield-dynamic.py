#!/usr/bin/env python
# coding: utf-8

import py11_akantu as aka
import subprocess

geometry_file = """
h1 = 1e-4;
h2 = 1e-3;
L = 32e-3;
H = 16e-3;
l = 4e-3;
Point(1) = {0, 0, 0, h1};
Point(2) = {L, 0, 0, h1};
Point(3) = {L, H/2, 0, h2};
Point(4) = {0, H/2, 0, h2};
Point(5) = {l, 0, 0, h1};

Point(6) =  {0, 0, 0, h1};
Point(7) =  {L, -H/2, 0, h2};
Point(8) =  {0, -H/2, 0, h2};


Line(1) = {1, 5};
Line(2) = {4, 1};
Line(3) = {3, 4};
Line(4) = {2, 3};
Line(5) = {5, 2};

Line Loop(1) = {2, 3, 4, 5, 1};
Plane Surface(1) = {1};

Line(6) =  {5, 6};
Line(7) =  {6, 8};
Line(8) =  {8, 7};
Line(9) =  {7, 2};
Line Loop(2) = {6, 7, 8, 9, -5};
Plane Surface(2) = {2};


Physical Surface(8) = {1,2};
Physical Line("left") = {2,7};
Physical Line("bottom") = {8};
Physical Line("top") = {3};
Physical Line("right") = {4,9};

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
    name = virtual
    rho = 1180.    # density
    E   = 3.09e9    # young's modulus
    nu  = 0.35  # poisson's ratio
    eta = 0.0
    finite_deformation = false
]

phasefield exponential [
  name = virtual
  E = 3.09e9
  nu = 0.35
  gc = 300.
  l0 = 0.1e-3
]

"""

with open('material.dat', 'w') as f:
    f.write(material_file)


# reading material file
aka.parseInput('material.dat')
# creating mesh
spatial_dimension = 2
mesh = aka.Mesh(spatial_dimension)
mesh.read('plate.msh')


model = aka.CouplerSolidPhaseField(mesh)

solid = model.getSolidMechanicsModel()
phase = model.getPhaseFieldModel()

# initializing the Solid Mechanics Model with implicit solver for static resolution
solid.initFull(_analysis_method=aka._static)
solver = solid.getNonLinearSolver('static')
solver.set('max_iterations', 100)
solver.set('threshold', 1e-10)
solver.set("convergence_type", aka.SolveConvergenceCriteria.residual)

# adding another solver dynamic/quasi-static resolution (explicit Newmark with lumped mass) 
solid.initNewSolver(aka._explicit_lumped_mass)

# initializing the PhaseField Model with linear implicit solver for static resolution
phase.initFull(_analysis_method=aka._static)

# initializing the PhaseField Model with Newton Raphson implicit solver for static resolution
phase.getNewSolver("nonlinear_static", aka.TimeStepSolverType.static,
                   aka.NonLinearSolverType.newton_raphson)
phase.setIntegrationScheme("nonlinear_static", "damage",
                           aka.IntegrationSchemeType.pseudo_time)

solver = phase.getNonLinearSolver('nonlinear_static')
solver.set('max_iterations', 100)
solver.set('threshold', 1e-3)
solver.set("convergence_type", aka.SolveConvergenceCriteria.solution)


# Initialization for bulk vizualisation
solid.setBaseName('plate')
solid.addDumpFieldVector('displacement')
solid.addDumpFieldVector('external_force')
solid.addDumpFieldVector('velocity')
solid.addDumpField('strain')
solid.addDumpField('stress')
solid.addDumpField('damage')
solid.addDumpField('blocked_dofs')


class FixedDamage (aka.DirichletFunctor):
    '''
        Fix the damage to 0 
    '''
    def __init__(self, axis):
        super().__init__(axis)
        self.axis = axis

    def __call__(self, node, flags, dam, coord):
        # sets the blocked dofs vector to true in the desired axis
        flags[int(self.axis)] = True
        dam[int(self.axis)] = 0.0


# Dirichlet
solid.applyBC(aka.FixedValue(0., aka._x), 'top')
solid.applyBC(aka.FixedValue(0., aka._x), 'bottom')

solid.applyBC(aka.FixedValue(0., aka._x), 'left')
solid.applyBC(aka.FixedValue(0., aka._x), 'right')


solid.applyBC(aka.FixedValue(0.06e-3, aka._y), 'top')
solid.applyBC(aka.FixedValue(-0.06e-3, aka._y), 'bottom')


solid.solveStep('static')
solid.dump()


# #### **Damped dynamics resolution**
solid.setTimeStep(solid.getStableTimeStep()*0.8)

# set maximum number of iteration
maxsteps = 1000
# solve using staggered scheme
for i in range(0, maxsteps):
    if i % 100 == 0:
        print('step {0}/{1}'.format(i, maxsteps))
    model.solve('explicit_lumped', '')
    if i % 100 == 0:
        model.dump()
