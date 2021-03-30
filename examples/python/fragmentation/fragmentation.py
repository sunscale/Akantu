#!/usr/bin/env python
# coding: utf-8

# import akantu
import subprocess
import akantu as aka
# import numpy for vector manipulation
import numpy as np
import matplotlib.pyplot as plt


# ### Setting up the *SolidMechanicsModelCohesive*
# We need to read again the material file and the mesh
# Create the Solid Mechanics Cohesive Model in Akantu

# reading material file
aka.parseInput('material.dat')
# creating mesh
spatial_dimension = 2
mesh = aka.Mesh(spatial_dimension)
mesh.read('plate.msh')

# creates the model
model = aka.SolidMechanicsModelCohesive(mesh)
model.initFull(_analysis_method=aka._static, _is_extrinsic=True)


# Initialize the solver
# configures the static solver
solver = model.getNonLinearSolver('static')
solver.set('max_iterations', 100)
solver.set('threshold', 1e-10)
solver.set("convergence_type", aka.SolveConvergenceCriteria.residual)

# Initilize a new solver (explicit Newmark with lumped mass)
model.initNewSolver(aka._explicit_lumped_mass)
# Dynamic insertion of cohesive elements
model.updateAutomaticInsertion()


# Implement Boundary and initial conditions
# Dirichlet Boundary condition
# model.applyBC(aka.FixedValue(0., aka._x), 'left')
model.applyBC(aka.FixedValue(0., aka._y), 'bottom')

model.getExternalForce()[:] = 0


# ### Generate paraview files
# Initialization for bulk vizualisation
model.setBaseName('plate')
model.addDumpFieldVector('displacement')
model.addDumpFieldVector('velocity')
model.addDumpFieldVector('external_force')
model.addDumpField('strain')
model.addDumpField('stress')
model.addDumpField('blocked_dofs')

# Initialization of vizualisation for Cohesive model
model.setBaseNameToDumper('cohesive elements', 'cohesive')
model.addDumpFieldVectorToDumper('cohesive elements', 'displacement')
model.addDumpFieldToDumper('cohesive elements', 'damage')
model.addDumpFieldVectorToDumper('cohesive elements', 'tractions')
model.addDumpFieldVectorToDumper('cohesive elements', 'opening')


# Custom Dirichlet Boundary Condition to impose constant velocity
# Boundary functor fixing the displacement as it is

class FixedDisplacement (aka.DirichletFunctor):
    '''
        Fix the displacement at its current value
    '''

    def __init__(self, axis, vel):
        super().__init__(axis)
        self.axis = axis
        self.time = 0
        self.vel = vel

    def set_time(self, t):
        self.time = t

    def __call__(self, node, flags, disp, coord):
        # sets the blocked dofs vector to true in the desired axis
        flags[int(self.axis)] = True
        disp[int(self.axis)] = self.vel*self.time


functor_r = FixedDisplacement(aka._x, 1e-1)
model.applyBC(functor_r, 'right')
functor_l = FixedDisplacement(aka._x, -1e-1)
model.applyBC(functor_l, 'left')

# Initial condition : velocity gradient:
#  in x = 0 we have -v and in x = L we have +v

nodes = model.getMesh().getNodes()
vel_field = np.zeros(nodes.shape)
vel_field[:, 0] = (2*nodes[:, 0]-L)/L*1e-1
model.getVelocity()[:] = vel_field


# ### Run the dynamical simulation
# Initialize data arrays
# Energies :
E_pot = []
E_kin = []
E_dis = []
E_rev = []
E_con = []

# Stress :
Stress = []

dt = model.getStableTimeStep()*0.1
# choose the timestep
model.setTimeStep(dt)
# set maximum number of iteration
maxsteps = 5000
# solve
for i in range(0, maxsteps):
    time = dt*i
    functor_r.set_time(time)
    # fix displacements of the right boundary
    model.applyBC(functor_r, 'right')

    functor_l.set_time(time)
    # fix displacements of the left boundary
    model.applyBC(functor_l, 'left')

    if i % 10 == 0:
        model.dump()
        model.dump('cohesive elements')
        pass
    if i % 50 == 0:
        print('step {0}/{1}'.format(i, maxsteps))
    model.checkCohesiveStress()
    model.solveStep('explicit_lumped')

    Ep = model.getEnergy("potential")
    Ek = model.getEnergy("kinetic")
    Ed = model.getEnergy("dissipated")
    Er = model.getEnergy("reversible")
    Ec = model.getEnergy("contact")

    E_pot.append(Ep)
    E_kin.append(Ek)
    E_dis.append(Ed)
    E_rev.append(Er)
    E_con.append(Ec)

    Stress_field = model.getMaterial(0).getStress(aka._triangle_3)
    Stress_mean = np.mean(Stress_field)
    Stress.append(Stress_mean)


# Use the fragment Manager
fragment_manager = aka.FragmentManager(model)
fragment_manager.computeAllData()
Nb_elem_per_frag = fragment_manager.getNbElementsPerFragment()
Nb_frag = fragment_manager.getNbFragment()
print('Nb_frag:', Nb_frag)
# Average number of elements per fragment
Nb_elem_mean = np.mean(Nb_elem_per_frag)
print('average Nb elem / fragment:', Nb_elem_mean)
# knowing the element size we can get the average fragment size
s_mean = Nb_elem_mean*l
print('average fragment size:', s_mean)


# ## Plots
# Plot stress as a function of time

Time = [i*dt for i in range(0, maxsteps)]

Stress_MPa = [x/10**6 for x in Stress]

plt.figure()
plt.plot(Time, Stress)
plt.xlabel("Time [s]")
plt.ylabel("Stress [Pa]")
plt.show()


# Plot the energies as a function of time
plt.figure()
plt.plot(Time, E_pot, label='Potential Energy')
plt.plot(Time, E_kin, label='Kinetic Energy')
plt.plot(Time, E_dis, label='Dissipated Energy')
plt.plot(Time, E_rev, label='Reversible Energy')
plt.plot(Time, E_con, label='Contact Energy')
plt.legend()
plt.show()
