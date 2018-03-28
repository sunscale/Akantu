#!/usr/bin/env python3

from __future__ import print_function
################################################################
import os
import subprocess
import numpy as np
import akantu
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


class LocalElastic:

    # declares all the internals
    def initMaterial(self, internals, params):
        self.E = params['E']
        self.nu = params['nu']
        self.rho = params['rho']
        # First Lame coefficient
        self.lame_lambda = self.nu * self.E / (
            (1 + self.nu) * (1 - 2 * self.nu))
        # Second Lame coefficient (shear modulus)
        self.lame_mu = self.E / (2 * (1 + self.nu))

    # declares all the internals
    @staticmethod
    def registerInternals():
        return ['potential']

    # declares all the internals
    @staticmethod
    def registerInternalSizes():
        return [1]

    # declares all the parameters that could be parsed
    @staticmethod
    def registerParam():
        return ['E', 'nu']

    # declares all the parameters that are needed
    def getPushWaveSpeed(self, params):
        return np.sqrt((self.lame_lambda + 2 * self.lame_mu) / self.rho)

    # compute small deformation tensor
    @staticmethod
    def computeEpsilon(grad_u):
        return 0.5 * (grad_u + np.einsum('aij->aji', grad_u))

    # constitutive law
    def computeStress(self, grad_u, sigma, internals, params):
        nquads = grad_u.shape[0]
        grad_u = grad_u.reshape((nquads, 2, 2))
        epsilon = self.computeEpsilon(grad_u)
        sigma = sigma.reshape((nquads, 2, 2))
        trace = np.trace(grad_u, axis1=1, axis2=2)
        sigma[:, :, :] = (
            np.einsum('a,ij->aij', trace,
                      self.lame_lambda * np.eye(2))
            + 2.*self.lame_mu * epsilon)

    # constitutive law tangent modulii
    def computeTangentModuli(self, grad_u, tangent, internals, params):
        n_quads = tangent.shape[0]
        tangent = tangent.reshape(n_quads, 3, 3)

        Miiii = self.lame_lambda + 2 * self.lame_mu
        Miijj = self.lame_lambda
        Mijij = self.lame_mu

        tangent[:, 0, 0] = Miiii
        tangent[:, 1, 1] = Miiii
        tangent[:, 0, 1] = Miijj
        tangent[:, 1, 0] = Miijj
        tangent[:, 2, 2] = Mijij

    # computes the energy density
    def getEnergyDensity(self, energy_type, energy_density,
                         grad_u, stress, internals, params):

        nquads = stress.shape[0]
        stress = stress.reshape(nquads, 2, 2)
        grad_u = grad_u.reshape((nquads, 2, 2))

        if energy_type != 'potential':
            raise RuntimeError('not known energy')

        epsilon = self.computeEpsilon(grad_u)

        energy_density[:, 0] = (
            0.5 * np.einsum('aij,aij->a', stress, epsilon))


################################################################
# main
################################################################

spatial_dimension = 2
akantu.parseInput('material.dat')

mesh_file = 'bar.msh'
max_steps = 250
time_step = 1e-3

# if mesh was not created the calls gmsh to generate it
if not os.path.isfile(mesh_file):
    ret = subprocess.call('gmsh -2 bar.geo bar.msh', shell=True)
    if ret != 0:
        raise Exception(
            'execution of GMSH failed: do you have it installed ?')

################################################################
# Initialization
################################################################
mesh = akantu.Mesh(spatial_dimension)
mesh.read(mesh_file)

mat = LocalElastic()
akantu.registerNewPythonMaterial(mat, "local_elastic")

model = akantu.SolidMechanicsModel(mesh)

model.initFull(_analysis_method=akantu._explicit_lumped_mass)
# model.initFull(_analysis_method=akantu._implicit_dynamic)

model.setBaseName("waves")
model.addDumpFieldVector("displacement")
model.addDumpFieldVector("acceleration")
model.addDumpFieldVector("velocity")
model.addDumpFieldVector("internal_force")
model.addDumpFieldVector("external_force")
model.addDumpField("strain")
model.addDumpField("stress")
model.addDumpField("blocked_dofs")

################################################################
# boundary conditions
################################################################

model.applyDirichletBC(FixedValue(0, akantu._x), "XBlocked")
model.applyDirichletBC(FixedValue(0, akantu._y), "YBlocked")

################################################################
# initial conditions
################################################################

displacement = model.getDisplacement()
nb_nodes = mesh.getNbNodes()
position = mesh.getNodes()

pulse_width = 1
A = 0.01
for i in range(0, nb_nodes):
    # Sinus * Gaussian
    x = position[i, 0] - 5.
    L = pulse_width
    k = 0.1 * 2 * np.pi * 3 / L
    displacement[i, 0] = A * \
        np.sin(k * x) * np.exp(-(k * x) * (k * x) / (L * L))

################################################################
# timestep value computation
################################################################
time_factor = 0.8
stable_time_step = model.getStableTimeStep() * time_factor

print("Stable Time Step = {0}".format(stable_time_step))
print("Required Time Step = {0}".format(time_step))

time_step = stable_time_step * time_factor

model.setTimeStep(time_step)

################################################################
# loop for evolution of motion dynamics
################################################################
model.assembleInternalForces()

print("step,step * time_step,epot,ekin,epot + ekin")
for step in range(0, max_steps + 1):

    model.solveStep()

    if step % 10 == 0:
        model.dump()

    epot = model.getEnergy('potential')
    ekin = model.getEnergy('kinetic')

    # output energy calculation to screen
    print("{0},{1},{2},{3},{4}".format(step, step * time_step,
                                       epot, ekin,
                                       (epot + ekin)))
