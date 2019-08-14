#!/usr/bin/env python3
import numpy as np
import akantu as aka


# ------------------------------------------------------------------------------
class LocalElastic(aka.Material):

    def __init__(self, model, _id):
        super().__init__(model, _id)
        super().registerParamReal('E',
                                  aka._pat_readable | aka._pat_parsable,
                                  'Youngs modulus')
        super().registerParamReal('nu',
                                  aka._pat_readable | aka._pat_parsable,
                                  'Poisson ratio')

    def initMaterial(self):
        nu = self.getReal('nu')
        E = self.getReal('E')
        self.mu = E / (2 * (1 + nu))
        self.lame_lambda = nu * E / (
            (1. + nu) * (1. - 2. * nu))
        # Second Lame coefficient (shear modulus)
        self.lame_mu = E / (2. * (1. + nu))
        super().initMaterial()

    # declares all the parameters that are needed
    def getPushWaveSpeed(self, element):
        rho = self.getReal('rho')
        return np.sqrt((self.lame_lambda + 2 * self.lame_mu) / rho)

    # compute small deformation tensor
    @staticmethod
    def computeEpsilon(grad_u):
        return 0.5 * (grad_u + np.einsum('aij->aji', grad_u))

    # constitutive law
    def computeStress(self, el_type, ghost_type):
        grad_u = self.getGradU(el_type, ghost_type)
        sigma = self.getStress(el_type, ghost_type)

        n_quads = grad_u.shape[0]
        grad_u = grad_u.reshape((n_quads, 2, 2))
        epsilon = self.computeEpsilon(grad_u)
        sigma = sigma.reshape((n_quads, 2, 2))
        trace = np.einsum('aii->a', grad_u)

        sigma[:, :, :] = (
            np.einsum('a,ij->aij', trace,
                      self.lame_lambda * np.eye(2))
            + 2. * self.lame_mu * epsilon)

    # constitutive law tangent modulii
    def computeTangentModuli(self, el_type, tangent_matrix, ghost_type):
        n_quads = tangent_matrix.shape[0]
        tangent = tangent_matrix.reshape(n_quads, 3, 3)

        Miiii = self.lame_lambda + 2 * self.lame_mu
        Miijj = self.lame_lambda
        Mijij = self.lame_mu

        tangent[:, 0, 0] = Miiii
        tangent[:, 1, 1] = Miiii
        tangent[:, 0, 1] = Miijj
        tangent[:, 1, 0] = Miijj
        tangent[:, 2, 2] = Mijij

    # computes the energy density
    def computePotentialEnergy(self, el_type):

        sigma = self.getStress(el_type)
        grad_u = self.getGradU(el_type)

        nquads = sigma.shape[0]
        stress = sigma.reshape(nquads, 2, 2)
        grad_u = grad_u.reshape((nquads, 2, 2))
        epsilon = self.computeEpsilon(grad_u)

        energy_density = self.getPotentialEnergy(el_type)
        energy_density[:, 0] = 0.5 * np.einsum('aij,aij->a', stress, epsilon)


# register material to the MaterialFactory
def allocator(_dim, unused, model, _id):
    return LocalElastic(model, _id)


mat_factory = aka.MaterialFactory.getInstance()
mat_factory.registerAllocator("local_elastic", allocator)

# ------------------------------------------------------------------------------
# main
# ------------------------------------------------------------------------------
spatial_dimension = 2
aka.parseInput('material.dat')

mesh_file = 'bar.msh'
max_steps = 250
time_step = 1e-3

# ------------------------------------------------------------------------------
# Initialization
# ------------------------------------------------------------------------------
mesh = aka.Mesh(spatial_dimension)
mesh.read(mesh_file)

# parse input file
aka.parseInput('material.dat')

model = aka.SolidMechanicsModel(mesh)
model.initFull(_analysis_method=aka._explicit_lumped_mass)

model.setBaseName("waves")
model.addDumpFieldVector("displacement")
model.addDumpFieldVector("acceleration")
model.addDumpFieldVector("velocity")
model.addDumpFieldVector("internal_force")
model.addDumpFieldVector("external_force")
model.addDumpField("strain")
model.addDumpField("stress")
model.addDumpField("blocked_dofs")

# ------------------------------------------------------------------------------
# boundary conditions
# ------------------------------------------------------------------------------
model.applyBC(aka.FixedValue(0, aka._x), "XBlocked")
model.applyBC(aka.FixedValue(0, aka._y), "YBlocked")

# ------------------------------------------------------------------------------
# initial conditions
# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
# timestep value computation
# ------------------------------------------------------------------------------
time_factor = 0.8
stable_time_step = model.getStableTimeStep() * time_factor

print("Stable Time Step = {0}".format(stable_time_step))
print("Required Time Step = {0}".format(time_step))

time_step = stable_time_step * time_factor

model.setTimeStep(time_step)

# ------------------------------------------------------------------------------
# loop for evolution of motion dynamics
# ------------------------------------------------------------------------------
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
