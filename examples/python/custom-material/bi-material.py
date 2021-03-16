import akantu as aka
import numpy as np


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

        # change it to have the initialize wrapped
        super().registerInternalReal('factor', 1)
        super().registerInternalReal('quad_coordinates', 2)

    def initMaterial(self):
        nu = self.getReal('nu')
        E = self.getReal('E')
        self.mu = E / (2 * (1 + nu))
        self.lame_lambda = nu * E / (
            (1. + nu) * (1. - 2. * nu))
        # Second Lame coefficient (shear modulus)
        self.lame_mu = E / (2. * (1. + nu))
        super().initMaterial()

        quad_coords = self.getInternalReal("quad_coordinates")
        factor = self.getInternalReal("factor")
        model = self.getModel()

        model.getFEEngine().computeIntegrationPointsCoordinates(
            quad_coords, self.getElementFilter())

        for elem_type in factor.elementTypes():
            factor = factor(elem_type)
            coords = quad_coords(elem_type)

            factor[:] = 1.
            factor[coords[:, 1] < 0.5] = .5

    # declares all the parameters that are needed
    def getPushWaveSpeed(self, params):
        return np.sqrt((self.lame_lambda + 2 * self.lame_mu) / self.rho)

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
        factor = self.getInternalReal('factor')(
            el_type, ghost_type).reshape(n_quads)
        epsilon = self.computeEpsilon(grad_u)
        sigma = sigma.reshape((n_quads, 2, 2))
        trace = np.einsum('aii->a', grad_u)

        sigma[:, :, :] = (
            np.einsum('a,ij->aij', trace,
                      self.lame_lambda * np.eye(2))
            + 2. * self.lame_mu * epsilon)

        sigma[:, :, :] = np.einsum('aij, a->aij', sigma, factor)

    # constitutive law tangent modulii
    def computeTangentModuli(self, el_type, tangent_matrix, ghost_type):
        n_quads = tangent_matrix.shape[0]
        tangent = tangent_matrix.reshape(n_quads, 3, 3)
        factor = self.getInternalReal('factor')(
            el_type, ghost_type).reshape(n_quads)

        Miiii = self.lame_lambda + 2 * self.lame_mu
        Miijj = self.lame_lambda
        Mijij = self.lame_mu

        tangent[:, 0, 0] = Miiii
        tangent[:, 1, 1] = Miiii
        tangent[:, 0, 1] = Miijj
        tangent[:, 1, 0] = Miijj
        tangent[:, 2, 2] = Mijij
        tangent[:, :, :] = np.einsum('aij, a->aij', tangent, factor)

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


# ------------------------------------------------------------------------------
# applies manually the boundary conditions
def applyBC(model):

    nbNodes = model.getMesh().getNbNodes()
    position = model.getMesh().getNodes()
    displacement = model.getDisplacement()
    blocked_dofs = model.getBlockedDOFs()

    width = 1.
    height = 1.
    epsilon = 1e-8
    for node in range(0, nbNodes):
        if((np.abs(position[node, 0]) < epsilon) or  # left side
           (np.abs(position[node, 0] - width) < epsilon)):  # right side
            blocked_dofs[node, 0] = True
            displacement[node, 0] = 0 * position[node, 0] + 0.

        if(np.abs(position[node, 1]) < epsilon):  # lower side
            blocked_dofs[node, 1] = True
            displacement[node, 1] = - 1.
        if(np.abs(position[node, 1] - height) < epsilon):  # upper side
            blocked_dofs[node, 1] = True
            displacement[node, 1] = 1.


# register the material to the material factory
def allocator(dim, option, model, id):
    return LocalElastic(model, id)


mat_factory = aka.MaterialFactory.getInstance()
mat_factory.registerAllocator("local_elastic", allocator)

# main parameters
spatial_dimension = 2
mesh_file = 'square.msh'

# read mesh
mesh = aka.Mesh(spatial_dimension)
mesh.read(mesh_file)

# parse input file
aka.parseInput('material.dat')

# init the SolidMechanicsModel
model = aka.SolidMechanicsModel(mesh)
model.initFull(_analysis_method=aka._static)

# configure the solver
solver = model.getNonLinearSolver()
solver.set("max_iterations", 2)
solver.set("threshold", 1e-3)
solver.set("convergence_type", aka.SolveConvergenceCriteria.solution)

# prepare the dumper
model.setBaseName("bimaterial")
model.addDumpFieldVector("displacement")
model.addDumpFieldVector("internal_force")
model.addDumpFieldVector("external_force")
model.addDumpField("strain")
model.addDumpField("stress")
# model.addDumpField("factor")
model.addDumpField("blocked_dofs")

# Boundary conditions
applyBC(model)

# solve the problem
model.solveStep()

# dump paraview files
model.dump()

epot = model.getEnergy('potential')
print('Potential energy: ' + str(epot))
