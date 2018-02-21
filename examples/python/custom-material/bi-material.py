from __future__ import print_function
# ------------------------------------------------------------- #
import akantu
import subprocess
import numpy as np
import time
# ------------------------------------------------------------- #


class LocalElastic:

    # declares all the internals
    def initMaterial(self, internals, params):
        self.E = params['E']
        self.nu = params['nu']
        self.rho = params['rho']
        # print(self.__dict__)
        # First Lame coefficient
        self.lame_lambda = self.nu * self.E / (
            (1. + self.nu) * (1. - 2. * self.nu))
        # Second Lame coefficient (shear modulus)
        self.lame_mu = self.E / (2. * (1. + self.nu))

        all_factor = internals['factor']
        all_quad_coords = internals['quad_coordinates']

        for elem_type in all_factor.keys():
            factor = all_factor[elem_type]
            quad_coords = all_quad_coords[elem_type]

            factor[:] = 1.
            factor[quad_coords[:, 1] < 0.5] = .5

    # declares all the internals
    @staticmethod
    def registerInternals():
        return ['potential', 'factor']

    # declares all the internals
    @staticmethod
    def registerInternalSizes():
        return [1, 1]

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
        n_quads = grad_u.shape[0]
        grad_u = grad_u.reshape((n_quads, 2, 2))
        factor = internals['factor'].reshape(n_quads)
        epsilon = self.computeEpsilon(grad_u)
        sigma = sigma.reshape((n_quads, 2, 2))
        trace = np.einsum('aii->a', grad_u)

        sigma[:, :, :] = (
            np.einsum('a,ij->aij', trace,
                      self.lame_lambda * np.eye(2))
            + 2. * self.lame_mu * epsilon)

        # print(sigma.reshape((n_quads, 4)))
        # print(grad_u.reshape((n_quads, 4)))
        sigma[:, :, :] = np.einsum('aij, a->aij', sigma, factor)

    # constitutive law tangent modulii
    def computeTangentModuli(self, grad_u, tangent, internals, params):
        n_quads = tangent.shape[0]
        tangent = tangent.reshape(n_quads, 3, 3)
        factor = internals['factor'].reshape(n_quads)

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


################################################################
# main
################################################################

# main paraeters
spatial_dimension = 2
mesh_file = 'square.msh'

# call gmsh to generate the mesh
ret = subprocess.call(['gmsh', '-2', 'square.geo', '-optimize', 'square.msh'])
if ret != 0:
    raise Exception(
        'execution of GMSH failed: do you have it installed ?')

time.sleep(1)

# read mesh
mesh = akantu.Mesh(spatial_dimension)
mesh.read(mesh_file)

# create the custom material
mat = LocalElastic()
akantu.registerNewPythonMaterial(mat, "local_elastic")

# init akantu
akantu.initialize('material.dat')

# init the SolidMechanicsModel
model = akantu.SolidMechanicsModel(mesh)
model.initFull(_analysis_method=akantu._static)

# configure the solver
solver = model.getNonLinearSolver()
solver.set("max_iterations", 2)
solver.set("threshold", 1e-3)
solver.set("convergence_type", akantu._scc_solution)

# prepare the dumper
model.setBaseName("bimaterial")
model.addDumpFieldVector("displacement")
model.addDumpFieldVector("internal_force")
model.addDumpFieldVector("external_force")
model.addDumpField("strain")
model.addDumpField("stress")
model.addDumpField("factor")
model.addDumpField("blocked_dofs")

# Boundary conditions
applyBC(model)

# solve the problem
model.solveStep()

# dump paraview files
model.dump()

# shutdown akantu
akantu.finalize()
