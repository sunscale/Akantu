#!/usr/bin/env python3

import patch_test_linear_fixture
import numpy as np
import akantu


# custom material (this patch test also checks for custom material features)
class LocalElastic:

    # declares all the internals
    def initMaterial(self, internals, params):
        self.E = params['E']
        self.nu = params['nu']
        self.rho = params['rho']
        print(self.__dict__)
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
            factor[quad_coords[:, 1] < 0.5] = 1.

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
        epsilon = self.computeEpsilon(grad_u)
        sigma = sigma.reshape((n_quads, 2, 2))
        trace = np.einsum('aii,aii->a', grad_u, grad_u)

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


class TestPatchTestSMMLinear(patch_test_linear_fixture.TestPatchTestLinear):

    plane_strain = True
    model_type = akantu.SolidMechanicsModel

    def __init__(self, *args, **kwargs):
        mat = LocalElastic()
        # try:
        #     akantu.registerNewPythonMaterial(mat, "local_elastic")
        # except Exception as e:
        #     print(e)
        #     pass
        super().__init__(*args, **kwargs)

    def initModel(self, method, material_file):
        super().initModel(method, material_file)

    def applyBC(self):
        super().applyBC()
        displacement = self.model.getDisplacement()
        self.applyBConDOFs(displacement)

    def checkAll(self):
        displacement = self.model.getDisplacement()
        mat = self.model.getMaterial(0)

        self.checkDOFs(displacement)
        self.checkGradient(mat.getGradU(self.elem_type), displacement)

        def foo(pstrain):
            nu = self.model.getMaterial(0).getParamReal("nu")
            E = self.model.getMaterial(0).getParamReal("E")

            strain = (pstrain + pstrain.transpose()) / 2.
            trace = strain.trace()

            _lambda = nu * E / ((1 + nu) * (1 - 2 * nu))
            mu = E / (2 * (1 + nu))

            if (not self.plane_strain):
                _lambda = nu * E / (1 - nu * nu)

            stress = np.zeros((self.dim, self.dim))

            if self.dim == 1:
                stress[0, 0] = E * strain[0, 0]
            else:
                stress[:, :] = (
                    _lambda * trace * np.eye(self.dim) + 2 * mu * strain)

            return stress

        self.checkResults(foo,
                          mat.getStress(self.elem_type),
                          displacement)
