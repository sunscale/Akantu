#!/usr/bin/env python3

# ------------------------------------------------------------------------------
__author__ = "Guillaume Anciaux"
__copyright__ = "Copyright (C) 2016-2018, EPFL (Ecole Polytechnique Fédérale" \
                " de Lausanne) Laboratory (LSMS - Laboratoire de Simulation" \
                " en Mécanique des Solides)"
__credits__ = ["Guillaume Anciaux"]
__license__ = "L-GPLv3"
__maintainer__ = "Guillaume Anciaux"
__email__ = "guillaume.anciaux@epfl.ch"
# ------------------------------------------------------------------------------

import patch_test_linear_fixture
import numpy as np
import akantu


# custom material (this patch test also checks for custom material features)
class LocalElastic(akantu.Material):

    def __init__(self, model, _id):
        super().__init__(model, _id)
        super().registerParamReal('E',
                                  akantu._pat_readable | akantu._pat_parsable,
                                  'Youngs modulus')
        super().registerParamReal('nu',
                                  akantu._pat_readable | akantu._pat_parsable,
                                  'Poisson ratio')

    # declares all the internals
    def initMaterial(self, internals, params):
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


class TestPatchTestSMMLinear(patch_test_linear_fixture.TestPatchTestLinear):

    plane_strain = True
    model_type = akantu.SolidMechanicsModel

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def initModel(self, method, material_file):
        # mat.__dict__['dim'] = self.dim
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
            nu = self.model.getMaterial(0).getReal("nu")
            E = self.model.getMaterial(0).getReal("E")

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
