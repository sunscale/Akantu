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
from mpi4py import MPI
import akantu
import unittest
import numpy as np


class TestPatchTestLinear(unittest.TestCase):

    alpha = np.array([[0.01, 0.02, 0.03, 0.04],
                      [0.05, 0.06, 0.07, 0.08],
                      [0.09, 0.10, 0.11, 0.12]])

    gradient_tolerance = 1e-13
    result_tolerance = 1e-13
    dofs_tolerance = 1e-15

    def __init__(self, test_name, elem_type_str, functor=None):

        self.test_name = test_name
        self.functor = functor
        self.elem_type = akantu.getElementTypes()[elem_type_str]
        self.elem_type_str = elem_type_str
        self.dim = akantu.Mesh.getSpatialDimension(self.elem_type)

        super().__init__(test_name)

    def _internal_call(self):
        self.functor(self)

    def __getattr__(self, key):
        if key == self.test_name:
            return self._internal_call

    def setUp(self):
        self.mesh = akantu.Mesh(self.dim, self.elem_type_str)
        self.mesh.read(str(self.elem_type_str) + ".msh")
        akantu.MeshUtils.buildFacets(self.mesh)
        self.mesh.createBoundaryGroupFromGeometry()
        self.model = self.model_type(self.mesh)

    def tearDown(self):
        del self.model
        del self.mesh

    def initModel(self, method, material_file):
        akantu.parseInput(material_file)
        akantu.setDebugLevel(akantu.dblError)
        self.model.initFull(method)
        self.applyBC()

        if method != akantu._static:
            self.model.setTimeStep(0.8 * self.model.getStableTimeStep())

    def applyBC(self):
        boundary = self.model.getBlockedDOFs()
        for eg in self.mesh.iterateElementGroups():
            nodes = eg.getNodeGroup().getNodes()
            boundary[nodes, :] = True

    def applyBConDOFs(self, dofs):
        coordinates = self.mesh.getNodes()

        for eg in self.mesh.iterateElementGroups():
            nodes = eg.getNodeGroup().getNodes().flatten()
            dofs[nodes] = self.setLinearDOF(dofs[nodes],
                                            coordinates[nodes])

    def prescribed_gradient(self, dof):
        gradient = self.alpha[:dof.shape[1], 1:self.dim + 1]
        return gradient

    def checkGradient(self, gradient, dofs):
        pgrad = self.prescribed_gradient(dofs).T
        gradient = gradient.reshape(gradient.shape[0], *pgrad.shape)
        diff = gradient[:] - pgrad
        norm = np.abs(pgrad).max()
        gradient_error = np.abs(diff).max() / norm
        self.assertAlmostEqual(0, gradient_error,
                               delta=self.gradient_tolerance)

    def checkResults(self, presult_func, results, dofs):
        presult = presult_func(self.prescribed_gradient(dofs)).flatten()
        remaining_size = np.prod(np.array(results.shape[1:]))
        results = results.reshape((results.shape[0], remaining_size))

        for result in results:
            diff = result - presult
            norm = np.abs(result).max()
            if norm == 0:
                result_error = np.abs(diff).max()
            else:
                result_error = np.abs(diff).max() / norm

            self.assertAlmostEqual(0., result_error,
                                   delta=self.result_tolerance)

    def setLinearDOF(self, dof, coord):
        nb_dofs = dof.shape[1]
        dof[:] = np.einsum('ik,ak->ai',
                           self.alpha[:nb_dofs, 1:self.dim + 1], coord)
        for i in range(0, nb_dofs):
            dof[:, i] += self.alpha[i, 0]

        return dof

    def checkDOFs(self, dofs):
        coordinates = self.mesh.getNodes()
        ref_dofs = np.zeros_like(dofs)
        self.setLinearDOF(ref_dofs, coordinates)
        diff = dofs - ref_dofs
        dofs_error = np.abs(diff).max()
        self.assertAlmostEqual(0., dofs_error, delta=self.dofs_tolerance)

    @classmethod
    def TYPED_TEST(cls, functor, label):
        for type_name, _type in akantu.getElementTypes().items():
            if type_name == "_point_1":
                continue

            method_name = type_name + '_' + label
            test_case = cls(method_name, type_name, functor)
            test_case.setUp()
            functor(test_case)
            test_case.tearDown()
