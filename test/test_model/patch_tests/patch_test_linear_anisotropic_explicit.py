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

from patch_test_linear_solid_mechanics_fixture import TestPatchTestSMMLinear
import akantu
import unittest
import numpy as np

# Stiffness tensor, rotated by hand
C = np.array([[[[112.93753505, 1.85842452538e-10, -4.47654358027e-10],
               [1.85847317471e-10, 54.2334345331, -3.69840984824],
               [-4.4764768395e-10, -3.69840984824, 56.848605217]],
              [[1.85847781609e-10, 25.429294233, -3.69840984816],
               [25.429294233, 3.31613847493e-10, -8.38797920011e-11],
               [-3.69840984816, -8.38804581349e-11, -1.97875715813e-10]],
              [[-4.47654358027e-10, -3.69840984816, 28.044464917],
               [-3.69840984816, 2.09374961813e-10, 9.4857455224e-12],
               [28.044464917, 9.48308098714e-12, -2.1367885239e-10]]],
             [[[1.85847781609e-10, 25.429294233, -3.69840984816],
               [25.429294233, 3.31613847493e-10, -8.38793479119e-11],
               [-3.69840984816, -8.38795699565e-11, -1.97876381947e-10]],
              [[54.2334345331, 3.31617400207e-10, 2.09372075233e-10],
               [3.3161562385e-10, 115.552705733, -3.15093728886e-10],
               [2.09372075233e-10, -3.15090176173e-10, 54.2334345333]],
              [[-3.69840984824, -8.38795699565e-11, 9.48219280872e-12],
               [-8.38795699565e-11, -3.1509195253e-10, 25.4292942335],
               [9.48441325477e-12, 25.4292942335, 3.69840984851]]],
             [[[-4.47653469848e-10, -3.69840984816, 28.044464917],
               [-3.69840984816, 2.09374073634e-10, 9.48752187924e-12],
               [28.044464917, 9.48552347779e-12, -2.1367885239e-10]],
              [[-3.69840984824, -8.3884899027e-11, 9.48219280872e-12],
               [-8.3884899027e-11, -3.150972816e-10, 25.4292942335],
               [9.48041645188e-12, 25.4292942335, 3.69840984851]],
              [[56.848605217, -1.97875493768e-10, -2.13681516925e-10],
               [-1.97877270125e-10, 54.2334345333, 3.69840984851],
               [-2.13683293282e-10, 3.69840984851, 112.93753505]]]])


def foo(self):

    self.initModel(
        akantu.SolidMechanicsModelOptions(akantu._explicit_lumped_mass),
        "material_anisotropic.dat")

    coordinates = self.mesh.getNodes()
    displacement = self.model.getDisplacement()

    # set the position of all nodes to the static solution
    self.setLinearDOF(displacement, coordinates)

    for s in range(0, 100):
        self.model.solveStep()

    ekin = self.model.getEnergy("kinetic")
    self.assertAlmostEqual(0, ekin, delta=1e-16)

    self.checkDisplacements()
    self.checkStrains()

    def foo(pstrain):
        strain = (pstrain + pstrain.transpose()) / 2.
        stress = np.zeros((self.dim, self.dim))

        for i in range(0, self.dim):
            for j in range(0, self.dim):
                stress[i, j] = 0
                for k in range(0, self.dim):
                    for l in range(0, self.dim):
                        stress[i, j] += C[i][j][k][l] * strain(k, l)
        return stress
    self.checkStresses(foo)


akantu.initialize()
suite = TestPatchTestSMMLinear.TYPED_TEST(foo, "AnisotropicExplicit")
unittest.TextTestRunner(verbosity=1).run(suite)
