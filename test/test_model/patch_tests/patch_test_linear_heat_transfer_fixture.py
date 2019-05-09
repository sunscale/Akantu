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
import akantu


class TestPatchTestHTMLinear(patch_test_linear_fixture.TestPatchTestLinear):

    model_type = akantu.HeatTransferModel

    def applyBC(self):
        super().applyBC()
        temperature = self.model.getTemperature()
        self.applyBConDOFs(temperature)

    def checkAll(self):
        temperature = self.model.getTemperature()
        C = self.model.getMatrix("conductivity")
        self.checkDOFs(temperature)
        self.checkGradient(self.model.getTemperatureGradient(self.elem_type),
                           temperature)

        self.prescribed_gradient(temperature)
        self.checkResults(lambda grad_T: C.dot(grad_T.T),
                          self.model.getKgradT(self.elem_type),
                          temperature)

    def initModel(self, method, material_file):
        super().initModel(method, material_file)

        if method != akantu._static:
            self.model.setTimeStep(0.5 * self.model.getStableTimeStep())
