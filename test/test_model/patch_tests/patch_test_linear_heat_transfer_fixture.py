#!/usr/bin/env python3

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
        C = self.model.getParamMatrix("conductivity")
        self.checkDOFs(temperature)
        self.checkGradient(self.model.getTemperatureGradient(self.elem_type),
                           temperature)

        self.prescribed_gradient(temperature)
        self.checkResults(lambda grad_T: C.dot(grad_T.T),
                          self.model.getKgradT(self.elem_type),
                          temperature)
