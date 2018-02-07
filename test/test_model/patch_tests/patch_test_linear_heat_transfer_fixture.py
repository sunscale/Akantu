#!/usr/bin/env python3

import patch_test_linear_fixture


class TestPatchTestHTMLinear(patch_test_linear_fixture.TestPatchTestLinear):

    def applyBC(self):
        super().applyBC()
        temperature = self.model.getTemperature()
        self.applyBConDOFs(temperature)

    def checkAll(self):
        temperature = self.model.getTemperature()
        C = self.model.getParamReal("conductivity")
        self.checkDOFs(temperature)
        self.checkGradient(self.model.getTemperatureGradient(self.type_elem),
                           temperature)

        self.checkResults(lambda grad_T: C * grad_T.T,
                          self.model.getKgradT(self.type_elem), temperature)
