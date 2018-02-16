#!/usr/bin/env python3

from patch_test_linear_heat_transfer_fixture import TestPatchTestHTMLinear
import akantu


def foo(self):

    self.initModel(
        akantu.HeatTransferModelOptions(akantu._explicit_lumped_mass),
        "heat_transfer_input.dat")

    coordinates = self.mesh.getNodes()
    temperature = self.model.getTemperature()
    # set the position of all nodes to the static solution
    self.setLinearDOF(temperature, coordinates)

    for s in range(0, 100):
        self.model.solveStep()

    self.checkAll()


TestPatchTestHTMLinear.TYPED_TEST(foo, "Explicit")
