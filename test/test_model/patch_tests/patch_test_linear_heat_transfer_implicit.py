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

from patch_test_linear_heat_transfer_fixture import TestPatchTestHTMLinear
import akantu


def foo(self):
    self.initModel(akantu.HeatTransferModelOptions(akantu._implicit_dynamic),
                   "heat_transfer_input.dat")

    coordinates = self.mesh.getNodes()
    temperature = self.model.getTemperature()
    #  set the position of all nodes to the static solution
    self.setLinearDOF(temperature, coordinates)

    for s in range(0, 100):
        self.model.solveStep()

    self.checkAll()


TestPatchTestHTMLinear.TYPED_TEST(foo, "Explicit")
