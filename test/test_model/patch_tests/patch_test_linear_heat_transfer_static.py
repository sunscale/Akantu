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
    self.initModel(akantu.HeatTransferModelOptions(akantu._static),
                   "heat_transfer_input.dat")

    solver = self.model.getNonLinearSolver()
    solver.set("max_iterations", 2)
    solver.set("threshold", 2e-4)
    solver.set("convergence_type", akantu._scc_residual)

    self.model.solveStep()

    self.checkAll()


TestPatchTestHTMLinear.TYPED_TEST(foo, "Static")
