#!/usr/bin/env python3

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


TestPatchTestHTMLinear.TYPED_TEST(foo, "Implicit")
