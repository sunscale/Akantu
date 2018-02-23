#!/usr/bin/env python3

from patch_test_linear_solid_mechanics_fixture import TestPatchTestSMMLinear
import akantu


def foo(self):

    filename = "material_check_stress_plane_stress.dat"
    if self.plane_strain:
        filename = "material_check_stress_plane_strain.dat"

    self.initModel(akantu.SolidMechanicsModelOptions(akantu._static), filename)

    solver = self.model.getNonLinearSolver()
    solver.set("max_iterations", 2)
    solver.set("threshold", 2e-4)
    solver.set("convergence_type", akantu._scc_residual)

    self.model.solveStep()

    self.checkAll()


TestPatchTestSMMLinear.TYPED_TEST(foo, "Static")
