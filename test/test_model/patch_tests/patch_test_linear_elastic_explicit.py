#!/usr/bin/env python3

from patch_test_linear_solid_mechanics_fixture import TestPatchTestSMMLinear
import akantu


def foo(self):

    filename = "material_check_stress_plane_stress.dat"
    if self.plane_strain:
        filename = "material_check_stress_plane_strain.dat"

    self.initModel(
        akantu.SolidMechanicsModelOptions(akantu._explicit_lumped_mass),
        filename)

    coordinates = self.mesh.getNodes()
    displacement = self.model.getDisplacement()
    # set the position of all nodes to the static solution
    self.setLinearDOF(displacement, coordinates)

    for s in range(0, 100):
            self.model.solveStep()

    ekin = self.model.getEnergy("kinetic")
    self.assertAlmostEqual(0, ekin, 1e-16)
    self.checkAll()


TestPatchTestSMMLinear.TYPED_TEST(foo, "Explicit")
