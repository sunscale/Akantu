#!/usr/bin/env python3

import patch_test_linear_fixture
import numpy as np
import akantu


class TestPatchTestSMMLinear(patch_test_linear_fixture.TestPatchTestLinear):

    plane_strain = True
    model_type = akantu.SolidMechanicsModel

    def applyBC(self):
        patch_test_linear_fixture.TestPatchTestLinear.applyBC(self)
        displacement = self.model.getDisplacement()
        self.applyBConDOFs(displacement)

    def checkAll(self):
        displacement = self.model.getDisplacement()
        mat = self.model.getMaterial(0)

        self.checkDOFs(displacement)
        self.checkGradient(mat.getGradU(self.elem_type), displacement)

        def foo(pstrain):
            nu = self.model.getMaterial(0).getParamReal("nu")
            E = self.model.getMaterial(0).getParamReal("E")

            strain = (pstrain + pstrain.transpose()) / 2.
            trace = strain.trace()

            _lambda = nu * E / ((1 + nu) * (1 - 2 * nu))
            mu = E / (2 * (1 + nu))

            if (not self.plane_strain):
                _lambda = nu * E / (1 - nu * nu)

            stress = np.zeros((self.dim, self.dim))

            if self.dim == 1:
                stress[0, 0] = E * strain[0, 0]
            else:
                stress[:, :] = (
                    _lambda * trace * np.eye(self.dim) + 2 * mu * strain)

            return stress

        self.checkResults(foo,
                          mat.getStress(self.elem_type),
                          displacement)


# template <typename T> struct invalid_plan_stress : std::true_type {}
# template <typename type, typename bool_c>
# struct invalid_plan_stress<std::tuple<type, bool_c>>
#     : aka::bool_constant<ElementClass<type::value>::getSpatialDimension() !=
#                              2 and
#                          not bool_c::value> {}
# 
# using true_false =
#     std::tuple<aka::bool_constant<true>, aka::bool_constant<false>>
# 
# template <typename T> using valid_types = aka::negation<invalid_plan_stress<T>>
# 
# using types = gtest_list_t<
#     tuple_filter_t<valid_types, cross_product_t<TestElementTypes, true_false>>>
# 
# TYPED_TEST_CASE(TestPatchTestSMMLinear, types)

#endif /* __AKANTU_PATCH_TEST_LINEAR_SOLID_MECHANICS_FIXTURE_HH__ */
