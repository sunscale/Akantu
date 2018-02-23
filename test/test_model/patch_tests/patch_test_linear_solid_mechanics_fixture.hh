/* -------------------------------------------------------------------------- */
#include "patch_test_linear_fixture.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_PATCH_TEST_LINEAR_SOLID_MECHANICS_FIXTURE_HH__
#define __AKANTU_PATCH_TEST_LINEAR_SOLID_MECHANICS_FIXTURE_HH__

/* -------------------------------------------------------------------------- */
template <typename tuple_>
class TestPatchTestSMMLinear
    : public TestPatchTestLinear<std::tuple_element_t<0, tuple_>,
                                 SolidMechanicsModel> {
  using parent =
      TestPatchTestLinear<std::tuple_element_t<0, tuple_>, SolidMechanicsModel>;

public:
  static constexpr bool plane_strain = std::tuple_element_t<1, tuple_>::value;

  void applyBC() override {
    parent::applyBC();
    auto & displacement = this->model->getDisplacement();
    this->applyBConDOFs(displacement);
  }

  void checkAll() {
    auto & displacement = this->model->getDisplacement();
    auto & mat = this->model->getMaterial(0);

    this->checkDOFs(displacement);
    this->checkGradient(mat.getGradU(this->type), displacement);
    this->checkResults(
        [&](const Matrix<Real> & pstrain) {
          Real nu = this->model->getMaterial(0).get("nu");
          Real E = this->model->getMaterial(0).get("E");

          auto strain = (pstrain + pstrain.transpose()) / 2.;
          auto trace = strain.trace();

          auto lambda = nu * E / ((1 + nu) * (1 - 2 * nu));
          auto mu = E / (2 * (1 + nu));

          if (not this->plane_strain) {
            lambda = nu * E / (1 - nu * nu);
          }

          decltype(strain) stress(this->dim, this->dim);

          if (this->dim == 1) {
            stress(0, 0) = E * strain(0, 0);
          } else {
            for (UInt i = 0; i < this->dim; ++i)
              for (UInt j = 0; j < this->dim; ++j)
                stress(i, j) =
                    (i == j) * lambda * trace + 2 * mu * strain(i, j);
          }

          return stress;
        },
        mat.getStress(this->type), displacement);
  }
};

template <typename tuple_>
constexpr bool TestPatchTestSMMLinear<tuple_>::plane_strain;

template <typename T> struct invalid_plan_stress : std::true_type {};
template <typename type, typename bool_c>
struct invalid_plan_stress<std::tuple<type, bool_c>>
    : aka::bool_constant<ElementClass<type::value>::getSpatialDimension() !=
                             2 and
                         not bool_c::value> {};

using true_false =
    std::tuple<aka::bool_constant<true>, aka::bool_constant<false>>;

template <typename T> using valid_types = aka::negation<invalid_plan_stress<T>>;

using types = gtest_list_t<
    tuple_filter_t<valid_types, cross_product_t<TestElementTypes, true_false>>>;

TYPED_TEST_CASE(TestPatchTestSMMLinear, types);

#endif /* __AKANTU_PATCH_TEST_LINEAR_SOLID_MECHANICS_FIXTURE_HH__ */
