/**
 * @file   patch_test_linear_solid_mechanics_fixture.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Jan 30 2018
 *
 * @brief  SolidMechanics patch tests fixture
 *
 *
 * Copyright (©) 2016-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

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

  void checkForces() {
    auto & mat = this->model->getMaterial(0);
    auto & internal_forces = this->model->getInternalForce();
    auto & external_forces = this->model->getExternalForce();
    auto dim = this->dim;

    Matrix<Real> sigma =
        make_view(mat.getStress(this->type), dim, dim).begin()[0];

    external_forces.clear();
    if (dim > 1) {
      for (auto & eg : this->mesh->iterateElementGroups()) {
        this->model->applyBC(BC::Neumann::FromHigherDim(sigma), eg.getName());
      }
    } else {
      external_forces(0) = -sigma(0, 0);
      external_forces(1) = sigma(0, 0);
    }

    Real force_norm_inf = -std::numeric_limits<Real>::max();

    Vector<Real> total_force(dim);
    total_force.clear();

    for (auto && f : make_view(internal_forces, dim)) {
      total_force += f;
      force_norm_inf = std::max(force_norm_inf, f.template norm<L_inf>());
    }

    EXPECT_NEAR(0, total_force.template norm<L_inf>() / force_norm_inf, 1e-9);

    for (auto && tuple : zip(make_view(internal_forces, dim),
                             make_view(external_forces, dim))) {
      auto && f_int = std::get<0>(tuple);
      auto && f_ext = std::get<1>(tuple);
      auto f = f_int + f_ext;
      EXPECT_NEAR(0, f.template norm<L_inf>() / force_norm_inf, 1e-9);
    }
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
    this->checkForces();
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

using model_types = gtest_list_t<
    tuple_filter_t<valid_types, cross_product_t<TestElementTypes, true_false>>>;

TYPED_TEST_SUITE(TestPatchTestSMMLinear, model_types);

#endif /* __AKANTU_PATCH_TEST_LINEAR_SOLID_MECHANICS_FIXTURE_HH__ */
