/**
 * @file   test_interpolate.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Sep 03 2010
 * @date last modification: Tue Nov 14 2017
 *
 * @brief  test of the fem class
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "test_fe_engine_fixture.hh"
/* -------------------------------------------------------------------------- */
using namespace akantu;

namespace {

TYPED_TEST(TestFEMFixture, InterpolateConstant) {
  const auto type = this->type;

  const auto & position = this->fem->getMesh().getNodes();

  Array<Real> const_val(position.size(), 2, "const_val");
  Array<Real> val_on_quad(this->nb_quadrature_points_total, 2, "val_on_quad");

  Vector<Real> value{1, 2};
  for (auto && const_ : make_view(const_val, 2)) {
    const_ = value;
  }

  // interpolate function on quadrature points
  this->fem->interpolateOnIntegrationPoints(const_val, val_on_quad, 2, type);

  for (auto && int_ : make_view(val_on_quad, 2)) {
    auto diff = (value - int_).template norm<L_inf>();
    EXPECT_NEAR(0, diff, 1e-14);
  }
}

// TYPED_TEST(TestFEMFixture, InterpolatePosition) {
//   const auto dim = this->dim;
//   const auto type = this->type;
//   const auto & position = this->fem->getMesh().getNodes();

//   Array<Real> coord_on_quad(this->nb_quadrature_points_total, dim,
//                             "coord_on_quad");

//   this->fem->interpolateOnIntegrationPoints(position, coord_on_quad, dim,
//   type);
// }

} // namespace
