/**
 * @file   test_inverse_map.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Fri Sep 03 2010
 * @date last modification: Thu Oct 15 2015
 *
 * @brief  test of the fem class
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
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

TYPED_TEST(TestFEMFixture, InverseMap) {
  const auto type = this->type;
  const auto dim = this->dim;

  Matrix<Real> quad = GaussIntegrationElement<type>::getQuadraturePoints();

  const auto & position = this->fem->getMesh().getNodes();

  /// get the quadrature points coordinates
  Array<Real> coord_on_quad(quad.cols() * this->nb_element, dim,
                            "coord_on_quad");

  this->fem->interpolateOnIntegrationPoints(position, coord_on_quad, dim, type);

  Vector<Real> natural_coords(dim);

  auto length = (this->upper - this->lower).template norm<L_inf>();

  for(auto && enum_ : enumerate(make_view(coord_on_quad, dim, quad.cols()))) {
    auto el = std::get<0>(enum_);
    const auto & quads_coords = std::get<1>(enum_);

    for (auto q : arange(quad.cols())) {
      Vector<Real> quad_coord = quads_coords(q);
      Vector<Real> ref_quad_coord = quad(q);
      this->fem->inverseMap(quad_coord, el, type, natural_coords);

      auto dis_normalized = ref_quad_coord.distance(natural_coords) / length;
      EXPECT_NEAR(0., dis_normalized, 5e-12);
    }
  }
}

} // namespace
