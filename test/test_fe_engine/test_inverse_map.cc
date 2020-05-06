/**
 * @file   test_inverse_map.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Fri Sep 03 2010
 * @date last modification: Mon Feb 19 2018
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

TYPED_TEST(TestFEMFixture, InverseMap) {
  this->fem->initShapeFunctions();

  Matrix<Real> quad =
      GaussIntegrationElement<TestFixture::type>::getQuadraturePoints();

  const auto & position = this->fem->getMesh().getNodes();

  /// get the quadrature points coordinates
  Array<Real> coord_on_quad(quad.cols() * this->nb_element, this->dim,
                            "coord_on_quad");

  this->fem->interpolateOnIntegrationPoints(position, coord_on_quad, this->dim,
                                            this->type);

  Vector<Real> natural_coords(this->dim);

  auto length = (this->upper - this->lower).template norm<L_inf>();

  for (auto && enum_ :
       enumerate(make_view(coord_on_quad, this->dim, quad.cols()))) {
    auto el = std::get<0>(enum_);
    const auto & quads_coords = std::get<1>(enum_);

    for (auto q : arange(quad.cols())) {
      Vector<Real> quad_coord = quads_coords(q);
      Vector<Real> ref_quad_coord = quad(q);
      this->fem->inverseMap(quad_coord, el, this->type, natural_coords);

      auto dis_normalized = ref_quad_coord.distance(natural_coords) / length;
      EXPECT_NEAR(0., dis_normalized, 3.5e-11);
    }
  }
}
