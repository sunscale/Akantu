/**
 * @file   test_integrate.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Peter Spijker <peter.spijker@epfl.ch>
 *
 * @date creation: Fri Sep 03 2010
 * @date last modification: Mon Feb 19 2018
 *
 * @brief  test of the fem class
 *
 * @section LICENSE
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
#include <cstdlib>
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;

TYPED_TEST(TestFEMFixture, IntegrateConstant) {
  this->fem->initShapeFunctions();

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

  // integrate function on elements
  Array<Real> int_val_on_elem(this->nb_element, 2, "int_val_on_elem");
  this->fem->integrate(val_on_quad, int_val_on_elem, 2, type);

  // get global integration value
  Vector<Real> sum{0., 0.};

  for (auto && int_ : make_view(int_val_on_elem, 2)) {
    sum += int_;
  }

  auto diff = (value - sum).template norm<L_inf>();
  EXPECT_NEAR(0, diff, 1e-14);
}
