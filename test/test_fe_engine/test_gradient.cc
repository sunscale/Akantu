/**
 * @file   test_gradient.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Peter Spijker <peter.spijker@epfl.ch>
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
 * @section DESCRIPTION
 *
 * This code is computing the gradient of a linear field and check that it gives
 * a constant result.  It also compute the gradient the  coordinates of the mesh
 * and check that it gives the identity
 *
 */

/* -------------------------------------------------------------------------- */
#include "test_fe_engine_fixture.hh"
/* -------------------------------------------------------------------------- */
#include <cstdlib>
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;

TYPED_TEST(TestFEMFixture, GradientPoly) {
  this->fem->initShapeFunctions();
  Real alpha[2][3] = {{13, 23, 31}, {11, 7, 5}};

  const auto dim = this->dim;
  const auto type = this->type;

  const auto & position = this->fem->getMesh().getNodes();
  Array<Real> const_val(this->fem->getMesh().getNbNodes(), 2, "const_val");
  for (auto && pair : zip(make_view(position, dim), make_view(const_val, 2))) {
    auto & pos = std::get<0>(pair);
    auto & const_ = std::get<1>(pair);

    const_.set(0.);

    for (UInt d = 0; d < dim; ++d) {
      const_(0) += alpha[0][d] * pos(d);
      const_(1) += alpha[1][d] * pos(d);
    }
  }

  /// compute the gradient
  Array<Real> grad_on_quad(this->nb_quadrature_points_total, 2 * dim,
                           "grad_on_quad");
  this->fem->gradientOnIntegrationPoints(const_val, grad_on_quad, 2, type);

  /// check the results
  for (auto && grad : make_view(grad_on_quad, 2, dim)) {
    for (UInt d = 0; d < dim; ++d) {
      EXPECT_NEAR(grad(0, d), alpha[0][d], 5e-13);
      EXPECT_NEAR(grad(1, d), alpha[1][d], 5e-13);
    }
  }
}

TYPED_TEST(TestFEMFixture, GradientPositions) {
  this->fem->initShapeFunctions();
  const auto dim = this->dim;
  const auto type = this->type;

  UInt nb_quadrature_points =
      this->fem->getNbIntegrationPoints(type) * this->nb_element;
  Array<Real> grad_coord_on_quad(nb_quadrature_points, dim * dim,
                                 "grad_coord_on_quad");

  const auto & position = this->mesh->getNodes();
  this->fem->gradientOnIntegrationPoints(position, grad_coord_on_quad, dim,
                                         type);

  auto I = Matrix<Real>::eye(UInt(dim));

  for (auto && grad : make_view(grad_coord_on_quad, dim, dim)) {
    auto diff = (I - grad).template norm<L_inf>();

    EXPECT_NEAR(0., diff, 2e-14);
  }
}
