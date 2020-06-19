/**
 * @file   test_fe_engine_precomputation_bernoulli_2.cc
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Jun 14 2010
 * @date last modification: Thu Jan 25 2018
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
#include "fe_engine.hh"
#include "integrator_gauss.hh"
#include "shape_structural.hh"
/* -------------------------------------------------------------------------- */
#include <cmath>
#include <functional>
#include <iostream>
/* -------------------------------------------------------------------------- */
using namespace akantu;

Matrix<Real> globalToLocalRotation(Real theta) {
  // clang-format off
  return {{ std::cos(theta), std::sin(theta), 0},
          {-std::sin(theta), std::cos(theta), 0},
          {               0,               0, 1}};
  // clang-format on
}

Vector<Real> axialReference(Real xq) {
  return {(1. - xq) / 2, 0, 0, (1. + xq) / 2, 0, 0};
}

Vector<Real> bendingReference(Real xq) {
  return {0,
          1. / 4. * Math::pow<2>(xq - 1) * (xq + 2),
          1. / 4. * Math::pow<2>(xq - 1) * (xq + 1),
          0,
          1. / 4. * Math::pow<2>(xq + 1) * (2 - xq),
          1. / 4. * Math::pow<2>(xq + 1) * (xq - 1)};
}

Vector<Real> bendingRotationReference(Real xq) {
  return {0, 3. / 4. * (xq * xq - 1), 1. / 4. * (3 * xq * xq - 2 * xq - 1),
          0, 3. / 4. * (1 - xq * xq), 1. / 4. * (3 * xq * xq + 2 * xq - 1)};
}

bool testBending(const Array<Real> & shape_functions, UInt shape_line_index,
                 std::function<Vector<Real>(Real)> reference) {
  Real xq = -1. / std::sqrt(3.);

  // Testing values for bending rotations shapes on quadrature points
  for (auto && N : make_view(shape_functions, 3, 6)) {
    auto Nt = N.transpose();
    Vector<Real> N_bending = Nt(shape_line_index);
    auto bending_reference = reference(xq);
    if (!Math::are_vector_equal(6, N_bending.storage(),
                                bending_reference.storage()))
      return false;
    xq *= -1;
  }

  std::cout.flush();
  return true;
}

int main(int argc, char * argv[]) {
  akantu::initialize(argc, argv);
  // debug::setDebugLevel(dblTest);

  constexpr ElementType type = _bernoulli_beam_2;
  UInt dim = ElementClass<type>::getSpatialDimension();

  Mesh mesh(dim);
  // creating nodes
  Vector<Real> node = {0, 0};
  mesh.getNodes().push_back(node);
  node = {3. / 5., 4. / 5.};
  mesh.getNodes().push_back(node);
  node = {2 * 3. / 5., 0};
  mesh.getNodes().push_back(node);

  mesh.addConnectivityType(type);
  auto & connectivity = mesh.getConnectivity(type);

  // creating elements
  Vector<UInt> elem = {0, 1};
  connectivity.push_back(elem);
  elem = {1, 2};
  connectivity.push_back(elem);
  elem = {0, 2};
  connectivity.push_back(elem);

  using FE = FEEngineTemplate<IntegratorGauss, ShapeStructural, _ek_structural>;
  using ShapeStruct = ShapeStructural<_ek_structural>;

  auto fem = std::make_unique<FE>(mesh, dim, "test_fem");

  fem->initShapeFunctions();

  auto & shape = dynamic_cast<const ShapeStruct &>(fem->getShapeFunctions());

  Array<Real> angles;
  angles.push_back(std::atan(4. / 3.));
  angles.push_back(-std::atan(4. / 3.));
  angles.push_back(0);

  /// Testing the rotation matrices
  for (auto && tuple : zip(make_view(shape.getRotations(type), 3, 3), angles)) {
    auto && rotation = std::get<0>(tuple);
    auto theta = std::get<1>(tuple);
    auto reference = globalToLocalRotation(theta);

    if (!Math::are_vector_equal(9, reference.storage(), rotation.storage()))
      return 1;
  }

  auto & shape_functions = shape.getShapes(type);

  if (!testBending(shape_functions, 0, axialReference))
    return 1;
  // if (!testBending(shape_functions, 1, bendingReference))
  //   return 1;
  // if (!testBending(shape_functions, 2, bendingRotationReference))
  //   return 1;

  std::cout.flush();

  finalize();

  return 0;
}
