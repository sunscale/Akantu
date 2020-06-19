/**
 * @file   test_fe_engine_precomputation_bernoulli_3.cc
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Jan 24 2018
 *
 * @brief  test of the fem class
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
#include "fe_engine.hh"
#include "integrator_gauss.hh"
#include "shape_structural.hh"
/* -------------------------------------------------------------------------- */
#include <cmath>
#include <functional>
#include <iostream>
/* -------------------------------------------------------------------------- */
using namespace akantu;

/**
 * Reference: p. 285, example 5.7 - A First Course in the Finite Elements Method
 * Logan, 6th Edition, 2016
 * ISBN-13: 978-1-305-63734-4
 */
Matrix<Real> rotationReference() {
  return {{3. / 13, 4. / 13, 12. / 13},
          {-4. / 5, 3. / 5, 0},
          {-36. / 65, -48. / 65, 5. / 13}};
}

int main(int argc, char * argv[]) {
  akantu::initialize(argc, argv);
  // debug::setDebugLevel(dblTest);

  constexpr ElementType type = _bernoulli_beam_3;
  UInt dim = ElementClass<type>::getSpatialDimension();

  Mesh mesh(dim);

  // Pushing nodes
  Vector<Real> node = {0, 0, 0};
  mesh.getNodes().push_back(node);
  node = {3, 4, 12};
  mesh.getNodes().push_back(node);

  // Pushing connectivity
  mesh.addConnectivityType(type);
  auto & connectivity = mesh.getConnectivity(type);
  Vector<UInt> elem = {0, 1};
  connectivity.push_back(elem);

  // Pushing normals
  auto & normals = mesh.registerElementalData<Real>("extra_normal")
                       .alloc(0, dim, type, _not_ghost);
  Vector<Real> normal = {-36. / 65, -48. / 65, 5. / 13};
  normals.push_back(normal);
  normals.push_back(normal);

  using FE = FEEngineTemplate<IntegratorGauss, ShapeStructural, _ek_structural>;
  using ShapeStruct = ShapeStructural<_ek_structural>;

  auto fem = std::make_unique<FE>(mesh, dim, "test_fem");

  fem->initShapeFunctions();

  auto & shape = dynamic_cast<const ShapeStruct &>(fem->getShapeFunctions());

  Matrix<Real> rot_ref = rotationReference();
  Matrix<Real> solution(6, 6);
  solution.block(rot_ref, 0, 0);
  solution.block(rot_ref, 3, 3);

  for (auto && rot : make_view(shape.getRotations(type), 6, 6)) {
    if (!Math::are_vector_equal(6 * 6, solution.storage(), rot.storage()))
      return 1;
  }

  /// TODO check shape functions and shape derivatives

  finalize();

  return 0;
}
