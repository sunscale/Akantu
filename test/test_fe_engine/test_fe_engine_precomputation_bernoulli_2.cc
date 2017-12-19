/**
 * @file   test_fe_engine_precomputation.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Jun 14 2010
 * @date last modification: Mon Jul 13 2015
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
#include "fe_engine.hh"
#include "shape_structural.hh"
#include "integrator_gauss.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */
using namespace akantu;

int main(int argc, char *argv[]) {
  akantu::initialize(argc, argv);

  // debug::setDebugLevel(dblTest);

  constexpr ElementType type = _bernoulli_beam_2;
  UInt dim = ElementClass<type>::getSpatialDimension();

  Mesh mesh(dim);
  Vector<Real> node = {0, 0};
  mesh.getNodes().push_back(node);
  node = {1, 1};
  mesh.getNodes().push_back(node);

  mesh.addConnectivityType(type);
  auto & connectivity = mesh.getConnectivity(type);
  Vector<UInt> elem = {0, 1};
  connectivity.push_back(elem);

  auto fem =
    std::make_unique<FEEngineTemplate<IntegratorGauss, ShapeStructural, _ek_structural>>(
          mesh, dim, "test_fem");

  fem->initShapeFunctions();

  std::cout << *fem << std::endl;

  finalize();

  return 0;
}
