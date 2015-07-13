/**
 * @file   test_interpolate.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 17 2011
 * @date last modification: Thu Jun 05 2014
 *
 * @brief  test of the fem class
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
#include "shape_lagrange.hh"
#include "integrator_gauss.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */
using namespace akantu;

int main(int argc, char *argv[]) {
  akantu::initialize(argc, argv);

  debug::setDebugLevel(dblTest);

  const ElementType type = TYPE;
  UInt dim = ElementClass<type>::getSpatialDimension();

  Mesh my_mesh(dim);

  std::stringstream meshfilename; meshfilename << type << ".msh";
  my_mesh.read(meshfilename.str());

  FEEngine *fem = new FEEngineTemplate<IntegratorGauss,ShapeLagrange>(my_mesh, dim, "my_fem");

  fem->initShapeFunctions();

  std::cout << *fem << std::endl;

  delete fem;
  finalize();

  return 0;
}
