/**
 * @file   test_interpolate.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
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
#include "fe_engine.hh"
#include "shape_lagrange.hh"
#include "integrator_gauss.hh"
/* -------------------------------------------------------------------------- */
#include <cstdlib>
#include <iostream>
/* -------------------------------------------------------------------------- */
using namespace akantu;

int main(int argc, char *argv[]) {
  akantu::initialize(argc, argv);

  debug::setDebugLevel(dblTest);
  const ElementType type = TYPE;
  UInt dim = ElementClass<type>::getSpatialDimension();

  Real eps = 3e-13;
  std::cout << "Epsilon : " << eps << std::endl;

  Mesh my_mesh(dim);

  std::stringstream meshfilename; meshfilename << type << ".msh";
  my_mesh.read(meshfilename.str());

  FEEngine *fem = new FEEngineTemplate<IntegratorGauss,ShapeLagrange>(my_mesh, dim, "my_fem");

  fem->initShapeFunctions();

  Array<Real> const_val(fem->getMesh().getNbNodes(), 2, "const_val");

  UInt nb_element = my_mesh.getNbElement(type);
  UInt nb_quadrature_points = fem->getNbIntegrationPoints(type) * nb_element;

  Array<Real> val_on_quad(nb_quadrature_points, 2, "val_on_quad");

  for (UInt i = 0; i < const_val.size(); ++i) {
    const_val.storage()[i * 2 + 0] = 1.;
    const_val.storage()[i * 2 + 1] = 2.;
  }

  fem->interpolateOnIntegrationPoints(const_val, val_on_quad, 2, type);

  std::cout << "Interpolation of array : " << const_val << std::endl;
  std::cout << "Gives on quads : " << val_on_quad << std::endl;

  // interpolate coordinates
  Array<Real> coord_on_quad(nb_quadrature_points, my_mesh.getSpatialDimension(), "coord_on_quad");

  fem->interpolateOnIntegrationPoints(my_mesh.getNodes(),
				     coord_on_quad,
				     my_mesh.getSpatialDimension(),
				     type);
  std::cout << "Interpolations of node coordinates : " << my_mesh.getNodes() << std::endl;
  std::cout << "Gives : " << coord_on_quad << std::endl;

  delete fem;
  finalize();

  return EXIT_SUCCESS;
}
