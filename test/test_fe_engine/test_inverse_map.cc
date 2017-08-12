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
#include "aka_common.hh"
#include "aka_iterators.hh"
#include "fe_engine.hh"
#include "integrator_gauss.hh"
#include "shape_lagrange.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {
  akantu::initialize(argc, argv);

  const auto type = TYPE;
  auto dim = ElementClass<type>::getSpatialDimension();

  Mesh my_mesh(dim);
  std::stringstream meshfilename;
  meshfilename << type << ".msh";

  my_mesh.read(meshfilename.str());

  const auto & lower = my_mesh.getLowerBounds();
  const auto & upper = my_mesh.getUpperBounds();
  UInt nb_elements = my_mesh.getNbElement(type);

  auto fem = std::make_unique<FEEngineTemplate<IntegratorGauss, ShapeLagrange>>(
      my_mesh, dim, "my_fem");

  fem->initShapeFunctions();
  Matrix<Real> quad = GaussIntegrationElement<type>::getQuadraturePoints();

  /// get the quadrature points coordinates
  Array<Real> coord_on_quad(quad.cols() * nb_elements,
                            my_mesh.getSpatialDimension(), "coord_on_quad");

  fem->interpolateOnIntegrationPoints(my_mesh.getNodes(), coord_on_quad,
                                      my_mesh.getSpatialDimension(), type);

  Vector<Real> natural_coords(dim);

  /// loop over the quadrature points
  auto it = coord_on_quad.begin_reinterpret(dim, quad.cols(), nb_elements);
  auto end = coord_on_quad.end_reinterpret(dim, quad.cols(), nb_elements);

  auto length = (upper - lower).norm<L_inf>();

  auto eps = 5e-12;
  UInt el = 0;
  for (; it != end; ++it, ++el) {
    const auto & quads_coords = *it;
    for (auto q : arange(quad.cols())) {
      Vector<Real> quad_coord = quads_coords(q);
      Vector<Real> ref_quad_coord = quad(q);
      fem->inverseMap(quad_coord, el, type, natural_coords);

      auto dis_normalized = ref_quad_coord.distance(natural_coords) / length;
      if (not(dis_normalized < eps)) {
        std::cout << "Real coordinates inversion test failed : "
                  << std::scientific << natural_coords << " - "
                  << ref_quad_coord << " / " << length << " = " << dis_normalized << std::endl;
        return 1;
      }

      // std::cout << "Real coordinates inversion : " << std::scientific
      //           << natural_coords << " = " << ref_quad_coord << " ("
      //           << dis_normalized << ")" << std::endl;
    }
  }

  std::cout << "inverse completed over " << nb_elements << " elements"
            << std::endl;

  finalize();

  return 0;
}
