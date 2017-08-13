/**
 * @file   test_gradient.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Peter Spijker <peter.spijker@epfl.ch>
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
 * @section DESCRIPTION
 *
 * This code is computing the gradient of a linear field and check that it gives
 * a constant result.  It also compute the gradient the  coordinates of the mesh
 * and check that it gives the identity
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

int main(int argc, char * argv[]) {
  akantu::initialize(argc, argv);
  debug::setDebugLevel(dblTest);

  const ElementType type = TYPE;
  UInt dim = ElementClass<type>::getSpatialDimension();

  Real eps = 1e-12;
  std::cout << "Epsilon : " << eps << std::endl;

  Mesh my_mesh(dim);

  std::stringstream meshfilename;
  meshfilename << type << ".msh";
  my_mesh.read(meshfilename.str());

  FEEngine * fem = new FEEngineTemplate<IntegratorGauss, ShapeLagrange>(
      my_mesh, dim, "my_fem");

  fem->initShapeFunctions();

  Real alpha[2][3] = {{13, 23, 31}, {11, 7, 5}};

  /// create the 2 component field
  const Array<Real> & position = fem->getMesh().getNodes();
  Array<Real> const_val(fem->getMesh().getNbNodes(), 2, "const_val");

  UInt nb_element = my_mesh.getNbElement(type);
  UInt nb_quadrature_points = fem->getNbIntegrationPoints(type) * nb_element;

  Array<Real> grad_on_quad(nb_quadrature_points, 2 * dim, "grad_on_quad");
  for (UInt i = 0; i < const_val.size(); ++i) {
    const_val(i, 0) = 0;
    const_val(i, 1) = 0;

    for (UInt d = 0; d < dim; ++d) {
      const_val(i, 0) += alpha[0][d] * position(i, d);
      const_val(i, 1) += alpha[1][d] * position(i, d);
    }
  }

  /// compute the gradient
  fem->gradientOnIntegrationPoints(const_val, grad_on_quad, 2, type);

  // std::cout << "Linear array on nodes : " << const_val << std::endl;
  // std::cout << "Gradient on quad : " << grad_on_quad << std::endl;

  /// check the results
  Array<Real>::matrix_iterator it = grad_on_quad.begin(2, dim);
  Array<Real>::matrix_iterator it_end = grad_on_quad.end(2, dim);
  for (; it != it_end; ++it) {
    for (UInt d = 0; d < dim; ++d) {
      Matrix<Real> & grad = *it;
      if (!(std::abs(grad(0, d) - alpha[0][d]) < eps) ||
          !(std::abs(grad(1, d) - alpha[1][d]) < eps)) {
        std::cout << "Error gradient is not correct " << (*it)(0, d) << " "
                  << alpha[0][d] << " (" << std::abs((*it)(0, d) - alpha[0][d])
                  << ")"
                  << " - " << (*it)(1, d) << " " << alpha[1][d] << " ("
                  << std::abs((*it)(1, d) - alpha[1][d]) << ")"
                  << " - " << d << std::endl;
        std::cout << *it << std::endl;
        exit(EXIT_FAILURE);
      }
    }
  }

  // compute gradient of coordinates
  Array<Real> grad_coord_on_quad(nb_quadrature_points, dim * dim,
                                 "grad_coord_on_quad");
  fem->gradientOnIntegrationPoints(my_mesh.getNodes(), grad_coord_on_quad,
                                   my_mesh.getSpatialDimension(), type);

  // std::cout << "Node positions : " << my_mesh.getNodes() << std::endl;
  // std::cout << "Gradient of nodes : " << grad_coord_on_quad << std::endl;

  Array<Real>::matrix_iterator itp = grad_coord_on_quad.begin(dim, dim);
  Array<Real>::matrix_iterator itp_end = grad_coord_on_quad.end(dim, dim);

  for (; itp != itp_end; ++itp) {
    for (UInt i = 0; i < dim; ++i) {
      for (UInt j = 0; j < dim; ++j) {
        if (!(std::abs((*itp)(i, j) - (i == j)) < eps)) {
          std::cout << *itp << std::endl;
          exit(EXIT_FAILURE);
        }
      }
    }
  }

  delete fem;
  finalize();

  return EXIT_SUCCESS;
}
