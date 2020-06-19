/**
 * @file   test_integrate.cc
 *
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 * @date creation: Fri Oct 30 2015
 * @date last modification: Fri Oct 30 2015
 *
 * @brief  test the integration of IGFEM elements
 *
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
#include "integrator_gauss_igfem.hh"
#include "shape_igfem.hh"
///#include "dumper_paraview.hh"
/* -------------------------------------------------------------------------- */
#include <cstdlib>
#include <fstream>
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;

bool integrate(const ElementType type);
void generateIGFEMMesh(const ElementType type, Mesh & mesh,
                       const std::string & filename);

int main(int argc, char * argv[]) {
  akantu::initialize(argc, argv);
  debug::setDebugLevel(dblTest);

  bool test_passed = false;

  /// integrate test of _igfem_triangle_4
  AKANTU_DEBUG_INFO("integrating _igfem_triangle_4...");
  test_passed = integrate(_igfem_triangle_4);
  if (!test_passed) {
    finalize();
    return EXIT_FAILURE;
  }

  /// integrate test of _igfem_triangle_5
  AKANTU_DEBUG_INFO("integrating _igfem_triangle_5...");
  test_passed = integrate(_igfem_triangle_5);
  if (!test_passed) {
    finalize();
    return EXIT_FAILURE;
  }

  finalize();
  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
bool integrate(const ElementType type) {

  bool correct_result = true;
  UInt dim = 2;
  std::stringstream mesh_info;
  mesh_info << "mesh_info" << type << ".txt";

  Mesh my_mesh(dim);
  generateIGFEMMesh(type, my_mesh, mesh_info.str());
  // DumperParaview dumper_igfem("mesh_igfem");
  // dumper_igfem.registerMesh(my_mesh, dim, _not_ghost, _ek_igfem);
  // dumper_igfem.dump();
  Real eps = 3e-13;
  std::cout << "Epsilon : " << eps << std::endl;

  FEEngine * fem =
      new FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_igfem>(
          my_mesh, dim, "my_fem");
  fem->initShapeFunctions();

  UInt nb_element = my_mesh.getNbElement(type);
  UInt nb_quadrature_points = fem->getNbIntegrationPoints(type) * nb_element;

  /// impose a constant field on the quadrature points with x = 1 and y = 2
  Array<Real> val_on_quad(nb_quadrature_points, dim, "val_on_quad");
  Array<Real>::vector_iterator val_it = val_on_quad.begin(dim);
  Vector<Real> values(dim);
  values(0) = 1.;
  values(1) = 2;
  for (UInt i = 0; i < val_on_quad.getSize(); ++i, ++val_it)
    *val_it = values;

  // integrate function on elements
  akantu::Array<akantu::Real> int_val_on_elem(nb_element, dim,
                                              "int_val_on_elem");
  fem->integrate(val_on_quad, int_val_on_elem, dim, type);

  // get global integration value
  Real value[2] = {0, 0};
  std::cout << "Val on quads : " << val_on_quad << std::endl;
  std::cout << "Integral on elements : " << int_val_on_elem << std::endl;

  for (UInt i = 0; i < fem->getMesh().getNbElement(type); ++i) {
    value[0] += int_val_on_elem.storage()[2 * i];
    value[1] += int_val_on_elem.storage()[2 * i + 1];
  }

  std::cout << "integral on the mesh of 1 is " << value[0] << " and of 2 is "
            << value[1] << std::endl;

  delete fem;

  if (!(std::abs(value[0] - 1.) < eps && std::abs(value[1] - 2.) < eps)) {
    std::cout << "|1 - " << value[0] << "| = " << std::abs(value[0] - 1.)
              << std::endl
              << "|2 - " << value[1] << "| = " << std::abs(value[1] - 2.)
              << std::endl;
    correct_result = false;
  }

  return correct_result;
}
