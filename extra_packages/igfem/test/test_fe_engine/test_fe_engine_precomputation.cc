/**
 * @file   test_fe_engine_percompuation.cc
 *
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 * @date creation: Fri Oct 30 2015
 * @date last modification: Fri Oct 30 2015
 *
 * @brief  test the fe-engine precomputations for IGFEM elements
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
/* -------------------------------------------------------------------------- */
#include <fstream>
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;

void precompute(const ElementType type);
void generateIGFEMMesh(const ElementType type, Mesh & mesh,
                       const std::string & filename);

int main(int argc, char * argv[]) {
  akantu::initialize(argc, argv);
  debug::setDebugLevel(dblTest);

  /// precompuation for _igfem_triangle_4
  AKANTU_DEBUG_INFO("precomputation for _igfem_triangle_4...");
  precompute(_igfem_triangle_4);

  /// precompuation for _igfem_triangle_5
  AKANTU_DEBUG_INFO("precomputation for _igfem_triangle_5...");
  precompute(_igfem_triangle_5);
  finalize();
  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
void precompute(const ElementType type) {

  UInt dim = 2;
  std::stringstream mesh_info;
  mesh_info << "mesh_info" << type << ".txt";

  Mesh my_mesh(dim);
  generateIGFEMMesh(type, my_mesh, mesh_info.str());

  FEEngine * fem =
      new FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_igfem>(
          my_mesh, dim, "my_fem");

  fem->initShapeFunctions();

  std::cout << *fem << std::endl;

  delete fem;
}
