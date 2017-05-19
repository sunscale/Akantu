/**
 * @file   test_sparse_matrix_profile.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Dec 13 2010
 * @date last modification: Sun Oct 19 2014
 *
 * @brief  test the profile generation of the SparseMatrix class
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
#include "dof_manager_default.hh"
#include "mesh.hh"
#include "sparse_matrix.hh"
/* -------------------------------------------------------------------------- */
#include <cstdlib>
/* -------------------------------------------------------------------------- */

using namespace akantu;

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  initialize(argc, argv);

  UInt spatial_dimension = 2;
  Mesh mesh(spatial_dimension);
  mesh.read("triangle.msh");

  UInt nb_nodes = mesh.getNbNodes();

  DOFManagerDefault dof_manager(mesh, "test_dof_manager");
  Array<Real> test_synchronize(nb_nodes, spatial_dimension, "Test vector");
  dof_manager.registerDOFs("test_synchronize", test_synchronize, _dst_nodal);

  auto & A = dof_manager.getNewMatrix("A", _symmetric);

  for (UInt i = 0; i < 10; ++i) {
    A.addToProfile(i, i);
  }

  A.addToProfile(0, 9);

  for (UInt i = 0; i < 10; ++i) {
    A.addToMatrix(i, i, i * 10);
  }
  A.addToMatrix(0, 9, 100);

  A.saveProfile("profile_hand.mtx");
  A.saveMatrix("matrix_hand.mtx");

  /* ------------------------------------------------------------------------ */
  finalize();

  return EXIT_SUCCESS;
}
