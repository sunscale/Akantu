/**
 * @file   test_sparse_matrix_product.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 17 2011
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  test the matrix vector product in parallel
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
#include <iostream>
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "dof_synchronizer.hh"
#include "element_synchronizer.hh"
#include "mesh.hh"
#include "mesh_partition_scotch.hh"
#include "sparse_matrix_aij.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  initialize(argc, argv);
  const UInt spatial_dimension = 2;
  const UInt nb_dof = 2;

  const auto & comm = Communicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  Mesh mesh(spatial_dimension);
  mesh.read("bar.msh");
  mesh.distribute();

  UInt nb_nodes = mesh.getNbNodes();
  DOFManagerDefault dof_manager(mesh, "test_dof_manager");

  Array<Real> test_synchronize(nb_nodes, nb_dof, "Test vector");
  dof_manager.registerDOFs("test_synchronize", test_synchronize, _dst_nodal);

  if (prank == 0)
    std::cout << "Creating a SparseMatrix" << std::endl;

  auto & A = dynamic_cast<SparseMatrixAIJ &>(
      dof_manager.getNewMatrix("A", _symmetric));

  Array<Real> dof_vector(nb_nodes, nb_dof, "vector");

  if (prank == 0)
    std::cout << "Filling the matrix" << std::endl;

  for (UInt i = 0; i < nb_nodes * nb_dof; ++i) {
    if (dof_manager.isLocalOrMasterDOF(i))
      A.add(i, i, 2.);
  }

  std::stringstream str;
  str << "Matrix_" << prank << ".mtx";
  A.saveMatrix(str.str());

  for (UInt n = 0; n < nb_nodes; ++n) {
    for (UInt d = 0; d < nb_dof; ++d) {
      dof_vector(n, d) = 1.;
    }
  }

  Array<Real> dof_vector_tmp(dof_vector);

  if (prank == 0)
    std::cout << "Computing x = A * x" << std::endl;
  A.matVecMul(dof_vector, dof_vector_tmp);
  dof_vector.copy(dof_vector_tmp);

  auto & sync =
      dynamic_cast<DOFManagerDefault &>(dof_manager).getSynchronizer();

  if (prank == 0)
    std::cout << "Gathering the results on proc 0" << std::endl;
  if (psize > 1) {
    if (prank == 0) {
      Array<Real> gathered;
      sync.gather(dof_vector, gathered);

      debug::setDebugLevel(dblTest);
      std::cout << gathered << std::endl;
      debug::setDebugLevel(dblWarning);
    } else {
      sync.gather(dof_vector);
    }
  } else {
    debug::setDebugLevel(dblTest);
    std::cout << dof_vector << std::endl;
    debug::setDebugLevel(dblWarning);
  }

  finalize();

  return 0;
}
