/**
 * @file   test_sparse_matrix_profile.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Dec 13 2010
 * @date last modification: Mon Sep 15 2014
 *
 * @brief  test the profile generation of the SparseMatrix class
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include <cstdlib>
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"

#include "sparse_matrix.hh"
#include "dof_synchronizer.hh"

#ifdef AKANTU_USE_SCOTCH
#include "mesh_partition_scotch.hh"
#endif

/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[]) {
  akantu::initialize(argc, argv);

  akantu::UInt spatial_dimension = 2;
  akantu::Mesh mesh(spatial_dimension);
  akantu::MeshIOMSH mesh_io;
  mesh_io.read("triangle.msh", mesh);

  akantu::SparseMatrix sparse_matrix_hand(10, akantu::_symmetric, "hand");

  for(akantu::UInt i = 0; i < 10; ++i) {
    sparse_matrix_hand.addToProfile(i, i);
  }

  sparse_matrix_hand.addToProfile(0,9);

  for(akantu::UInt i = 0; i < 10; ++i) {
    sparse_matrix_hand.addToMatrix(i, i, i*10);
  }
  sparse_matrix_hand.addToMatrix(0,9, 100);

  sparse_matrix_hand.saveProfile("profile_hand.mtx");
  sparse_matrix_hand.saveMatrix("matrix_hand.mtx");


  /* ------------------------------------------------------------------------ */
  akantu::UInt nb_nodes = mesh.getNbNodes();
  akantu::SparseMatrix sparse_matrix(nb_nodes * spatial_dimension, akantu::_symmetric, "mesh");
  akantu::DOFSynchronizer dof_synchronizer(mesh, spatial_dimension);
  dof_synchronizer.initGlobalDOFEquationNumbers();

  sparse_matrix.buildProfile(mesh, dof_synchronizer, spatial_dimension);
  sparse_matrix.saveProfile("profile.mtx");
  /* ------------------------------------------------------------------------ */

#ifdef AKANTU_USE_SCOTCH
  mesh_io.write("triangle_breorder.msh", mesh);

  akantu::MeshPartition * partition = new akantu::MeshPartitionScotch(mesh, spatial_dimension);
  partition->reorder();
  delete partition;

  akantu::DOFSynchronizer dof_synchronizer_re(mesh, spatial_dimension);
  dof_synchronizer_re.initGlobalDOFEquationNumbers();

  sparse_matrix.buildProfile(mesh, dof_synchronizer_re, spatial_dimension);
  sparse_matrix.saveProfile("profile_reorder.mtx");
  mesh_io.write("triangle_reorder.msh", mesh);
#endif

  akantu::finalize();

  return EXIT_SUCCESS;
}
