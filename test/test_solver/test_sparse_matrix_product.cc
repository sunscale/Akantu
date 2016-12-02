/**
 * @file   test_sparse_matrix_product.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 17 2011
 * @date last modification: Sun Oct 19 2014
 *
 * @brief  test the matrix vector product in parallel
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
#include <iostream>
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_partition_scotch.hh"
#include "element_synchronizer.hh"
#include "sparse_matrix.hh"
#include "dof_synchronizer.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
  initialize(argc, argv);
  const UInt spatial_dimension = 2;
  const UInt nb_dof = 3;

  StaticCommunicator & comm = akantu::StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  Mesh mesh(spatial_dimension);

  /* ------------------------------------------------------------------------ */
  /* Parallel initialization                                                  */
  /* ------------------------------------------------------------------------ */
  MeshPartition * partition;

  if(prank == 0) {
    mesh.read("bar.msh");

    std::cout << "Partitioning mesh..." << std::endl;
    partition = new akantu::MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
    ElementSynchronizer::createDistributedSynchronizerMesh(mesh, partition);
    delete partition;
  } else {
    ElementSynchronizer::createDistributedSynchronizerMesh(mesh, NULL);
  }

  UInt nb_nodes = mesh.getNbNodes();
  UInt nb_global_node = mesh.getNbGlobalNodes();

  // Array<Int>  equation_number(nb_nodes, nb_dof);
  // for (UInt n = 0; n < nb_nodes; ++n) {
  //   UInt real_n = mesh.getNodeGlobalId(n);
  //   bool is_local_node = !(mesh.isPureGhostNode(n));
  //   for (UInt d = 0; d < nb_dof; ++d) {
  //     UInt global_eq_num = (is_local_node ? real_n : nb_global_node + real_n) * nb_dof + d;
  //     equation_number(n, d) = global_eq_num;
  //   }
  // }

  if (prank == 0) std::cout << "Creating a SparseMatrix" << std::endl;
  SparseMatrix sparse_matrix(nb_global_node * nb_dof, _symmetric, "matrix");


  DOFSynchronizer dof_synchronizer(mesh, nb_dof);
  dof_synchronizer.initGlobalDOFEquationNumbers();
  sparse_matrix.buildProfile(mesh, dof_synchronizer, nb_dof);

  Array<Real> dof_vector(nb_nodes, nb_dof, "vector");

  if (prank == 0) std::cout << "Filling the matrix" << std::endl;
  UInt nz = sparse_matrix.getNbNonZero();
  const Array<Int> & irn = sparse_matrix.getIRN();
  const Array<Int> & jcn = sparse_matrix.getJCN();
  for (UInt i = 0; i < nz; ++i) {
    sparse_matrix.addToMatrix(irn(i) - 1, jcn(i) - 1, 1.);
  }

  std::stringstream str; str << "Matrix_" << prank << ".mtx";
  sparse_matrix.saveMatrix(str.str());

  for (UInt n = 0; n < nb_nodes; ++n) {
    for (UInt d = 0; d < nb_dof; ++d) {
      dof_vector(n,d) = 1.;
    }
  }

  if (prank == 0) std::cout << "Computing x = A * x" << std::endl;
  dof_vector *= sparse_matrix;

  if (prank == 0) std::cout << "Lumping the matrix" << std::endl;
  Array<Real> lumped(0,nb_dof);
  sparse_matrix.lump(lumped);

  if (prank == 0) std::cout << "Gathering the results on proc 0" << std::endl;
  if(psize > 1) {
    const_cast<DOFSynchronizer &>(sparse_matrix.getDOFSynchronizer()).initScatterGatherCommunicationScheme();
    if(prank == 0) {
      Array<Real> gathered(0, nb_dof);
      Array<Real> lump_gathered(0, nb_dof);
      sparse_matrix.getDOFSynchronizer().gather(dof_vector, 0, &gathered);
      sparse_matrix.getDOFSynchronizer().gather(lumped, 0, &lump_gathered);
      debug::setDebugLevel(dblTest);
      std::cout << gathered << std::endl;
      std::cout << lump_gathered << std::endl;
      debug::setDebugLevel(dblWarning);
    } else {
      sparse_matrix.getDOFSynchronizer().gather(dof_vector, 0);
      sparse_matrix.getDOFSynchronizer().gather(lumped, 0);
    }
  } else {
    debug::setDebugLevel(dblTest);
    std::cout << dof_vector << std::endl;
    std::cout << lumped << std::endl;
    debug::setDebugLevel(dblWarning);
  }

  finalize();

  return 0;
}
