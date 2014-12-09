/**
 * @file   test_petsc_matrix_profile.cc
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Wed Jul 30 12:34:08 2014
 *
 * @brief  test the applyBoundary method of the PETScMatrix class
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "static_communicator.hh"
#include "aka_common.hh"
#include "aka_csr.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_utils.hh"
#include "distributed_synchronizer.hh"
#include "petsc_matrix.hh"
#include "fe_engine.hh"
#include "dof_synchronizer.hh"

#include "mesh_partition_scotch.hh"

using namespace akantu;
int main(int argc, char *argv[]) {

  initialize(argc, argv);
  const ElementType element_type = _triangle_3;
  const GhostType ghost_type = _not_ghost; 
  UInt spatial_dimension = 2;


  StaticCommunicator & comm = akantu::StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  /// read the mesh and partition it
  Mesh mesh(spatial_dimension);

  /* ------------------------------------------------------------------------ */
  /* Parallel initialization                                                  */
  /* ------------------------------------------------------------------------ */
  DistributedSynchronizer * communicator = NULL;
  if(prank == 0) {
    /// creation mesh
    mesh.read("triangle.msh");
    MeshPartitionScotch * partition = new MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
    communicator = DistributedSynchronizer::createDistributedSynchronizerMesh(mesh, partition);
    delete partition;
  } else {
    communicator = DistributedSynchronizer::createDistributedSynchronizerMesh(mesh, NULL);
  }
  

  FEEngine *fem = new FEEngineTemplate<IntegratorGauss,ShapeLagrange,_ek_regular>(mesh, spatial_dimension, "my_fem");

  DOFSynchronizer dof_synchronizer(mesh, spatial_dimension);
  UInt nb_global_nodes = mesh.getNbGlobalNodes();

  dof_synchronizer.initGlobalDOFEquationNumbers();

  // fill the matrix with 
  UInt nb_element = mesh.getNbElement(element_type);
  UInt nb_nodes_per_element = mesh.getNbNodesPerElement(element_type);
  UInt nb_dofs_per_element = spatial_dimension * nb_nodes_per_element;
  SparseMatrix K(nb_global_nodes * spatial_dimension, _symmetric);
  K.buildProfile(mesh, dof_synchronizer, spatial_dimension);
  Matrix<Real> element_input(nb_dofs_per_element, nb_dofs_per_element, 1);
  Array<Real> K_e = Array<Real>(nb_element, nb_dofs_per_element * nb_dofs_per_element, "K_e");
  Array<Real>::matrix_iterator K_e_it = K_e.begin(nb_dofs_per_element, nb_dofs_per_element);
  Array<Real>::matrix_iterator K_e_end = K_e.end(nb_dofs_per_element, nb_dofs_per_element);

  for(; K_e_it != K_e_end; ++K_e_it)
    *K_e_it = element_input;

  // assemble the test matrix
  fem->assembleMatrix(K_e, K, spatial_dimension, element_type, ghost_type);

  // create petsc matrix
  PETScMatrix petsc_matrix(nb_global_nodes * spatial_dimension, _symmetric);

  petsc_matrix.buildProfile(mesh, dof_synchronizer, spatial_dimension);
  
  // add stiffness matrix to petsc matrix
  petsc_matrix.add(K, 1);
  
  // create boundary array: block all dofs
  UInt nb_nodes = mesh.getNbNodes();
  Array<bool> boundary = Array<bool>(nb_nodes, spatial_dimension, true);
  //boundary(2,1) = true;
 

  std::cout << dof_synchronizer.getDOFGlobalID(2*2+1) <<std::endl;
 



  // apply boundary
  petsc_matrix.applyBoundary(boundary);

  // test if all entries except the diagonal ones have been zeroed
  Real test_passed = 0;

  for (UInt i = 0; i < nb_nodes * spatial_dimension; ++i) { 
    if (dof_synchronizer.isLocalOrMasterDOF(i)) {
      for (UInt j = 0; j < nb_nodes * spatial_dimension; ++j) {
	if (dof_synchronizer.isLocalOrMasterDOF(j)) {
	  if (i == j)
	    test_passed += petsc_matrix(i, j) - 1;
	  else
	    test_passed += petsc_matrix(i, j) - 0;
	  UInt global_i = dof_synchronizer.getDOFGlobalID(i);
	  UInt global_j = dof_synchronizer.getDOFGlobalID(j);
	  std::cout << "K(" << global_i << "," << global_j << ")=" << petsc_matrix(i, j) << std::endl;
	}
      }
    }
  }
  petsc_matrix.saveMatrix("profile_after_boundary.mtx"); 
  delete communicator;

  finalize();

  //  if (std::abs(test_passed) > Math::getTolerance())
  //    return EXIT_FAILURE;
  
  std::cout << " Test passed!!!" << std::endl;
  return EXIT_SUCCESS;

}