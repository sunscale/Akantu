/**
 * @file   test_petsc_matrix_diagonal.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 * @date creation: Mon Oct 13 2014
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  test the connectivity is correctly represented in the PETScMatrix
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_csr.hh"
#include "communicator.hh"
#include "dof_synchronizer.hh"
#include "dumper_paraview.hh"
#include "element_synchronizer.hh"
#include "fe_engine.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_partition_scotch.hh"
#include "mesh_utils.hh"
#include "petsc_matrix.hh"

using namespace akantu;
int main(int argc, char * argv[]) {

  initialize(argc, argv);
  const ElementType element_type = _triangle_3;
  const GhostType ghost_type = _not_ghost;
  UInt spatial_dimension = 2;

  const auto & comm = akantu::Communicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  /// read the mesh and partition it
  Mesh mesh(spatial_dimension);

  /* ------------------------------------------------------------------------ */
  /* Parallel initialization                                                  */
  /* ------------------------------------------------------------------------ */
  ElementSynchronizer * communicator = NULL;
  if (prank == 0) {
    /// creation mesh
    mesh.read("triangle.msh");
    MeshPartitionScotch * partition =
        new MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
    communicator =
        ElementSynchronizer::createDistributedSynchronizerMesh(mesh, partition);
    delete partition;
  } else {
    communicator =
        ElementSynchronizer::createDistributedSynchronizerMesh(mesh, NULL);
  }

  // DumperParaview mesh_dumper("mesh_dumper");
  // mesh_dumper.registerMesh(mesh, spatial_dimension, _not_ghost);
  // mesh_dumper.dump();

  /// initialize the FEEngine and the dof_synchronizer
  FEEngine * fem =
      new FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_regular>(
          mesh, spatial_dimension, "my_fem");

  DOFSynchronizer dof_synchronizer(mesh, spatial_dimension);
  UInt nb_global_nodes = mesh.getNbGlobalNodes();

  dof_synchronizer.initGlobalDOFEquationNumbers();

  // construct an Akantu sparse matrix, build the profile and fill the matrix
  // for the given mesh
  UInt nb_element = mesh.getNbElement(element_type);
  UInt nb_nodes_per_element = mesh.getNbNodesPerElement(element_type);
  UInt nb_dofs_per_element = spatial_dimension * nb_nodes_per_element;
  SparseMatrix K_akantu(nb_global_nodes * spatial_dimension, _unsymmetric);
  K_akantu.buildProfile(mesh, dof_synchronizer, spatial_dimension);
  /// use as elemental matrices a matrix with values equal to 1 every where
  Matrix<Real> element_input(nb_dofs_per_element, nb_dofs_per_element, 1.);
  Array<Real> K_e =
      Array<Real>(nb_element, nb_dofs_per_element * nb_dofs_per_element, "K_e");
  Array<Real>::matrix_iterator K_e_it =
      K_e.begin(nb_dofs_per_element, nb_dofs_per_element);
  Array<Real>::matrix_iterator K_e_end =
      K_e.end(nb_dofs_per_element, nb_dofs_per_element);

  for (; K_e_it != K_e_end; ++K_e_it)
    *K_e_it = element_input;

  // assemble the test matrix
  fem->assembleMatrix(K_e, K_akantu, spatial_dimension, element_type,
                      ghost_type);

  /// construct a PETSc matrix
  PETScMatrix K_petsc(nb_global_nodes * spatial_dimension, _unsymmetric);
  /// build the profile of the PETSc matrix for the mesh of this example
  K_petsc.buildProfile(mesh, dof_synchronizer, spatial_dimension);
  /// add an Akantu sparse matrix to a PETSc sparse matrix
  K_petsc.add(K_akantu, 1);

  /// check to how many elements each node is connected
  CSR<Element> node_to_elem;

  MeshUtils::buildNode2Elements(mesh, node_to_elem, spatial_dimension);

  /// test the diagonal of the PETSc matrix: the diagonal entries
  /// of the PETSc matrix correspond to the number of elements
  /// connected to the node of the dof. Note: for an Akantu matrix this is only
  /// true for the serial case
  Real error = 0.;
  /// loop over all diagonal values of the matrix
  for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
    for (UInt j = 0; j < spatial_dimension; ++j) {
      UInt dof = i * spatial_dimension + j;
      /// for PETSc matrix only DOFs on the processor and be accessed
      if (dof_synchronizer.isLocalOrMasterDOF(dof)) {
        UInt global_dof = dof_synchronizer.getDOFGlobalID(dof);
        std::cout << "Number of elements connected: "
                  << node_to_elem.getNbCols(i) << std::endl;
        std::cout << "K_petsc(" << global_dof << "," << global_dof
                  << ")=" << K_petsc(dof, dof) << std::endl;
        error += std::abs(K_petsc(dof, dof) - node_to_elem.getNbCols(i));
      }
    }
  }

  if (error > Math::getTolerance()) {
    std::cout << "error in the stiffness matrix!!!" << std::endl;
    finalize();
    return EXIT_FAILURE;
  }

  delete communicator;

  finalize();

  return EXIT_SUCCESS;
}
