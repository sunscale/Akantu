/**
 * @file   test_petsc_matrix_profile_parallel.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 * @date creation: Mon Oct 13 2014
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  test the profile generation of the PETScMatrix class in parallel
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
#include <cstdlib>
#include <fstream>
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_csr.hh"
#include "communicator.hh"
#include "dof_synchronizer.hh"
#include "element_synchronizer.hh"
#include "fe_engine.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_utils.hh"
#include "sparse_matrix_petsc.hh"
/// #include "dumper_paraview.hh"
#include "mesh_partition_scotch.hh"

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
    mesh.read("square.msh");
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

  // dump mesh in paraview
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

  /// save the profile
  K_petsc.saveMatrix("profile_parallel.txt");
  /// print the matrix to screen
  if (prank == 0) {
    std::ifstream profile;
    profile.open("profile_parallel.txt");
    std::string current_line;
    while (getline(profile, current_line))
      std::cout << current_line << std::endl;
    profile.close();
  }

  delete communicator;

  finalize();

  return EXIT_SUCCESS;
}
