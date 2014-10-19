/**
 * @file   test_petsc_matrix_profile.cc
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Wed Jul 30 12:34:08 2014
 *
 * @brief  test the profile generation of the PETScMatrix class
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
#include "mesh.hh"
#include "mesh_io.hh"

#include "petsc_matrix.hh"
#include "dof_synchronizer.hh"

#include "mesh_partition_scotch.hh"

using namespace akantu;
int main(int argc, char *argv[]) {

  initialize(argc, argv);
  UInt spatial_dimension = 2;

  StaticCommunicator & comm = akantu::StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  /// read the mesh and partition it
  Mesh mesh(spatial_dimension);
  akantu::MeshPartition * partition = NULL;

  if(prank == 0) {
    /// creation mesh
    mesh.read("triangle.msh");
    partition = new MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
  }

  
  DOFSynchronizer dof_synchronizer(mesh, spatial_dimension);
  UInt nb_global_nodes = mesh.getNbGlobalNodes();

  PETScMatrix petsc_matrix(nb_global_nodes * spatial_dimension, _symmetric);
  
  dof_synchronizer.initGlobalDOFEquationNumbers();

  petsc_matrix.resize(dof_synchronizer);
  petsc_matrix.buildProfile(mesh, dof_synchronizer, spatial_dimension);

  petsc_matrix.saveMatrix("profile.mtx");

}
