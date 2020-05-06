/**
 * @file   test_sparse_matrix_profile_parallel.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sun Sep 12 2010
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  test the sparse matrix class in parallel
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

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "communicator.hh"
#include "mesh.hh"
#include "mesh_io_msh.hh"
#include "mesh_partition_scotch.hh"
#include "solver_mumps.hh"
#include "sparse_matrix.hh"

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/* Main                                                                       */
/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  akantu::initialize(argc, argv);

  int dim = 2;
  //#ifdef AKANTU_USE_IOHELPER
  // akantu::ElementType type = akantu::_triangle_6;
  //#endif //AKANTU_USE_IOHELPER
  akantu::Mesh mesh(dim);

  //  akantu::debug::setDebugLevel(akantu::dblDump);

  akantu::StaticCommunicator * comm =
      akantu::Communicator::getStaticCommunicator();
  akantu::Int psize = comm->getNbProc();
  akantu::Int prank = comm->whoAmI();

  akantu::UInt n = 0;

  /* ------------------------------------------------------------------------ */
  /* Parallel initialization                                                  */
  /* ------------------------------------------------------------------------ */
  akantu::Communicator * communicator;
  if (prank == 0) {
    akantu::MeshIOMSH mesh_io;
    mesh_io.read("triangle.msh", mesh);
    akantu::MeshPartition * partition =
        new akantu::MeshPartitionScotch(mesh, dim);

    //    partition->reorder();
    mesh_io.write("triangle_reorder.msh", mesh);

    n = mesh.getNbNodes();

    partition->partitionate(psize);
    communicator =
        akantu::Communicator::createCommunicatorDistributeMesh(mesh, partition);
    delete partition;
  } else {
    communicator =
        akantu::Communicator::createCommunicatorDistributeMesh(mesh, NULL);
  }

  std::cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA "
            << mesh.getNbGlobalNodes() << std::endl;

  akantu::SparseMatrix sparse_matrix(mesh, akantu::_symmetric, 2, "mesh");
  sparse_matrix.buildProfile();

  akantu::Solver * solver = new akantu::SolverMumps(sparse_matrix);

  if (prank == 0) {
    for (akantu::UInt i = 0; i < n; ++i) {
      solver->getRHS().storage()[i] = 1.;
    }
  }

  akantu::debug::setDebugLevel(akantu::dblDump);
  solver->initialize();

  std::stringstream sstr;
  sstr << "profile_" << prank << ".mtx";
  sparse_matrix.saveProfile(sstr.str());

  akantu::finalize();

  return EXIT_SUCCESS;
}
