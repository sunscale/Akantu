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
#include "aka_common.hh"
#include "dof_synchronizer.hh"
#include "element_synchronizer.hh"
#include "mesh.hh"
#include "mesh_accessor.hh"
#include "mesh_partition_scotch.hh"
#include "sparse_matrix_aij.hh"
#include "sparse_solver_mumps.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;

/* -------------------------------------------------------------------------- */
void genMesh(Mesh & mesh, UInt nb_nodes);
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  initialize(argc, argv);
  const UInt spatial_dimension = 1;
  const UInt nb_global_dof = 100;
  StaticCommunicator & comm =
      akantu::StaticCommunicator::getStaticCommunicator();
  // Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  Mesh mesh(spatial_dimension);

  if (prank == 0)
    genMesh(mesh, nb_global_dof);

  mesh.distribute();

  UInt nb_nodes = mesh.getNbNodes();

  DOFManagerDefault dof_manager(mesh, "test_dof_manager");

  Array<Real> x(nb_nodes);
  dof_manager.registerDOFs("x", x, _dst_nodal);

  Array<UInt> local_equation_number(nb_nodes);
  dof_manager.getEquationsNumbers("x", local_equation_number);

  auto & A = dof_manager.getNewMatrix("A", _symmetric);

  Array<Real> b(nb_nodes);

  for (UInt i = 0; i < nb_nodes; ++i) {
    if (dof_manager.isLocalOrMasterDOF(i)) {
      UInt eqn = local_equation_number(i);
      A.addToMatrix(eqn, eqn, 1. / (eqn+1));
    }
  }

  std::stringstream str;
  str << "Matrix_" << prank << ".mtx";
  A.saveMatrix(str.str());

  for (UInt n = 0; n < nb_nodes; ++n) {
    b(n) = 1.;
  }

  SparseSolverMumps solver(dof_manager, "A");

  solver.solve(x, b);

  auto & sync =
      dynamic_cast<DOFManagerDefault &>(dof_manager).getSynchronizer();

  if (prank == 0) {
    Array<Real> x_gathered(dof_manager.getSystemSize());
    sync.gather(x, x_gathered);

    debug::setDebugLevel(dblTest);
    std::cout << x_gathered << std::endl;
    debug::setDebugLevel(dblWarning);

    UInt d = 1.;
    for(auto x : x_gathered) {
      if ( std::abs(x - d)/d > 1e-15 ) AKANTU_EXCEPTION("Error in the solution");
      ++d;
    }
  } else {
    sync.gather(x);
  }

  finalize();

  return 0;
}

/* -------------------------------------------------------------------------- */
void genMesh(Mesh & mesh, UInt nb_nodes) {
  MeshAccessor mesh_accessor(mesh);
  Array<Real> & nodes = mesh_accessor.getNodes();
  Array<UInt> & conn = mesh_accessor.getConnectivity(_segment_2);

  nodes.resize(nb_nodes);

  for (UInt n = 0; n < nb_nodes; ++n) {
    nodes(n, _x) = n * (1. / (nb_nodes - 1));
  }

  conn.resize(nb_nodes - 1);
  for (UInt n = 0; n < nb_nodes - 1; ++n) {
    conn(n, 0) = n;
    conn(n, 1) = n + 1;
  }
}
