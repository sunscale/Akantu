/**
 * @file   test_model_solver.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Apr 13 2016
 * @date last modification: Thu Feb 01 2018
 *
 * @brief  Test default dof manager
 *
 * @section LICENSE
 *
 * Copyright (©) 2016-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_random_generator.hh"
#include "dof_manager.hh"
#include "dof_synchronizer.hh"
#include "mesh.hh"
#include "mesh_accessor.hh"
#include "model_solver.hh"
#include "non_linear_solver.hh"
#include "sparse_matrix.hh"

/* -------------------------------------------------------------------------- */
#include "test_model_solver_my_model.hh"
/* -------------------------------------------------------------------------- */
#include <cmath>
/* -------------------------------------------------------------------------- */

using namespace akantu;

static void genMesh(Mesh & mesh, UInt nb_nodes);
static void printResults(MyModel & model, UInt nb_nodes);

Real F = -10;

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  initialize(argc, argv);

  UInt prank = Communicator::getStaticCommunicator().whoAmI();

  std::cout << std::setprecision(7);

  ID dof_manager_type = "default";
#if defined(DOF_MANAGER_TYPE)
  dof_manager_type = DOF_MANAGER_TYPE;
#endif

  UInt global_nb_nodes = 100;
  Mesh mesh(1);

  RandomGenerator<UInt>::seed(1);

  if (prank == 0) {
    genMesh(mesh, global_nb_nodes);
  }

  // std::cout << prank << RandGenerator<Real>::seed() << std::endl;
  mesh.distribute();

  MyModel model(F, mesh, false, dof_manager_type);

  model.getNewSolver("static", TimeStepSolverType::_static,
                     NonLinearSolverType::_newton_raphson);
  model.setIntegrationScheme("static", "disp",
                             IntegrationSchemeType::_pseudo_time);

  NonLinearSolver & solver = model.getDOFManager().getNonLinearSolver("static");
  solver.set("max_iterations", 2);

  model.solveStep();

  printResults(model, global_nb_nodes);

  finalize();
  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
void genMesh(Mesh & mesh, UInt nb_nodes) {
  MeshAccessor mesh_accessor(mesh);
  Array<Real> & nodes = mesh_accessor.getNodes();
  Array<UInt> & conn = mesh_accessor.getConnectivity(_segment_2);

  nodes.resize(nb_nodes);
  mesh_accessor.setNbGlobalNodes(nb_nodes);

  for (UInt n = 0; n < nb_nodes; ++n) {
    nodes(n, _x) = n * (1. / (nb_nodes - 1));
  }

  conn.resize(nb_nodes - 1);
  for (UInt n = 0; n < nb_nodes - 1; ++n) {
    conn(n, 0) = n;
    conn(n, 1) = n + 1;
  }
  mesh_accessor.makeReady();
}

/* -------------------------------------------------------------------------- */
void printResults(MyModel & model, UInt /*nb_nodes*/) {
  // if (model.mesh.isDistributed()) {
  //   UInt prank = model.mesh.getCommunicator().whoAmI();
  //   auto & sync = dynamic_cast<DOFManagerDefault &>(model.getDOFManager())
  //       .getSynchronizer();

  //   if (prank == 0) {
  //     Array<Real> global_displacement(nb_nodes);
  //     Array<Real> global_forces(nb_nodes);
  //     Array<bool> global_blocked(nb_nodes);

  //     sync.gather(model.forces, global_forces);
  //     sync.gather(model.displacement, global_displacement);
  //     sync.gather(model.blocked, global_blocked);

  //     auto force_it = global_forces.begin();
  //     auto disp_it = global_displacement.begin();
  //     auto blocked_it = global_blocked.begin();

  //     std::cout << "node"
  //               << ", " << std::setw(8) << "disp"
  //               << ", " << std::setw(8) << "force"
  //               << ", " << std::setw(8) << "blocked" << std::endl;

  //     UInt node = 0;
  //     for (; disp_it != global_displacement.end();
  //          ++disp_it, ++force_it, ++blocked_it, ++node) {
  //       std::cout << node << ", " << std::setw(8) << *disp_it << ", "
  //                 << std::setw(8) << *force_it << ", " << std::setw(8)
  //                 << *blocked_it << std::endl;

  //       std::cout << std::flush;
  //     }
  //   } else {
  //     sync.gather(model.forces);
  //     sync.gather(model.displacement);
  //     sync.gather(model.blocked);
  //   }
  // } else {
  auto force_it = model.forces.begin();
  auto disp_it = model.displacement.begin();
  auto blocked_it = model.blocked.begin();

  std::cout << "node"
            << ", " << std::setw(8) << "disp"
            << ", " << std::setw(8) << "force"
            << ", " << std::setw(8) << "blocked" << std::endl;

  UInt node = 0;
  for (; disp_it != model.displacement.end();
       ++disp_it, ++force_it, ++blocked_it, ++node) {
    std::cout << node << ", " << std::setw(8) << *disp_it << ", "
              << std::setw(8) << *force_it << ", " << std::setw(8)
              << *blocked_it << std::endl;

    std::cout << std::flush;
  }
  //  }
}
