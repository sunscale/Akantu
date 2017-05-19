/**
 * @file   test_dof_manager_default.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Feb 24 12:28:44 2016
 *
 * @brief  Test default dof manager
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
#include "dof_manager.hh"
#include "dof_synchronizer.hh"
#include "mesh.hh"
#include "mesh_accessor.hh"
#include "model_solver.hh"
#include "sparse_matrix.hh"
/* -------------------------------------------------------------------------- */
#include "test_model_solver_my_model.hh"
/* -------------------------------------------------------------------------- */
#include <fstream>
#include <cmath>
/* -------------------------------------------------------------------------- */

using namespace akantu;

static void genMesh(Mesh & mesh, UInt nb_nodes);
static void printResults(MyModel & model, UInt nb_nodes);

Real F = -10;

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  initialize(argc, argv);

  UInt prank = StaticCommunicator::getStaticCommunicator().whoAmI();

  std::cout << std::setprecision(7);

  UInt global_nb_nodes = 11;
  Mesh mesh(1);

  if (prank == 0)
    genMesh(mesh, global_nb_nodes);

  mesh.distribute();

  MyModel model(F, mesh, false);

  model.getNewSolver("static", _tsst_static, _nls_newton_raphson);
  model.setIntegrationScheme("static", "disp", _ist_pseudo_time);

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

  for (UInt n = 0; n < nb_nodes; ++n) {
    nodes(n, _x) = n * (1. / (nb_nodes - 1));
  }

  conn.resize(nb_nodes - 1);
  for (UInt n = 0; n < nb_nodes - 1; ++n) {
    conn(n, 0) = n;
    conn(n, 1) = n + 1;
  }
}

/* -------------------------------------------------------------------------- */
void printResults(MyModel & model, UInt nb_nodes) {
  UInt prank = StaticCommunicator::getStaticCommunicator().whoAmI();
  auto & sync = dynamic_cast<DOFManagerDefault &>(model.getDOFManager())
    .getSynchronizer();

  if (prank == 0) {
    Array<Real> global_displacement(nb_nodes);
    Array<Real> global_forces(nb_nodes);
    Array<bool> global_blocked(nb_nodes);

    sync.gather(model.forces, global_forces);

    auto force_it = global_forces.begin();
    bool a = std::abs(global_forces(nb_nodes - 1) - F) > 1e-10;
    while(a) {}

    sync.gather(model.displacement, global_displacement);
    sync.gather(model.blocked, global_blocked);

    auto disp_it = global_displacement.begin();
    auto blocked_it = global_blocked.begin();

    std::cout << "node"
              << ", " << std::setw(8) << "disp"
              << ", " << std::setw(8) << "force"
              << ", " << std::setw(8) << "blocked" << std::endl;

    UInt node = 0;
    for (; disp_it != global_displacement.end();
         ++disp_it, ++force_it, ++blocked_it, ++node) {
      std::cout << node << ", " << std::setw(8) << *disp_it << ", "
                << std::setw(8) << *force_it << ", " << std::setw(8)
                << *blocked_it << std::endl;

      std::cout << std::flush;
    }
  } else {
    sync.gather(model.forces);
    sync.gather(model.displacement);
    sync.gather(model.blocked);
  }
}
