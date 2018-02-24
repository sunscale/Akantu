/**
 * @file   test_model_solver_dynamic.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Apr 13 2016
 * @date last modification: Tue Feb 20 2018
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
#include "communicator.hh"
#include "data_accessor.hh"
#include "dof_manager.hh"
#include "dof_manager_default.hh"
#include "element_synchronizer.hh"
#include "mesh.hh"
#include "mesh_accessor.hh"
#include "model_solver.hh"
#include "non_linear_solver.hh"
#include "sparse_matrix.hh"
#include "synchronizer_registry.hh"
/* -------------------------------------------------------------------------- */
#include "test_model_solver_my_model.hh"
/* -------------------------------------------------------------------------- */
#include <fstream>
/* -------------------------------------------------------------------------- */

#ifndef EXPLICIT
#define EXPLICIT true
#endif

using namespace akantu;

static void genMesh(Mesh & mesh, UInt nb_nodes);

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  initialize(argc, argv);

  UInt prank = Communicator::getStaticCommunicator().whoAmI();
  UInt global_nb_nodes = 201;
  UInt max_steps = 200;
  Real time_step = 0.001;
  Mesh mesh(1);
  Real F = -9.81;
  bool _explicit = EXPLICIT;

  if (prank == 0)
    genMesh(mesh, global_nb_nodes);

  mesh.distribute();

  MyModel model(F, mesh, _explicit);

  if (!_explicit) {
    model.getNewSolver("dynamic", _tsst_dynamic, _nls_newton_raphson);
    model.setIntegrationScheme("dynamic", "disp", _ist_trapezoidal_rule_2,
                               IntegrationScheme::_displacement);
  } else {
    model.getNewSolver("dynamic", _tsst_dynamic_lumped, _nls_lumped);
    model.setIntegrationScheme("dynamic", "disp", _ist_central_difference,
                               IntegrationScheme::_acceleration);
  }

  model.setTimeStep(time_step);

  // #if EXPLICIT == true
  //   std::ofstream output("output_dynamic_explicit.csv");
  // #else
  //   std::ofstream output("output_dynamic_implicit.csv");
  // #endif

  if (prank == 0) {
    std::cout << std::scientific;
    std::cout << std::setw(14) << "time"
              << "," << std::setw(14) << "wext"
              << "," << std::setw(14) << "epot"
              << "," << std::setw(14) << "ekin"
              << "," << std::setw(14) << "total" << std::endl;
  }
  Real wext = 0.;

  model.getDOFManager().clearResidual();
  model.assembleResidual();

  Real epot = 0; // model.getPotentialEnergy();
  Real ekin = 0; // model.getKineticEnergy();
  Real einit = ekin + epot;
  Real etot = ekin + epot - wext - einit;
  if (prank == 0) {
    std::cout << std::setw(14) << 0. << "," << std::setw(14) << wext << ","
              << std::setw(14) << epot << "," << std::setw(14) << ekin << ","
              << std::setw(14) << etot << std::endl;
  }

#if EXPLICIT == false
  NonLinearSolver & solver =
      model.getDOFManager().getNonLinearSolver("dynamic");

  solver.set("max_iterations", 20);
#endif

  for (UInt i = 1; i < max_steps + 1; ++i) {
    model.solveStep("dynamic");

#if EXPLICIT == false
// int nb_iter = solver.get("nb_iterations");
// Real error = solver.get("error");
// bool converged = solver.get("converged");
// if (prank == 0)
// std::cerr << error << " " << nb_iter << " -> " << converged << std::endl;
#endif

    epot = model.getPotentialEnergy();
    ekin = model.getKineticEnergy();
    wext += model.getExternalWorkIncrement();
    etot = ekin + epot - wext - einit;
    if (prank == 0) {
      std::cout << std::setw(14) << time_step * i << "," << std::setw(14)
                << wext << "," << std::setw(14) << epot << "," << std::setw(14)
                << ekin << "," << std::setw(14) << etot << std::endl;
    }

    if (std::abs(etot) > 1e-1) {
      AKANTU_ERROR("The total energy of the system is not conserved!");
    }
  }

  // output.close();
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
