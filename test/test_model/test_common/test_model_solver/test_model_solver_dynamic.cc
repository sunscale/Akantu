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
#include "element_group.hh"
#include "mesh.hh"
#include "mesh_accessor.hh"
#include "non_linear_solver.hh"
/* -------------------------------------------------------------------------- */
#include "dumpable_inline_impl.hh"
#include "dumper_element_partition.hh"
#include "dumper_iohelper_paraview.hh"
/* -------------------------------------------------------------------------- */
#include "test_model_solver_my_model.hh"
/* -------------------------------------------------------------------------- */
#include <fstream>
/* -------------------------------------------------------------------------- */

#ifndef EXPLICIT
#define EXPLICIT true
#endif

using namespace akantu;

class Sinusoidal : public BC::Dirichlet::DirichletFunctor {
public:
  Sinusoidal(MyModel & model, Real amplitude, Real pulse_width, Real t)
      : model(model), A(amplitude), k(2 * M_PI / pulse_width),
        t(t), v{std::sqrt(model.E / model.rho)} {}

  void operator()(UInt n, Vector<bool> & /*flags*/, Vector<Real> & disp,
                  const Vector<Real> & coord) const {
    auto x = coord(_x);
    model.velocity(n, _x) = k * v * A * sin(k * (x - v * t));
    disp(_x) = A * cos(k * (x - v * t));
  }

private:
  MyModel & model;
  Real A{1.};
  Real k{2 * M_PI};
  Real t{1.};
  Real v{1.};
};

static void genMesh(Mesh & mesh, UInt nb_nodes);

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  initialize(argc, argv);

  UInt prank = Communicator::getStaticCommunicator().whoAmI();
  UInt global_nb_nodes = 201;
  UInt max_steps = 400;
  Real time_step = 0.001;
  Mesh mesh(1);
  Real F = -9.81;
  bool _explicit = EXPLICIT;
  const Real pulse_width = 0.2;
  const Real A = 0.01;

  ID dof_manager_type = "default";
#if defined(DOF_MANAGER_TYPE)
  dof_manager_type = DOF_MANAGER_TYPE;
#endif

  if (prank == 0)
    genMesh(mesh, global_nb_nodes);

  mesh.distribute();

  // mesh.makePeriodic(_x);

  MyModel model(F, mesh, _explicit, dof_manager_type);

  model.forces.zero();
  model.blocked.zero();

  model.applyBC(Sinusoidal(model, A, pulse_width, 0.), "all");
  model.applyBC(BC::Dirichlet::FlagOnly(_x), "border");

  if (!_explicit) {
    model.getNewSolver("dynamic", TimeStepSolverType::_dynamic,
                       NonLinearSolverType::_newton_raphson);
    model.setIntegrationScheme("dynamic", "disp",
                               IntegrationSchemeType::_trapezoidal_rule_2,
                               IntegrationScheme::_displacement);
  } else {
    model.getNewSolver("dynamic", TimeStepSolverType::_dynamic_lumped,
                       NonLinearSolverType::_lumped);
    model.setIntegrationScheme("dynamic", "disp",
                               IntegrationSchemeType::_central_difference,
                               IntegrationScheme::_acceleration);
  }

  model.setTimeStep(time_step);

  if (prank == 0) {
    std::cout << std::scientific;
    std::cout << std::setw(14) << "time"
              << "," << std::setw(14) << "wext"
              << "," << std::setw(14) << "epot"
              << "," << std::setw(14) << "ekin"
              << "," << std::setw(14) << "total"
              << "," << std::setw(14) << "max_disp"
              << "," << std::setw(14) << "min_disp" << std::endl;
  }
  Real wext = 0.;

  model.getDOFManager().zeroResidual();
  model.assembleResidual();

  Real epot = 0; // model.getPotentialEnergy();
  Real ekin = 0; // model.getKineticEnergy();
  Real einit = ekin + epot;
  Real etot = ekin + epot - wext - einit;

  Real max_disp = 0., min_disp = 0.;
  for (auto && disp : model.displacement) {
    max_disp = std::max(max_disp, disp);
    min_disp = std::min(min_disp, disp);
  }

  if (prank == 0) {
    std::cout << std::setw(14) << 0. << "," << std::setw(14) << wext << ","
              << std::setw(14) << epot << "," << std::setw(14) << ekin << ","
              << std::setw(14) << etot << "," << std::setw(14) << max_disp
              << "," << std::setw(14) << min_disp << std::endl;
  }

#if EXPLICIT == false
  NonLinearSolver & solver =
      model.getDOFManager().getNonLinearSolver("dynamic");

  solver.set("max_iterations", 20);
#endif

  auto && dumper = std::make_shared<DumperParaview>("dynamic", "./paraview");
  mesh.registerExternalDumper(dumper, "dynamic", true);
  mesh.addDumpMesh(mesh);

  mesh.addDumpFieldExternalToDumper("dynamic", "displacement",
                                    model.displacement);
  mesh.addDumpFieldExternalToDumper("dynamic", "velocity", model.velocity);
  mesh.addDumpFieldExternalToDumper("dynamic", "forces", model.forces);
  mesh.addDumpFieldExternalToDumper("dynamic", "internal_forces",
                                    model.internal_forces);
  mesh.addDumpFieldExternalToDumper("dynamic", "acceleration",
                                    model.acceleration);

  mesh.dump();

  for (UInt i = 1; i < max_steps + 1; ++i) {
    model.applyBC(Sinusoidal(model, A, pulse_width, time_step * (i - 1)),
                  "border");

    model.solveStep("dynamic");
    mesh.dump();

    epot = model.getPotentialEnergy();
    ekin = model.getKineticEnergy();
    wext += model.getExternalWorkIncrement();
    etot = ekin + epot - wext - einit;

    Real max_disp = 0., min_disp = 0.;
    for (auto && disp : model.displacement) {
      max_disp = std::max(max_disp, disp);
      min_disp = std::min(min_disp, disp);
    }

    if (prank == 0) {
      std::cout << std::setw(14) << time_step * i << "," << std::setw(14)
                << wext << "," << std::setw(14) << epot << "," << std::setw(14)
                << ekin << "," << std::setw(14) << etot << "," << std::setw(14)
                << max_disp << "," << std::setw(14) << min_disp << std::endl;
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

  auto & all = mesh.createNodeGroup("all_nodes");

  for (UInt n = 0; n < nb_nodes; ++n) {
    nodes(n, _x) = n * (1. / (nb_nodes - 1));
    all.add(n);
  }

  mesh.createElementGroupFromNodeGroup("all", "all_nodes");

  conn.resize(nb_nodes - 1);
  for (UInt n = 0; n < nb_nodes - 1; ++n) {
    conn(n, 0) = n;
    conn(n, 1) = n + 1;
  }

  Array<UInt> & conn_points = mesh_accessor.getConnectivity(_point_1);
  conn_points.resize(2);

  conn_points(0, 0) = 0;
  conn_points(1, 0) = nb_nodes - 1;

  auto & border = mesh.createElementGroup("border", 0);
  border.add({_point_1, 0, _not_ghost}, true);
  border.add({_point_1, 1, _not_ghost}, true);

  mesh_accessor.makeReady();
}
