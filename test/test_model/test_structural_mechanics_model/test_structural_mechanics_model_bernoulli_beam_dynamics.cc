/**
 * @file   test_structural_mechanics_model_bernoulli_beam_dynamics.cc
 *
 * @author Sébastien Hartmann <sebastien.hartmann@epfl.ch>
 *
 * @date creation: Mon Jul 07 2014
 * @date last modification: Wed Feb 03 2016
 *
 * @brief  Test for _bernouilli_beam in dynamic
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "mesh_accessor.hh"
#include "non_linear_solver_newton_raphson.hh"
#include "structural_mechanics_model.hh"
/* -------------------------------------------------------------------------- */
#include <fstream>
#include <iostream>
#include <limits>
/* -------------------------------------------------------------------------- */

using namespace akantu;

/* -------------------------------------------------------------------------- */
static Real analytical_solution(Real time, Real L, Real rho, Real E,
                                __attribute__((unused)) Real A, Real I,
                                Real F) {
  Real omega = M_PI * M_PI / L / L * sqrt(E * I / rho);
  Real sum = 0.;
  UInt i = 5;
  for (UInt n = 1; n <= i; n += 2) {
    sum += (1. - cos(n * n * omega * time)) / pow(n, 4);
  }

  return 2. * F * pow(L, 3) / pow(M_PI, 4) / E / I * sum;
}

// load
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  initialize(argc, argv);

  /* ------------------------------------------------------------------------ */
  UInt spatial_dimension = 3;
  Real total_length = 20.;

  //auto method = _static;
  auto method = _implicit_dynamic;
  UInt nb_element = 10;
  /* ------------------------------------------------------------------------ */

  Mesh beams(spatial_dimension);
  debug::setDebugLevel(dblWarning);

  ElementType type = _bernoulli_beam_2;
  if (spatial_dimension == 3) {
    type = _bernoulli_beam_3;
  }
  /* ------------------------------------------------------------------------ */
  // Mesh
  //  UInt nb_element = 8;
  // nb_element *= 2;

  UInt nb_nodes = nb_element + 1;
  Real length = total_length / nb_element;

  MeshAccessor mesh_accessor(beams);
  auto & nodes = mesh_accessor.getNodes();
  nodes.resize(nb_nodes);

  beams.addConnectivityType(type);
  auto & connectivity = mesh_accessor.getConnectivity(type);
  connectivity.resize(nb_element);

  beams.initNormals();
  auto & normals = mesh_accessor.getData<Real>("extra_normal", type, _not_ghost,
                                               spatial_dimension);
  normals.resize(nb_element);

  for (UInt i = 0; i < nb_nodes; ++i) {
    nodes(i, _x) = i * length;
    nodes(i, _y) = 0;
    if (spatial_dimension == 3) {
      nodes(i, _z) = 0;
    }
  }

  for (UInt i = 0; i < nb_element; ++i) {
    connectivity(i, 0) = i;
    connectivity(i, 1) = i + 1;

    if (spatial_dimension == 3) {
      normals(i, _x) = 0;
      normals(i, _y) = 0;
      normals(i, _z) = 1;
    }
  }

  mesh_accessor.makeReady();

  /* ------------------------------------------------------------------------ */
  // Materials
  StructuralMechanicsModel model(beams);

  const Real F = 1;

  StructuralMaterial mat1;
  mat1.E = 1e9;
  mat1.rho = 1;

  mat1.A = 1;

  // Real heigth = .3;
  // Real width  = .1;
  mat1.I = 1;
  mat1.Iz = 1;
  mat1.Iy = 1;
  mat1.GJ = 1;

  model.addMaterial(mat1);

  /* ------------------------------------------------------------------------ */
  // Forces
  // model.initFull(_analysis_method = _static);
  model.initFull(_analysis_method = method);

  auto & forces = model.getExternalForce();
  auto & displacement = model.getDisplacement();
  auto & boundary = model.getBlockedDOFs();

  forces.zero();
  displacement.zero();
  UInt node_to_print = -1;
  for (UInt i = 0; i < nb_element; ++i) {
    model.getElementMaterial(type)(i, 0) = 0;
  }

  // for (UInt n = 0; n < nb_nodes; ++n) {
  //   Real x = position(n, _x);
  //   if (Math::are_float_equal(x, total_length/2)) {
  //     forces(n, _y) = F;
  //     node_to_print = n;
  //   }
  // }

  std::ofstream pos;
  pos.open("position.csv");
  if (not pos.good()) {
    AKANTU_ERROR("Cannot open file \"position.csv\"");
  }

  pos << "id,time,position,solution" << std::endl;
  /* --------------------------------------------------------------------------
   */
  // Boundary conditions
  for (UInt d = 0; d < spatial_dimension; ++d) {
    boundary(0, d) = true;
  }
  for (UInt d = 0; d < spatial_dimension - 1; ++d) {
    boundary(nb_nodes - 1, d + 1) = true;
  }

  node_to_print = nb_nodes / 2 + 1;
  forces(node_to_print - 1, _y) = F;

  /* ------------------------------------------------------------------------ */
  // "Solve"
  Real time = 0;
  // model.assembleStiffnessMatrix();
  // model.assembleMass();
  model.addDumpFieldVector("displacement");
  model.addDumpFieldVector("rotation");
  model.addDumpFieldVector("force");
  model.addDumpFieldVector("momentum");
  model.addDumpFieldVector("internal_force");
  model.addDumpFieldVector("internal_momentum");
  model.addDumpField("stress");
  if (method == _implicit_dynamic) {
    model.addDumpFieldVector("velocity");
    model.addDumpFieldVector("acceleration");
  }
  model.dump();

  Real time_step = 1e-5;
  model.setTimeStep(time_step);

  model.assembleMatrix("K");
  model.getDOFManager().getMatrix("K").saveMatrix("K.mtx");

  model.assembleMatrix("M");
  model.getDOFManager().getMatrix("M").saveMatrix("M.mtx");

  if (method == _static) {
    model.solveStep();
    model.dump();

    model.assembleResidual();

    const auto & stress = model.getStress(type);
    const auto & internal_force = model.getInternalForce();
    std::cout << "θ_a = " << displacement(0, spatial_dimension) << "\n"
              << "θ_b = " << displacement(nb_nodes - 1, spatial_dimension) << "\n"
              << "R_a = " << internal_force(0, _y) << "\n"
              << "R_b = " << internal_force(nb_nodes - 1, _y) << "\n"
              << "w = " << displacement(node_to_print, _y) << "\n"
              << "M = " << stress(stress.size() / 2, 1) << "\n";

    return 0;
  }

  auto & solver = model.getNonLinearSolver();
  solver.set("max_iterations", 100);
  solver.set("threshold", 1e-8);
  solver.set("convergence_type", SolveConvergenceCriteria::_solution);

  std::cout << "Time  |   Mid-Span Displacement" << std::endl;
  /// time loop
  for (UInt s = 1; time < 0.01606; ++s) {
    try {
      model.solveStep();
      model.getDOFManager().getMatrix("J").saveMatrix("J.mtx");
    } catch (debug::NLSNotConvergedException & e) {
      std::cerr << "Did not converge after " << e.niter << " error: " << e.error
                << " > " << e.threshold << std::endl;
      return 1;
    }
    pos << s << "," << time << "," << displacement(node_to_print, 1) << ","
        << analytical_solution(s * time_step, total_length, mat1.rho, mat1.E,
                               mat1.A, mat1.I, F)
        << std::endl;
    //    pos << s << "," << time << "," << displacement(node_to_print, 1) <<
    //    "," << analytical_solution(s*time_step) << std::endl;

    Int n_iter = solver.get("nb_iterations");

    time += time_step;
    if (s % 100 == 0) {
      std::cout << std::setw(5) << time << "  |   " << std::setw(10)
                << displacement(node_to_print, 1) << " |  " << n_iter
                << std::endl;
      model.dump(time, s);
    }
  }

  pos.close();

  return EXIT_SUCCESS;
}
