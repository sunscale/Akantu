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
const Real F = 0.5e4;

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  initialize(argc, argv);

  /* ------------------------------------------------------------------------ */
  UInt spatial_dimension = 2;
  // auto method = _static;
  auto method = _implicit_dynamic;
  UInt nb_element = 100;
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
  //nb_element *= 2;

  UInt nb_nodes = nb_element + 1;
  Real total_length = 2.;
  Real length = total_length / nb_element;

  MeshAccessor mesh_accessor(beams);
  auto & nodes = mesh_accessor.getNodes();
  nodes.resize(nb_nodes);

  beams.addConnectivityType(type);
  auto & connectivity = mesh_accessor.getConnectivity(type);
  connectivity.resize(nb_element);

  beams.initNormals();
  auto & normals =
      mesh_accessor.getData<Real>("extra_normal", type, _not_ghost, spatial_dimension);
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

  StructuralMaterial mat1;
  mat1.E = 120e6;
  mat1.rho = 1000;

  mat1.A = 0.003;

  // Real heigth = .3;
  // Real width  = .1;
  mat1.I = 0.000225;
  mat1.Iz = 0.000225;
  mat1.Iy = 0.000025;
  mat1.GJ = 0.0000790216;

  model.addMaterial(mat1);

  /* ------------------------------------------------------------------------ */
  // Forces
  //model.initFull(_analysis_method = _static);
  model.initFull(_analysis_method = method);

  const auto & position = beams.getNodes();
  auto & forces = model.getExternalForce();
  auto & displacement = model.getDisplacement();
  auto & boundary = model.getBlockedDOFs();

  forces.zero();
  displacement.zero();
  //  boundary.zero();
  // model.getElementMaterial(type)(i,0) = 0;
  // model.getElementMaterial(type)(i,0) = 1;
  UInt node_to_print = -1;
  for (UInt i = 0; i < nb_element; ++i) {
    model.getElementMaterial(type)(i, 0) = 0;
  }

  for (UInt n = 0; n < nb_nodes; ++n) {
    Real x = position(n, _x);
    //    Real y = position(n, 1);

    if (Math::are_float_equal(x, total_length/2)) {
      forces(n, _y) = F;
      node_to_print = n;
    }
  }

  std::ofstream pos;
  pos.open("position.csv");
  if (not pos.good()) {
    AKANTU_ERROR("Cannot open file \"position.csv\"");
  }

  pos << "id,time,position,solution" << std::endl;
  /* --------------------------------------------------------------------------
   */
  // Boundary conditions
  // true ~ displacement is blocked
  boundary(0, _x) = true;
  boundary(0, _y) = true;

  boundary(nb_nodes - 1, _y) = true;

  if(spatial_dimension == 3) {
    boundary(0, _z) = true;
    boundary(nb_nodes - 1, _z) = true;
  }

  /* ------------------------------------------------------------------------ */
  // "Solve"
  Real time = 0;
  // model.assembleStiffnessMatrix();
  // model.assembleMass();
  model.addDumpFieldVector("displacement");
  model.addDumpFieldVector("force");
  model.addDumpFieldVector("internal_force");
  if (method == _implicit_dynamic) {
    model.addDumpFieldVector("velocity");
    model.addDumpFieldVector("acceleration");
  }
  model.dump();

  Real time_step = 1e-4;
  model.setTimeStep(time_step);

  model.assembleMatrix("K");
  model.getDOFManager().getMatrix("K").saveMatrix("K.mtx");

  model.assembleMatrix("M");
  model.getDOFManager().getMatrix("M").saveMatrix("M.mtx");

  auto & solver = model.getNonLinearSolver();
  solver.set("max_iterations", 100);
  solver.set("threshold", 1e-8);
  solver.set("convergence_type", SolveConvergenceCriteria::_solution);

  if (method == _static) {
    model.solveStep();
    model.dump();
    return 0;
  }

  std::cout << "Time  |   Mid-Span Displacement" << std::endl;
  /// time loop
  for (UInt s = 1; time < 0.64; ++s) {
    try {
      model.solveStep();
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
      std::cout << std::setw(5) << time << "  |   "
                << std::setw(10) <<displacement(node_to_print, 1) << " |  " << n_iter
                << std::endl;
      model.dump(time, s);
    }

  }

  pos.close();

  return EXIT_SUCCESS;
}
