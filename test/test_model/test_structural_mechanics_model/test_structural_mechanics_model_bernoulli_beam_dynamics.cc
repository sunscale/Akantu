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
#include <fstream>
#include <limits>

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "material.hh"
#include "mesh.hh"
#include "mesh_accessor.hh"
#include "mesh_io.hh"
#include "mesh_io_msh_struct.hh"
#include "structural_mechanics_model.hh"

#include <iostream>
using namespace akantu;
/* -------------------------------------------------------------------------- */

#define TYPE _bernoulli_beam_2

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
  Mesh beams(2);
  debug::setDebugLevel(dblWarning);
  const ElementType type = _bernoulli_beam_2;

  /* --------------------------------------------------------------------------
   */
  // Mesh
  UInt nb_element = 9;
  UInt nb_nodes = nb_element + 1;
  Real total_length = 10.;
  Real length = total_length / nb_element;
  Real heigth = 1.;

  MeshAccessor mesh_accessor(beams);
  auto & nodes = mesh_accessor.getNodes();
  nodes.resize(nb_nodes);

  beams.addConnectivityType(type);
  auto & connectivity = mesh_accessor.getConnectivity(type);
  connectivity.resize(nb_element);

  beams.initNormals();
  auto & normals =
      mesh_accessor.getData<Real>("extra_normal", type, _not_ghost, 2);
  normals.resize(nb_element);

  for (UInt i = 0; i < nb_nodes; ++i) {
    nodes(i, 0) = i * length;
    nodes(i, 1) = 0;
  }

  for (UInt i = 0; i < nb_element; ++i) {
    connectivity(i, 0) = i;
    connectivity(i, 1) = i + 1;
  }

  mesh_accessor.makeReady();

  /* --------------------------------------------------------------------------
   */
  // Materials

  StructuralMechanicsModel model(beams);

  StructuralMaterial mat1;
  mat1.E = 120e6;
  mat1.rho = 1000;
  mat1.A = heigth;
  mat1.I = heigth * heigth * heigth / 12;
  model.addMaterial(mat1);

  /* --------------------------------------------------------------------------
   */
  // Forces
  // model.initFull();
  model.initFull(_analysis_method = _implicit_dynamic);

  const auto & position = beams.getNodes();
  auto & forces = model.getExternalForce();
  auto & displacement = model.getDisplacement();
  auto & boundary = model.getBlockedDOFs();

  UInt node_to_print = -1;

  forces.zero();
  displacement.zero();
  //  boundary.zero();
  // model.getElementMaterial(type)(i,0) = 0;
  // model.getElementMaterial(type)(i,0) = 1;
  for (UInt i = 0; i < nb_element; ++i) {
    model.getElementMaterial(type)(i, 0) = 0;
  }

  for (UInt n = 0; n < nb_nodes; ++n) {
    Real x = position(n, _x);
    //    Real y = position(n, 1);

    if (Math::are_float_equal(x, total_length / 2.)) {
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
  // model.computeForcesFromFunction<type>(load, _bft_traction)
  /* --------------------------------------------------------------------------
   */
  // Boundary conditions
  // true ~ displacement is blocked
  boundary(0, 0) = true;
  boundary(0, 1) = true;
  boundary(nb_nodes - 1, 1) = true;
  /* --------------------------------------------------------------------------
   */
  // "Solve"

  Real time = 0;
  // model.assembleStiffnessMatrix();
  // model.assembleMass();
  model.addDumpFieldVector("displacement");
  model.addDumpFieldVector("force");
  model.addDumpFieldVector("velocity");
  model.addDumpFieldVector("acceleration");
  model.dump();
  Real time_step = 1e-4;
  model.setTimeStep(time_step);

  std::cout << "Time  |   Mid-Span Displacement" << std::endl;

  /// time loop
  for (UInt s = 1; time < 0.64; ++s) {

    model.solveStep();

    pos << s << "," << time << "," << displacement(node_to_print, 1) << ","
        << analytical_solution(s * time_step, total_length, mat1.rho, mat1.E,
                               mat1.A, mat1.I, F)
        << std::endl;
    //    pos << s << "," << time << "," << displacement(node_to_print, 1) <<
    //    "," << analytical_solution(s*time_step) << std::endl;

    time += time_step;
    if (s % 100 == 0)
      std::cout << time << "  |   " << displacement(node_to_print, 1)
                << std::endl;
    model.dump();
  }

  model.getDOFManager().getMatrix("K").saveMatrix("K.mtx");
  model.getDOFManager().getMatrix("M").saveMatrix("M.mtx");


  pos.close();

  finalize();

  return EXIT_SUCCESS;
}
