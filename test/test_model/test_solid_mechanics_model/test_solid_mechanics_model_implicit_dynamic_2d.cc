/**
 * @file   test_solid_mechanics_model_implicit_dynamic_2d.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed May 11 2011
 * @date last modification: Sun Oct 19 2014
 *
 * @brief  test of the dynamic implicit code
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

/* -------------------------------------------------------------------------- */
#include <fstream>
#include <limits>
/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"
#include "static_communicator.hh"

/* -------------------------------------------------------------------------- */

using namespace akantu;

/* -------------------------------------------------------------------------- */
#define bar_length 10.
#define bar_height 1.
#define bar_depth 1.
const ElementType TYPE = _triangle_6;
// const ElementType TYPE = _tetrahedron_10;

const UInt spatial_dimension = 2;
Real time_step = 1e-4;

const Real F = 0.5e4;
const Real L = bar_length;
const Real I = bar_depth * bar_height * bar_height * bar_height / 12.;
const Real E = 12e7;
const Real rho = 1000;
const Real m = rho * bar_height * bar_depth;

static Real w(UInt n) {
  return n * n * M_PI * M_PI / (L * L) * sqrt(E * I / m);
}

static Real analytical_solution(Real time) {
  return 2 * F * L * L * L / (pow(M_PI, 4) * E * I) *
         ((1. - cos(w(1) * time)) + (1. - cos(w(3) * time)) / 81. +
          (1. - cos(w(5) * time)) / 625.);
}

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  debug::setDebugLevel(dblError);
  initialize("material_implicit_dynamic.dat", argc, argv);

  Mesh mesh(spatial_dimension);

  const auto & comm = StaticCommunicator::getStaticCommunicator();
  Int prank = comm.whoAmI();

  if (prank == 0) {
    mesh.read("beam_2d_quad.msh");
  }

  mesh.distribute();

  SolidMechanicsModel model(mesh);

  /// model initialization
  model.initFull(SolidMechanicsModelOptions(_implicit_dynamic));
  Material & mat = model.getMaterial(0);
  mat.setParam("E", E);
  mat.setParam("rho", rho);

  // boundary conditions
  const Array<Real> & position = mesh.getNodes();
  Array<bool> & boundary = model.getBlockedDOFs();
  Array<Real> & force = model.getForce();
  Array<Real> & displacment = model.getDisplacement();

  // initial conditions
  UInt node_to_print = UInt(-1);
  bool print_node = false;

  Array<UInt> node_to_displace;
  for (UInt n = 0; n < mesh.getNbNodes(); ++n) {
    Real x = position(n, 0);
    Real y = position(n, 1);
    Real z = 0;
    if (spatial_dimension == 3)
      z = position(n, 2);

    if (Math::are_float_equal(x, 0.) &&
        Math::are_float_equal(y, bar_height / 2.)) {
      boundary(n, 0) = true;
      boundary(n, 1) = true;
      if (spatial_dimension == 3 && Math::are_float_equal(z, bar_depth / 2.))
        boundary(n, 2) = true;
    }

    if (Math::are_float_equal(x, bar_length) &&
        Math::are_float_equal(y, bar_height / 2.)) {
      boundary(n, 1) = true;
      if (spatial_dimension == 3 && Math::are_float_equal(z, bar_depth / 2.))
        boundary(n, 2) = true;
    }

    if (Math::are_float_equal(x, bar_length / 2.) &&
        Math::are_float_equal(y, bar_height / 2.)) {

      if (spatial_dimension < 3 || (spatial_dimension == 3 &&
                                    Math::are_float_equal(z, bar_depth / 2.))) {
        force(n, 1) = F;
        if (mesh.isLocalOrMasterNode(n)) {
          print_node = true;
          node_to_print = n;
          std::cout << "I, proc " << prank + 1 << " handle the print of node "
                    << n << "(" << x << ", " << y << ", " << z << ")"
                    << std::endl;
        }
      }
    }
  }

  model.setTimeStep(time_step);

  std::stringstream out;
  out << "position-" << TYPE << "_" << std::scientific << time_step << ".csv";

  model.setBaseName("implicit_dynamic");
  model.addDumpField("displacement");
  model.addDumpField("velocity");
  model.addDumpField("acceleration");
  model.addDumpField("external_force");
  model.addDumpField("internal_force");
  model.addDumpField("stress");
  model.addDumpField("strain");

  std::ofstream pos;
  if (print_node) {
    pos.open(out.str().c_str());
    if (!pos.good()) {
      std::cerr << "Cannot open file " << out.str() << std::endl;
      exit(EXIT_FAILURE);
    }
    pos << "id,time,position,solution" << std::endl;
  }

  Real time = 0;

    /// time loop
  for (UInt s = 1; time < 0.0062; ++s) {
    model.solveStep();
    if (prank == 0)
      std::cout << "passing step " << s << " " << s * time_step << "s\r"
                << std::flush;

    if (print_node)
      pos << s << "," << s * time_step << "," << displacment(node_to_print, 1)
          << "," << analytical_solution(s * time_step) << std::endl;

    // if (s % 10 == 0) {
    //   model.computeStresses();
    // }

    time += time_step;
  }

  if (print_node)
    pos.close();

  finalize();

  return EXIT_SUCCESS;
}
