/**
 * @file   test_cohesive_intrinsic_impl.cc
 *
 * @author Seyedeh Mohadeseh Taheri Mousavi <mohadeseh.taherimousavi@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Mon Jul 09 2012
 * @date last modification: Fri Jul 07 2017
 *
 * @brief  Test for cohesive elements
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
#include "non_linear_solver.hh"
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */
#include <fstream>
#include <iostream>
#include <limits>
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {
  initialize("material.dat", argc, argv);

  debug::setDebugLevel(dblError);

  const UInt spatial_dimension = 2;
  const ElementType type = _triangle_6;

  Mesh mesh(spatial_dimension);
  mesh.read("implicit.msh");

  CohesiveElementInserter inserter(mesh);
  inserter.setLimit(_y, 0.9, 1.1);
  inserter.insertIntrinsicElements();

  //  mesh.write("implicit_cohesive.msh");

  SolidMechanicsModelCohesive model(mesh);

  /// model initialization
  model.initFull(SolidMechanicsModelCohesiveOptions(_static));

  /// boundary conditions
  Array<bool> & boundary = model.getBlockedDOFs();
  UInt nb_nodes = mesh.getNbNodes();
  Array<Real> & position = mesh.getNodes();
  Array<Real> & displacement = model.getDisplacement();
  const ElementType type_facet = mesh.getFacetType(type);

  for (UInt n = 0; n < nb_nodes; ++n) {

    if (std::abs(position(n, 1)) < Math::getTolerance()) {
      boundary(n, 1) = true;
      displacement(n, 1) = 0.0;
    }

    if ((std::abs(position(n, 0)) < Math::getTolerance()) &&
        (position(n, 1) < 1.1)) {
      boundary(n, 0) = true;
      displacement(n, 0) = 0.0;
    }
    if ((std::abs(position(n, 0) - 1) < Math::getTolerance()) &&
        (std::abs(position(n, 1) - 1) < Math::getTolerance())) {
      boundary(n, 0) = true;
      displacement(n, 0) = 0.0;
    }

    if (std::abs(position(n, 1) - 2) < Math::getTolerance()) {
      boundary(n, 1) = true;
    }
  }

  model.setBaseName("intrinsic_impl");
  model.addDumpField("displacement");
  // model.addDumpField("mass"        );
  model.addDumpField("velocity");
  model.addDumpField("acceleration");
  model.addDumpField("force");
  model.addDumpField("residual");
  //  model.addDumpField("damage"      );
  model.addDumpField("stress");
  model.addDumpField("strain");
  model.dump();

  const MaterialCohesive & mat_coh =
      dynamic_cast<const MaterialCohesive &>(model.getMaterial(1));

  ElementType type_cohesive = FEEngine::getCohesiveElementType(type_facet);

  const Array<Real> & opening = mat_coh.getOpening(type_cohesive);
  // const Array<Real> & traction = mat_coh.getTraction(type_cohesive);

  model.assembleInternalForces();
  const Array<Real> & residual = model.getInternalForce();

  UInt max_step = 1000;
  Real increment = 3. / max_step;
  Real error_tol = 10e-6;

  std::ofstream fout;
  fout.open("output");

  auto & solver = model.getNonLinearSolver();
  solver.set("max_iterations", 100);
  solver.set("threshold", 1e-5);
  solver.set("convergence_type", SolveConvergenceCriteria::_residual);

  /// Main loop
  for (UInt nstep = 0; nstep < max_step; ++nstep) {
    for (UInt n = 0; n < nb_nodes; ++n) {
      if (std::abs(position(n, 1) - 2) < Math::getTolerance()) {
        displacement(n, 1) += increment;
      }
    }

    model.solveStep();

    //    model.dump();

    Real resid = 0;
    for (UInt n = 0; n < nb_nodes; ++n) {
      if (std::abs(position(n, 1) - 2.) / 2. < Math::getTolerance()) {
        resid += residual(n, 1);
      }
    }

    Real analytical = exp(1) * std::abs(opening(0, 1)) *
                      exp(-std::abs(opening(0, 1)) / 0.5) / 0.5;

    // the residual force is comparing with the theoretical value of the
    // cohesive law
    error_tol = std::abs((std::abs(resid) - analytical) / analytical);

    fout << nstep << " " << -resid << " " << analytical << " " << error_tol
         << std::endl;

    if (error_tol > 1e-3) {
      std::cout << "Relative error: " << error_tol << std::endl;
      std::cout << "Test failed!" << std::endl;
      return EXIT_FAILURE;
    }
  }

  model.dump();

  fout.close();

  finalize();

  std::cout << "Test passed!" << std::endl;
  return EXIT_SUCCESS;
}
