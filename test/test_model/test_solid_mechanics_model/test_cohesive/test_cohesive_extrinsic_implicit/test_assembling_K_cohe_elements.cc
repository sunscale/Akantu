/**
 * @file   test_assembling_K_cohe_elements.cc
 *
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 *
 * @date creation: Fri May 15 2015
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Test to check the correct matrix assembling for cohesive elements
 * with degenerated nodes
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "dof_manager.hh"
#include "non_linear_solver.hh"
#include "solid_mechanics_model_cohesive.hh"
#include "sparse_matrix.hh"
/* -------------------------------------------------------------------------- */
#include <fstream>
#include <iostream>
#include <limits>
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char * argv[]) {
  initialize("material.dat", argc, argv);

  debug::setDebugLevel(dblWarning);

  const UInt spatial_dimension = 2;
  Real increment = 0.004;
  bool passed = true;
  Real tol = 1.0e-13;

  Mesh mesh(spatial_dimension);
  mesh.read("quadrangle.msh");

  SolidMechanicsModelCohesive model(mesh);

  /// model initialization
  model.initFull(SolidMechanicsModelCohesiveOptions(_static, true));

  ///  CohesiveElementInserter
  model.getElementInserter().setLimit(_y, -0.001, 0.001);
  model.updateAutomaticInsertion();

  Array<bool> & boundary = model.getBlockedDOFs();
  Array<Real> & position = mesh.getNodes();
  Array<Real> & displacement = model.getDisplacement();
  // SparseMatrix & K_test = model.getStiffnessMatrix();
  Array<Real> K_verified(0, 3, "K_matrix_verified");
  Array<Real> K_test(0, 3, "K_matrix_test");

  /// load the verified stiffness matrix
  Vector<Real> tmp(3);
  UInt nb_lines;

  std::ifstream infile("K_matrix_verified.dat");
  std::string line;
  if (!infile.good())
    AKANTU_ERROR("Cannot open file K_matrix_verified.dat");
  else {
    for (UInt i = 0; i < 2; ++i) {
      getline(infile, line);
      std::stringstream sstr_data(line);
      if (i == 1) {
        sstr_data >> tmp(0) >> tmp(1) >> tmp(2);
        nb_lines = tmp(2);
      }
    }

    for (UInt i = 0; i < nb_lines; ++i) {
      getline(infile, line);
      std::stringstream sstr_data(line);
      sstr_data >> tmp(0) >> tmp(1) >> tmp(2);
      K_verified.push_back(tmp);
    }
  }
  infile.close();

  /// impose boundary conditions
  for (UInt n = 0; n < mesh.getNbNodes(); ++n) {
    if (position(n, 1) < -0.99) {
      boundary(n, 1) = true;
      boundary(n, 0) = true;
    }
    if (position(n, 1) > 0.99 && position(n, 0) < -0.99)
      boundary(n, 1) = true;
  }

  /// solve step
  for (UInt n = 0; n < mesh.getNbNodes(); ++n) {
    if (position(n, 1) > 0.99 && position(n, 0) < -0.99)
      displacement(n, 1) += increment;
  }

  auto & solver = model.getNonLinearSolver();
  solver.set("max_iterations", 10);
  solver.set("threshold", 1e-13);

  model.solveStep();

  model.getDOFManager().getMatrix("K").saveMatrix("K_matrix_test.dat");

  /// load the stiffness matrix to be tested
  std::ifstream infile2("K_matrix_test.dat");
  if (!infile2.good())
    AKANTU_ERROR("Cannot open file K_matrix_test.dat");
  else {
    for (UInt i = 0; i < 2; ++i) {
      getline(infile2, line);
      std::stringstream sstr_data(line);
      if (i == 1) {
        sstr_data >> tmp(0) >> tmp(1) >> tmp(2);
        nb_lines = tmp(2);
      }
    }

    for (UInt i = 0; i < nb_lines; ++i) {
      getline(infile2, line);
      std::stringstream sstr_data(line);
      sstr_data >> tmp(0) >> tmp(1) >> tmp(2);
      K_test.push_back(tmp);
    }
  }
  infile2.close();

  for (UInt i = 0; i < K_verified.size(); ++i) {
    for (UInt j = 0; j < K_test.size(); ++j) {
      if ((K_test(j, 0) == K_verified(i, 0)) &&
          (K_test(j, 1) == K_verified(i, 1))) {
        if (std::abs(K_verified(i, 2)) < tol) {
          if (std::abs(K_test(j, 2)) > tol)
            passed = false;
        } else {
          Real ratio = (std::abs(K_test(j, 2) - K_verified(i, 2))) /
                       (std::abs(K_verified(i, 2)));
          if (ratio > tol)
            passed = false;
        }
      }
    }
  }

  finalize();

  if (passed)
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;
}
