/**
 * @file   matrix_assembling_cohesive_elements.cc
 *
 * @author Mauro Corrado
 *
 * @date creation: Thu April 30 2015
 *
 * @brief  Test to check the correct matrix assembling for cohesive elements
 *         with degenerated nodes
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include <limits>
#include <fstream>
#include <iostream>

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model_cohesive.hh"
#include "material_cohesive.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[]) {
  initialize("material.dat", argc, argv);

  debug::setDebugLevel(dblWarning);

  const UInt spatial_dimension = 2;

  Mesh mesh(spatial_dimension);
  mesh.read("quadrangle.msh");

  SolidMechanicsModelCohesive model(mesh);

  /// model initialization
  model.initFull(SolidMechanicsModelCohesiveOptions(_static, true));

  //  CohesiveElementInserter inserter(mesh);
  model.limitInsertion(_y, -0.001, 0.001);
  model.updateAutomaticInsertion();

  /// boundary conditions
  Array<bool> & boundary = model.getBlockedDOFs();
  Array<Real> & position = mesh.getNodes();
  Array<Real> & displacement = model.getDisplacement();

  for (UInt n = 0; n < mesh.getNbNodes(); ++n) {
    if (position(n,1) < -0.99){
      boundary(n,1) = true;
      boundary(n,0) = true;
    }
    if (position(n,1) > 0.99 && position(n,0) < -0.99)
      boundary(n,1) = true;
  }

  Real increment = 0.004;
  Real error;

  for (UInt n = 0; n < mesh.getNbNodes(); ++n) {
    if (position(n,1) > 0.99 && position(n,0) < -0.99)
      displacement(n,1) += increment;
  }

  model.solveStepCohesive<_scm_newton_raphson_tangent, _scc_increment>(1e-13, error, 10);

  /// save the matrix
  model.getStiffnessMatrix().saveMatrix("K_matrix_test");

  finalize();

  return EXIT_SUCCESS;
}
