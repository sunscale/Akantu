/**
 * @file   test_material_thermal.cc
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Thu Jan 16 2014
 * @date last modification: Fri Sep 19 2014
 *
 * @brief  test of the class akantu::MaterialThermal
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
#include <iostream>

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"

/* -------------------------------------------------------------------------- */
using namespace akantu;

int main(int argc, char *argv[]) {
  debug::setDebugLevel(dblWarning);
  initialize("material_thermal.dat", argc, argv);

  Math::setTolerance(1.e-13);

  Mesh mesh(2);
  mesh.read("square.msh");

  SolidMechanicsModel model(mesh);
  model.initFull(SolidMechanicsModelOptions(_static));

  mesh.computeBoundingBox();
  const Vector<Real> & min = mesh.getLowerBounds();
  const Vector<Real> & max = mesh.getUpperBounds();

  Array<Real> & pos = mesh.getNodes();
  Array<bool> & boundary = model.getBlockedDOFs();
  Array<Real> & disp = model.getDisplacement();

  for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
    if (Math::are_float_equal(pos(i, 0), min(0))) {
      boundary(i, 0) = true;
    }

    if (Math::are_float_equal(pos(i, 1), min(1))) {
      boundary(i, 1) = true;
    }
  }

  model.setBaseName("test_material_thermal");
  model.addDumpField("displacement");
  model.addDumpField("strain");
  model.addDumpField("stress");
  model.addDumpField("delta_T");

  model.solveStatic<_scm_newton_raphson_tangent_modified, _scc_increment>(1e-10, 2);

  for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
    if (Math::are_float_equal(pos(i, 0), max(0)) && Math::are_float_equal(pos(i, 1), max(1))) {
      if (!Math::are_float_equal(disp(i, 0), 1.0) || !Math::are_float_equal(disp(i, 1), 1.0)) {
	AKANTU_DEBUG_ERROR("Test not passed");
        return EXIT_FAILURE;
      }
    }
  }

  model.dump();

  finalize();

  std::cout << "Test passed" << std::endl;
  return EXIT_SUCCESS;
}



