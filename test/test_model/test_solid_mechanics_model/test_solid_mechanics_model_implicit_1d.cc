/**
 * @file   test_solid_mechanics_model_implicit_1d.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Aug 09 2010
 * @date last modification: Sun Oct 19 2014
 *
 * @brief  test of traction in implicit
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
#include <limits>
#include <fstream>

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model.hh"

/* -------------------------------------------------------------------------- */
using namespace akantu;

int main(int argc, char *argv[])
{
  initialize("material.dat", argc, argv);
  UInt spatial_dimension = 1;

  Mesh mesh(spatial_dimension);
  mesh.read("segment1.msh");

  SolidMechanicsModel model(mesh);
  model.initFull(_analysis_method = _static);

  std::cout << model.getMaterial(0) << std::endl;

  /// boundary conditions
  model.getBlockedDOFs()(0,0) = true;
  model.getForce()(1,0) = 1000;

  model.setBaseName("implicit_1d");
  model.addDumpField("displacement");
  model.addDumpField("velocity"    );
  model.addDumpField("acceleration");
  model.addDumpFieldVector("external_force");
  model.addDumpFieldVector("internal_force");
  model.addDumpField("stress"      );
  model.addDumpField("strain"      );

  debug::setDebugLevel(dblInfo);

  model.dump();

  model.solveStep();

  model.dump();

  finalize();

  return EXIT_SUCCESS;
}
