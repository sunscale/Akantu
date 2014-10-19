/**
 * @file   test_pbc_tweak.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Wed Feb 09 2011
 * @date last modification: Fri Jun 13 2014
 *
 * @brief  test of internal facet extraction
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_utils.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

int main(int argc, char *argv[])
{
  int dim = 3;

  initialize("material.dat", argc, argv);
  debug::setDebugLevel(akantu::dblInfo);

  Mesh mesh(dim);
  mesh.read("cube.msh");

  SolidMechanicsModel model(mesh);
  /* -------------------------------------------------------------------------- */
  model.initFull();
  /* -------------------------------------------------------------------------- */
  model.setPBC(1,1,1);
  model.initPBC();
  model.assembleMassLumped();
  /* -------------------------------------------------------------------------- */

  model.setBaseName("test-pbc-tweak");
  model.addDumpField("mass");
  model.dump();

  finalize();

  return EXIT_SUCCESS;
}
