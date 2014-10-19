/**
 * @file   test_purify_mesh.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Jul 30 2012
 * @date last modification: Tue Nov 06 2012
 *
 * @brief  Test the purifyMesh function from MeshUtils
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

#include "mesh.hh"
#include "mesh_utils.hh"
#include "mesh_io.hh"

using namespace akantu;

int main(int argc, char *argv[]) {
  akantu::initialize(argc, argv);
  
  Mesh mesh(2);

  MeshIOMSH mesh_io;
  mesh_io.read("purify_mesh.msh", mesh);
  
  MeshUtils::purifyMesh(mesh);

  mesh_io.write("purify_mesh_after.msh", mesh);

  if(mesh.getNbNodes() != 21)
    AKANTU_DEBUG_ERROR("The purified mesh does not contain the good number of nodes.");

  if(mesh.getNbElement(_quadrangle_8) != 4)
    AKANTU_DEBUG_ERROR("The purified mesh does not contain the good number of element.");

  
  akantu::finalize();
  
  return EXIT_SUCCESS;
}

