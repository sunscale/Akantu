/**
 * @file   test_embedded_interface_model.cc
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date   Thu 19 Mar 2015
 *
 * @brief  test of the class EmbeddedInterfaceModel
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#include "aka_common.hh"
#include "embedded_interface_model.hh"

using namespace akantu;

int main(int argc, char * argv[]) {
  debug::setDebugLevel(dblWarning);
  initialize("embedded_element.dat", argc, argv);

  const UInt dim = 2;

  Mesh mesh(dim);
  mesh.read("triangle.msh");

  EmbeddedInterfaceModel model(mesh, dim);

  const Real height = 0.5;

  model.initFull(SolidMechanicsModelOptions(_static));

  if (model.getInterfaceMesh().getNbElement(_segment_2) != 1)
    return EXIT_FAILURE;

  if (model.getInterfaceMesh().getSpatialDimension() != 2)
    return EXIT_FAILURE;

  model.assembleStiffnessMatrix();

  SparseMatrix & K = model.getStiffnessMatrix();

  Math::setTolerance(1e-8);

  // Testing the assembled stiffness matrix
  if (!Math::are_float_equal(K(0, 0), 1. - height) ||
      !Math::are_float_equal(K(0, 2), height - 1.) ||
      !Math::are_float_equal(K(2, 0), height - 1.) ||
      !Math::are_float_equal(K(2, 2), 1. - height))
    return EXIT_FAILURE;

  model.updateResidual();

  return EXIT_SUCCESS;
}
