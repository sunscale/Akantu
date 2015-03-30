/**
 * @file   reinforced_concrete
 *
 * @author Lucas Frerot
 *
 * @date   lun. 09 29 16:03:10 2014
 *
 * @brief  Reinforced concrete simulation
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#include <iostream>

#include "aka_common.hh"
#include "embedded_interface_model.hh"

using namespace akantu;

int main (int argc, char * argv[]) {
  debug::setDebugLevel(dblWarning);
  initialize("material.dat", argc, argv);

  UInt dim = 2;
  Math::setTolerance(1e-7);
  
  // Mesh here is a 1x1 patch
  Mesh mesh(dim);
  mesh.read("embedded_mesh.msh");

  EmbeddedInterfaceModel model(mesh, dim);
  model.initFull(SolidMechanicsModelOptions(_static));

  Array<Real> & nodes  = mesh.getNodes();
  Array<Real> & forces = model.getForce();
  Array<bool> & bound  = model.getBlockedDOFs();

  forces(2, 0) = -250;
  forces(5, 0) = -500;
  forces(8, 0) = -250;

  for (UInt i = 0 ; i < mesh.getNbNodes() ; i++) {
    if (Math::are_float_equal(nodes(i, 0), 0.))
      bound(i, 0) = true;
    if (Math::are_float_equal(nodes(i, 1), 0.))
      bound(i, 1) = true;
  }

  // Assemble the global stiffness matrix
  model.assembleStiffnessMatrix();
  model.updateResidual();

  model.getStiffnessMatrix().saveMatrix("matrix_test");

  model.solveStatic<_scm_newton_raphson_tangent_not_computed, _scc_residual>(1e-7, 1);
  model.updateResidual();

  Real pot_energy = model.getEnergy("potential");

  if (std::abs(pot_energy - 7.37343e-06) > 1e-5)
    return EXIT_FAILURE;

  finalize();
  return 0;
}