/**
 * @file   boundary_functions.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Tue Dec 02 2014
 * @date last modification: Fri Feb 23 2018
 *
 * @brief
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
#include "boundary_functions.hh"
#include "communicator.hh"
#include "element_group.hh"
#include "node_group.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
Real integrateResidual(const std::string & sub_boundary_name,
                       const SolidMechanicsModel & model, UInt dir) {
  Real int_res = 0.;

  const Mesh & mesh = model.getMesh();
  const Array<Real> & residual = model.getInternalForce();

  const ElementGroup & boundary = mesh.getElementGroup(sub_boundary_name);
  for (auto & node : boundary.getNodeGroup().getNodes()) {
    bool is_local_node = mesh.isLocalOrMasterNode(node);
    if (is_local_node) {
      int_res += residual(node, dir);
    }
  }

  mesh.getCommunicator().allReduce(int_res, SynchronizerOperation::_sum);
  return int_res;
}

/* -------------------------------------------------------------------------- */
void boundaryFix(Mesh & mesh,
                 const std::vector<std::string> & sub_boundary_names) {
  for (auto && eg : sub_boundary_names) {
    if (mesh.elementGroupExists(eg)) {
      mesh.createElementGroup(eg, mesh.getSpatialDimension() -
                                      1); // empty element group
    }
  }
}

} // namespace akantu
