/**
 * @file   boundary_functions.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Tue Dec 02 2014
 * @date last modification: Fri Jan 22 2016
 *
 * @brief  
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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

#include "boundary_functions.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
Real integrateResidual(const std::string & sub_boundary_name,
		       const SolidMechanicsModel & model,
		       UInt dir) {
  Real int_res = 0.;
  
  const Mesh & mesh = model.getMesh();
  const Array<Real> & residual = model.getResidual();

  // do not need try catch, as all subboundaries should be everywhere. 
  //  try {
  const ElementGroup & boundary = mesh.getElementGroup(sub_boundary_name);
  ElementGroup::const_node_iterator nit  = boundary.node_begin();
  ElementGroup::const_node_iterator nend = boundary.node_end();
  for (; nit != nend; ++nit) {
    bool is_local_node = mesh.isLocalOrMasterNode(*nit);
    if (is_local_node) {
      int_res += residual(*nit, dir);
    }
  }
  // } catch(debug::Exception e) {
  //   // AKANTU_DEBUG_ERROR("Error computing integrateResidual. Cannot get SubBoundary: " 
  //   // 		       << sub_boundary_name << " [" << e.what() << "]");
  // }

  StaticCommunicator::getStaticCommunicator().allReduce(&int_res, 1, _so_sum);
  return int_res;
}

/* -------------------------------------------------------------------------- */
void boundaryFix(Mesh & mesh, 
		 const std::vector<std::string> & sub_boundary_names) {
  
   std::vector<std::string>::const_iterator it  = sub_boundary_names.begin();
  std::vector<std::string>::const_iterator end = sub_boundary_names.end();

  for (; it != end; ++it) {
    if (mesh.element_group_find(*it) == mesh.element_group_end()) {
      mesh.createElementGroup(*it,mesh.getSpatialDimension()-1); // empty element group
    }
  }
}


__END_AKANTU__
