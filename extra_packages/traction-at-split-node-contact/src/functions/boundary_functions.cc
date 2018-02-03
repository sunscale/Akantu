/**
 * @file   boundary_functions.cc
 *
 *
 *
 * @brief  
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

#include "boundary_functions.hh"

namespace akantu {

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


} // namespace akantu
