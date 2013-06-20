#include "boundary_functions.hh"

__BEGIN_SIMTOOLS__

/* -------------------------------------------------------------------------- */
Real integrateResidual(const std::string & sub_boundary_name,
		       const SolidMechanicsModel & model,
		       UInt dir) {
  Real int_res = 0.;
  
  const Mesh & mesh = model.getMesh();
  const Array<Real> & residual = model.getResidual();

  // do not need try catch, as all subboundaries should be everywhere. 
  //  try {
  const SubBoundary & boundary = mesh.getSubBoundary(sub_boundary_name);
  SubBoundary::nodes_const_iterator nit  = boundary.nodes_begin();
  SubBoundary::nodes_const_iterator nend = boundary.nodes_end();
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
  
  const Array<UInt> empty_node_vector;
  Boundary & boundary = mesh.getBoundary();

  std::vector<std::string>::const_iterator it  = sub_boundary_names.begin();
  std::vector<std::string>::const_iterator end = sub_boundary_names.end();

  for (; it != end; ++it) {
    if (boundary.find(*it) == boundary.end()) {
      boundary.createSubBoundaryFromNodeGroup(*it, empty_node_vector);
    }
  }
}


__END_SIMTOOLS__
