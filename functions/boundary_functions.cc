#include "boundary_functions.hh"

__BEGIN_SIMTOOLS__

/* -------------------------------------------------------------------------- */
Real integrateResidual(const std::string & sub_boundary_name,
		       const SolidMechanicsModel & model,
		       UInt dir) {
  Real int_res = 0.;
  
  const Mesh & mesh = model.getMesh();
  const Array<Real> & residual = model.getResidual();

  try {
    const SubBoundary & boundary = mesh.getSubBoundary(sub_boundary_name);
    SubBoundary::nodes_const_iterator nit  = boundary.nodes_begin();
    SubBoundary::nodes_const_iterator nend = boundary.nodes_end();
    for (; nit != nend; ++nit) {
      bool is_local_node = mesh.isLocalOrMasterNode(*nit);
      if (is_local_node) {
	int_res += residual(*nit, dir);
      }
    }
  } catch(debug::Exception e) {
    // AKANTU_DEBUG_ERROR("Error computing integrateResidual. Cannot get SubBoundary: " 
    // 		       << sub_boundary_name << " [" << e.what() << "]");
  }

  StaticCommunicator::getStaticCommunicator().allReduce(&int_res, 1, _so_sum);
  return int_res;
}

__END_SIMTOOLS__
