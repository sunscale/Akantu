/**
 * @file   boundary_functions.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 *
 * @brief  functions for boundaries
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
// akantu
#include "aka_common.hh"
#include "solid_mechanics_model.hh"

namespace akantu {

Real integrateResidual(const std::string & sub_boundary_name,
		       const SolidMechanicsModel & model,
		       UInt dir);

/// this is a fix so that all subboundaries exist on all procs
void boundaryFix(Mesh & mesh,
		 const std::vector<std::string> & sub_boundary_names);

} // namespace akantu
