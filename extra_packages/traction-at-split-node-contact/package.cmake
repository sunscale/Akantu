#===============================================================================
# @file   package.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Tue Dec 02 2014
# @date last modification: Wed Feb 03 2016
#
# @brief 
#
# @section LICENSE
#
# Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================
package_declare(traction_at_split_node_contact
  DESCRIPTION "The super contact of David"
  DEPENDS iohelper)

package_declare_sources(traction_at_split_node_contact
  common/synchronized_array.cc
  common/parameter_reader.cc
  common/manual_restart.cc
  functions/boundary_functions.cc
  ntn_contact/ntn_base_contact.cc
  ntn_contact/ntn_contact.cc
  ntn_contact/ntrf_contact.cc
  ntn_contact/mIIasym_contact.cc
  ntn_contact/ntn_base_friction.cc
  ntn_contact/friction_regularisations/ntn_fricreg_no_regularisation.cc
  ntn_contact/friction_regularisations/ntn_fricreg_rubin_ampuero.cc
  ntn_contact/friction_regularisations/ntn_fricreg_simplified_prakash_clifton.cc
  ntn_contact/ntn_initiation_function.cc

  tasn_contact.hh

  # headers
  common/synchronized_array.hh
  common/parameter_reader.hh
  common/manual_restart.hh

  functions/boundary_functions.hh
  functions/node_filter.hh

  boundary_conditions/force_based_dirichlet.hh
  boundary_conditions/spring_bc.hh
  boundary_conditions/inclined_flat_dirichlet.hh

  ntn_contact/ntn_base_contact.hh
  ntn_contact/ntn_contact.hh
  ntn_contact/ntrf_contact.hh
  ntn_contact/mIIasym_contact.hh
  ntn_contact/ntn_base_friction.hh

  ntn_contact/friction_regularisations/ntn_fricreg_no_regularisation.hh
  ntn_contact/friction_regularisations/ntn_fricreg_rubin_ampuero.hh
  ntn_contact/friction_regularisations/ntn_fricreg_simplified_prakash_clifton.hh

  ntn_contact/friction_laws/ntn_friclaw_coulomb.hh
  ntn_contact/friction_laws/ntn_friclaw_coulomb_tmpl.hh
  ntn_contact/friction_laws/ntn_friclaw_linear_slip_weakening.hh
  ntn_contact/friction_laws/ntn_friclaw_linear_slip_weakening_tmpl.hh
  ntn_contact/friction_laws/ntn_friclaw_linear_slip_weakening_no_healing.hh
  ntn_contact/friction_laws/ntn_friclaw_linear_slip_weakening_no_healing_tmpl.hh
  ntn_contact/friction_laws/ntn_friclaw_linear_cohesive.hh
  ntn_contact/friction_laws/ntn_friclaw_linear_cohesive_tmpl.hh

  ntn_contact/ntn_friction.hh
  ntn_contact/ntn_friction_tmpl.hh
  ntn_contact/ntrf_friction.hh
  ntn_contact/ntrf_friction_tmpl.hh

  ntn_contact/ntn_initiation_function.hh

  # inlines
  common/synchronized_array_inline_impl.hh
  ntn_contact/ntn_base_contact_inline_impl.hh
  )
