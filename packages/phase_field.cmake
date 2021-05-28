#===============================================================================
# @file   phase_field.cmake
#
# @author Mohit Pundir <mohit.pundir@epfl.ch>
#
# @date creation: Sun Sep 30 2018
# @date last modification: Sun Sep 30 2018
#
# @brief  package description for phase field model
#
# @section LICENSE
#
# Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
# Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
# Solides)
#
# Akantu is free  software: you can redistribute it and/or  modify it under the
# terms  of the  GNU Lesser  General Public  License as  published by  the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
# details.
#
# You should  have received  a copy  of the GNU  Lesser General  Public License
# along with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================

package_declare(phase_field
  DEPENDS model_couplers
  DESCRIPTION "Use Phase Field package of Akantu")

package_declare_sources(phase_field
  model/phase_field/phasefield.cc
  model/phase_field/phasefield.hh
  model/phase_field/phasefield_inline_impl.cc
  model/phase_field/phasefield_selector.hh
  model/phase_field/phasefield_selector_tmpl.hh

  model/phase_field/phasefields/phasefield_exponential.hh
  model/phase_field/phasefields/phasefield_exponential.cc
  
  model/phase_field/phase_field_model.cc
  model/phase_field/phase_field_model.hh
  model/phase_field/phase_field_model_inline_impl.cc
  )

package_declare_documentation_files(phase_field
  #
  )

package_declare_documentation(phase_field
  "This package activates the phase field model within Akantu. "
  )
