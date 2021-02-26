#===============================================================================
# @file   model_couplers.cmake
#
# @author Mohit Pundir <mohit.pundir@epfl.ch>
#
# @date creation: Sun Sep 30 2018
# @date last modification: Sun Sep 28 2018
#
# @brief  package description for model couplers
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

package_declare(model_couplers ADVANCED
  DESCRIPTION "Use Model Couplers package of Akantu")

package_declare_sources(model_couplers 
  model/model_couplers/coupler_solid_phasefield.hh
  model/model_couplers/coupler_solid_phasefield.cc
  )

package_declare_documentation_files(model_couplers
  #
  )

package_declare_documentation(model_couplers
  "This package activates the modle couplers within Akantu. "
  "It has no additional dependencies."
  )
