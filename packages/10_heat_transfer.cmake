#===============================================================================
# @file   10_heat_transfer.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
#
# @date creation: Mon Nov 21 2011
# @date last modification: Thu Jun 12 2014
#
# @brief  package description for heat transfer
#
# @section LICENSE
#
# Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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

option(AKANTU_HEAT_TRANSFER "Use Heat Transfer package of Akantu" OFF)

set(AKANTU_HEAT_TRANSFER_FILES
  model/heat_transfer/heat_transfer_model.cc
  model/heat_transfer/heat_transfer_model.hh
  model/heat_transfer/heat_transfer_model_inline_impl.cc
  )

set(AKANTU_HEAT_TRANSFER_DOC
  manual/manual-heattransfermodel.tex
  )

set(AKANTU_HEAT_TRANSFER_TESTS
   test_heat_transfer_model_cube3d
   test_heat_transfer_model_cube3d_pbc
   test_heat_transfer_model_square2d_pbc
   test_heat_transfer_model_square2d
   test_heat_transfer_model_cube3d_istropic_conductivity
   )

set(AKANTU_HEAT_TRANSFER_MANUAL_FILES
  manual-heattransfermodel.tex
)

set(AKANTU_HEAT_TRANSFER_DOCUMENTATION 
"
This package activates the heat transfer model within Akantu. 
It has no additional dependencies.
")