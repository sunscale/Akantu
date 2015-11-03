#===============================================================================
# @file   80_mpi.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Nov 21 2011
# @date last modification: Sat Jun 14 2014
#
# @brief  package description for mpi
#
# @section LICENSE
#
# Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

package_declare(MPI EXTERNAL
  DESCRIPTION "Add MPI support in akantu"
  EXTRA_PACKAGE_OPTIONS PREFIX MPI_C MPI
  DEPENDS scotch)

package_declare_sources(MPI
  synchronizer/mpi_type_wrapper.hh
  synchronizer/static_communicator_mpi.cc
  synchronizer/static_communicator_mpi_inline_impl.hh
  synchronizer/static_communicator_mpi.hh
  )

get_cmake_property(_all_vars VARIABLES)
foreach(_var ${_all_vars})
  if(_var MATCHES "^MPI_.*")
    mark_as_advanced(${_var})
  endif()
endforeach()

package_declare_documentation(MPI
  "This is a meta package providing access to MPI."
  ""
  "Under Ubuntu (14.04 LTS) the installation can be performed using the commands:"
  "\\begin{command}"
  "  > sudo apt-get install libopenmpi-dev"
  "\\end{command}"
  ""
  "Under Mac OS X the installation requires the following steps:"
  "\\begin{command}"
  "  > sudo port install mpich-devel"
  "\\end{command}"
  )
