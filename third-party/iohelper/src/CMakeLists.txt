#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
#
# @date creation: Thu Nov 24 2011
# @date last modification: Wed Nov 13 2013
#
# @brief  main iohelper configuration
#
# @section LICENSE
#
# Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# IOHelper is free  software: you can redistribute it and/or  modify it under the
# terms  of the  GNU Lesser  General Public  License as  published by  the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# IOHelper is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
# details.
#
# You should  have received  a copy  of the GNU  Lesser General  Public License
# along with IOHelper. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================

#===============================================================================
# List of source files
#===============================================================================
set(IOHELPER_COMMON_SRC
  dumper_lammps.cc
  dumper.cc
  dumper_paraview.cc
  dumper_text.cc
  paraview_helper.cc
#  src/reader_restart.cpp
)

set(IOHELPER_COMMON_HEADERS
  field_inline_impl.hh
  dumper_restart.hh
  field_interface.hh
  visitor.hh
  field.hh
  variable_inline_impl.hh
  variable_interface.hh
  variable.hh
  container_array.hh
  dumper_paraview.hh
  dumper_text.hh
  paraview_helper.tcc
  iohelper_common.hh
  io_helper.hh
  dumper.hh
  file_manager.hh
  paraview_helper.hh
  dumper_C_wrapper.h
  base64.hh
  dumper_lammps.hh
  base64_reader.hh
)

#===============================================================================
# Library creation rule
#===============================================================================
add_library(iohelper ${IOHELPER_COMMON_SRC})

# link library with other libraries
target_link_libraries(iohelper ${IOHELPER_EXTERNAL_LIBS})
target_include_directories(iohelper
  PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
            $<INSTALL_INTERFACE:include/iohelper>
  PUBLIC ${ZLIB_INCLUDE_DIR} ${ZLIB_INCLUDE_DIRS}
  )


set_target_properties(iohelper PROPERTIES PUBLIC_HEADER "${IOHELPER_COMMON_HEADERS}")
set_property(TARGET iohelper PROPERTY CXX_STANDARD 11)

export(TARGETS iohelper
  FILE "${CMAKE_BINARY_DIR}/IOHelperLibraryDepends.cmake")
export(PACKAGE IOHelper)

#===============================================================================
# Install rules
#===============================================================================
# Tweak for when IOHelper is a subproject
if(NOT IOHELPER_TARGETS_EXPORT)
  set(IOHELPER_TARGETS_EXPORT IOHelperLibraryDepends)
endif()

install(TARGETS iohelper
  EXPORT ${IOHELPER_TARGETS_EXPORT}
  ARCHIVE DESTINATION lib COMPONENT lib
  LIBRARY DESTINATION lib COMPONENT lib
  RUNTIME DESTINATION lib COMPONENT lib
  PUBLIC_HEADER DESTINATION include/iohelper COMPONENT dev
  )

# Install the export set for use with the install-tree
if("${IOHELPER_TARGETS_EXPORT}" STREQUAL "IOHelperLibraryDepends")
  install(EXPORT IOHelperLibraryDepends DESTINATION  lib/iohelper
    COMPONENT dev)
endif()


set(IOHELPER_INCLUDE_DIRS ${IOHelper_SOURCE_DIR}/src CACHE INTERNAL "Internal include directorie to link with IOHelper as a subproject")
