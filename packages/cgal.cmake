#===============================================================================
# @file   80_cgal.cmake
#
# @author Lucas Frérot <lucas.frerot@epfl.ch>
#
# @date creation: Thu Feb 19 2015
# @date last modification: Tue Feb 24 2015
#
# @brief  package description for CGAL
#
# @section LICENSE
#
# Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

package_declare(CGAL EXTERNAL
  DESCRIPTION "Add CGAL support in akantu"
  COMPILE_FLAGS "-frounding-math"
  )

package_declare_sources(CGAL
  geometry/mesh_abstract_constructor.hh
  geometry/mesh_abstract_constructor.cc
  geometry/mesh_tree_constructor.hh
  geometry/mesh_tree_constructor_tmpl.hh
  geometry/mesh_tree_container.hh
  geometry/mesh_tree_container.cc
  geometry/tree_type_helper.hh
  )

package_boost_component_needed(system)
## Adding CGAL library
#find_package(CGAL COMPONENTS Core)
#if (NOT CGAL_FOUND)
#message(STATUS "This project requires the CGAL library, and will not be compiled.")
#return()  
#endif()
