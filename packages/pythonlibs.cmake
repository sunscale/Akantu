#===============================================================================
# @file   pythonlibs.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @brief  package description for the python library
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
set(Python_ADDITIONAL_VERSIONS 2.7)

package_declare(PythonLibs EXTERNAL DESCRIPTION "Akantu's python interface"
  DEPENDS numpy
  EXTRA_PACKAGE_OPTIONS PREFIX PYTHON FOUND PYTHONLIBS_FOUND
  )

package_declare_sources(Pythonlibs
  python/python_functor.cc
  python/python_functor.hh
  python/python_functor_inline_impl.cc
  model/boundary_condition_python_functor.hh
  model/boundary_condition_python_functor.cc
  model/solid_mechanics/materials/material_python/material_python.cc
  model/solid_mechanics/materials/material_python/material_python.hh
  )