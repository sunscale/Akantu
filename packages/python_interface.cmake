#===============================================================================
# @file   python_interface.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Tue Nov 29 2011
# @date last modification: Fri Jan 22 2016
#
# @brief  package description for the python interface
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

package_declare(python_interface
  DESCRIPTION "Akantu's python interface"
  DEPENDS PythonLibs)

package_declare_sources(python_interface
  python/python_functor.cc
  python/python_functor.hh
  python/python_functor_inline_impl.cc
  model/boundary_condition_python_functor.hh
  model/boundary_condition_python_functor.cc
  model/solid_mechanics/materials/material_python/material_python.cc
  model/solid_mechanics/materials/material_python/material_python.hh
  )


package_set_package_system_dependency(python_interface deb-src swig3.0)

package_declare_documentation(python_interface
  "This package enables the python interface of Akantu. It relies on swig3.0 to generate the code"
  ""
  "Under Ubuntu (14.04 LTS) the installation can be performed using the commands:"
  "\\begin{command}"
  "  > sudo apt-get install swig3.0"
  "\\end{command}"
  ""
)

package_declare_documentation_files(python_interface
  manual-python.tex
  )
