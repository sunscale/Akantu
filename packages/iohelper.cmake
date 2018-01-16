#===============================================================================
# @file   iohelper.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Tue Nov 29 2011
# @date last modification: Mon Jan 18 2016
#
# @brief  package description for iohelper
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
package_declare(IOHelper EXTERNAL
  DESCRIPTION "Add IOHelper support in akantu"
  SYSTEM OFF third-party/cmake/iohelper.cmake
  DEFAULT ON)

set(_version "1.1.1")
package_add_third_party_script_variable(IOHelper
  IOHELPER_VERSION ${_version})
package_add_third_party_script_variable(IOHelper
  IOHELPER_GIT "https://c4science.ch/source/iohelper.git")
package_add_third_party_script_variable(IOHelper
  IOHELPER_ARCHIVE "iohelper_${_version}.tar.gz")
