#===============================================================================
# @file   CTestConfig.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Thu Jan 06 2011
# @date last modification: Thu Aug 21 2014
#
# @brief  configuration for ctest
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

## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
## # The following are required to uses Dart and the Cdash dashboard
##   ENABLE_TESTING()
##   INCLUDE(CTest)
set(CTEST_PROJECT_NAME "Akantu")
set(CTEST_NIGHTLY_START_TIME "06:10:00 EST")
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")

set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "lsmssrv1.epfl.ch")
set(CTEST_DROP_LOCATION "/CDash/submit.php?project=Akantu")
set(CTEST_DROP_SITE_CDASH TRUE)
