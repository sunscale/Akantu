#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Fri Jan 30 2015
#
# @brief  Build the documentation
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

if(AKANTU_DOCUMENTATION_DEVELOPER_MANUAL)
  add_subdirectory(dev-doc)
endif()


file(GLOB_RECURSE __all_files "*.*")
set(_all_files)
foreach(_file ${__all_files})
  if("${_file}" MATCHES  "${PROJECT_SOURCE_DIR}/doc/manual")
    file(RELATIVE_PATH __file ${PROJECT_SOURCE_DIR} ${_file})
    list(APPEND _all_files ${__file})
  endif()
endforeach()

set(AKANTU_MANUAL_FILES)
set(_akantu_manual_not_found_files)
foreach(_pkg ${${_project}_PACKAGE_SYSTEM_PACKAGES_ON})
  string(TOUPPER "${_pkg}" _pkg)
  foreach(_f ${AKANTU_${_pkg}_MANUAL_FILES})
    #check if some file are registered but not present
    list(FIND _all_files doc/manual/${_f} _ret)
    if(_ret EQUAL -1)
      list(APPEND _akantu_manual_not_found_files "${_pkg}: doc/manual/${_f} ")
    else()
      list(APPEND AKANTU_MANUAL_FILES doc/manual/${_f})
    endif()
  endforeach()
endforeach()

if(_akantu_manual_not_found_files)
  message("")
  message("******************************************************")
  message("There are files registered in packages but not present")
  message("******************************************************")
  foreach(_file ${_akantu_manual_not_found_files})
    message(${_file})
  endforeach()
  message("******************************************************")
  message(FATAL_ERROR "abort")
endif()

################################################################
#construct the list of files to exclude because not registered
################################################################

set(_akantu_doc_exclude_files)
foreach(_file ${_all_files})
  list(FIND AKANTU_MANUAL_FILES ${_file} _ret)
  if(_ret EQUAL -1)
    list(APPEND _akantu_doc_exclude_files /${_file})
  endif()
endforeach()

list(REMOVE_ITEM _akantu_doc_exclude_files /doc/manual/CMakeLists.txt)

set(AKANTU_DOC_EXCLUDE_FILES ${_akantu_doc_exclude_files}
  CACHE INTERNAL "Documentation files to excluse from Akantu Package" FORCE)

if (AKANTU_DOCUMENTATION_MANUAL)
  add_subdirectory(manual)
endif()
