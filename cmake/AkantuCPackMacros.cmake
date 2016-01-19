#===============================================================================
# @file   AkantuCPack.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Wed Jan 18 2015
#
# @brief  macros to help for cpack
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

# ==============================================================================
# let start the ugly things to list what should not be in the package
# ==============================================================================
function(generate_cpack_ignore_list ignore_list)
  message(STATUS "Generating CPack ignore list...")

  # get a full list of all the files include in the source folder
  file(GLOB _first_level "${PROJECT_SOURCE_DIR}/*")
  set(_all_files "${PROJECT_SOURCE_DIR}" ${_first_level})
  foreach(_path ${_first_level})
    if(IS_DIRECTORY "${_path}" AND
        NOT _path MATCHES "build.*" AND
        NOT _path MATCHES "\\.git.*")
      file(GLOB_RECURSE _second_level "${_path}/*")
      list(APPEND _all_files ${_second_level})
    endif()
  endforeach()

  set(_dirs)
  foreach(_file ${_all_files})
    get_filename_component(_dir ${_file} DIRECTORY)
    list(APPEND _dirs ${_dir})
  endforeach()
  list(REMOVE_DUPLICATES _dirs)
  list(APPEND _all_files ${_dirs})

  # getting list of all the files that should be in the source package
  package_get_files_for_package(_all_package_files)
  package_get_all_activated_packages(_activated_packages)
  foreach(_pkg ${_activated_packages})
    _package_get_filename(${_pkg} _file_name)

    _package_get_source_files(${_pkg} _srcs _pub_hdrs _priv_hdrs)
    _package_get_variable(TESTS_FILES ${_pkg} _tests_files)
    # adding the source directory
    set(_need_source_folder)
    foreach(_src ${_srcs} ${_pub_hdrs} ${_priv_hdrs} ${_tests_files})
      list(APPEND _need_source_folder ${PROJECT_SOURCE_DIR}/${_src})
    endforeach()

    _package_get_documentation_files(${_pkg} _doc_files)
    # adding the manual directory
    set(_all_docs)
    _package_get_manual_folder(${_pkg} _doc_folder)
    foreach(_file ${_doc_files})
      list(APPEND _all_docs ${_doc_folder}/${_file})
    endforeach()

    _package_get_variable(EXAMPLES_FILES ${_pkg} _examples_files)
    _package_get_variable(EXTRA_FILES ${_pkg} _extra_files)

    # split the set in 2 for debug reasons
    set(_package_files
      ${_file_name} ${_need_source_folder}
      ${_examples_files} ${_extra_files} ${_all_docs}
      )

    list(APPEND _all_package_files ${_package_files})
  endforeach()

  # generate ignore list
  set(_ignore_list)
  set(_keep_dirs)
  set(_ignore_dirs)
  foreach(_file ${_all_files})
    set(_found FALSE)
    foreach(_pkg_file ${_all_package_files})
      if(IS_DIRECTORY "${_file}" AND
          _pkg_file MATCHES "${_file}/.*")
        set(_found TRUE)
        list(APPEND _keep_dirs "${_file}")
        break()
      elseif(_pkg_file MATCHES "${_file}")
        set(_found TRUE)
        break()
      endif()
    endforeach()

    if(NOT _found)
      list(APPEND _ignore_list "${_file}")
      if(IS_DIRECTORY "${_file}")
        list(APPEND _ignore_dirs "${_file}")
      endif()
    endif()
  endforeach()

  # Save CMakeLists.txt in folder that are kept
  foreach(_dir ${_keep_dirs})
    if(EXISTS "${_dir}/CMakeLists.txt")
      list(REMOVE_ITEM _ignore_list "${_dir}/CMakeLists.txt")
    endif()
  endforeach()

  set(_tmp_ignore_list ${_ignore_list})
  set(_ignore_list)
  # cleaning the ignore
  foreach(_file ${_tmp_ignore_list})
    set(_found FALSE)
    foreach(_dir ${_ignore_dirs})
      if(_file MATCHES "${_dir}/.*")
        set(_found TRUE)
        break()
      endif()
    endforeach()

    if(NOT _found)
      if(IS_DIRECTORY ${_file})
        list(APPEND _ignore_list "${_file}/")
      else()
        list(APPEND _ignore_list "${_file}")
      endif()
    endif()
  endforeach()

  list(SORT _ignore_list)
  list(REMOVE_DUPLICATES _ignore_list)

  set(${ignore_list} ${_ignore_list} PARENT_SCOPE)
endfunction()
# ==============================================================================
# Let's hope that after all this the list is complete and not to heavy for cpack
# ==============================================================================
