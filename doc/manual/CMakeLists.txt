#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Sep 03 2010
# @date last modification: Tue Jan 19 2016
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

#------------------------------------------------------------------------------#
function(get_doc_label pkg_name label)
  string(REPLACE "_" "-" _pkg_tex_label_dash ${pkg_name})
  string(TOLOWER "${_pkg_tex_label_dash}" _pkg_tex_label)
  set(TOLOWER "${_pkg_tex_label_dash}" _pkg_tex_label)
  set(${label} "pkg:${_pkg_tex_label}" PARENT_SCOPE)
endfunction()

function(get_doc_package_name pkg_name pkg_tex)
  _package_get_option_name(${pkg_name} _opt_name)
  string(REPLACE "_" "\\_" _pkg_tex "${_opt_name}")
  set(${pkg_tex} ${_pkg_tex} PARENT_SCOPE)
endfunction()

#------------------------------------------------------------------------------#
function(generate_package_dependency_tex_doc pkg_name FILENAME LEVEL)
  _package_get_option_name(${pkg_name} _opt_name)
  string(REPLACE "_" "\\_" _pkg_tex "${_opt_name}")
  get_doc_label(${pkg_name} _pkg_tex_label)

  file(APPEND ${FILENAME} "\\AkantuPackageNameWithLabel{${_pkg_tex}}{${_pkg_tex_label}}{${LEVEL}}\\allowbreak\\xspace")

  math(EXPR _sub_level "${LEVEL}+1")

  _package_get_dependencies(${pkg_name} PRIVATE _private_dependencies)
  _package_get_dependencies(${pkg_name} INTERFACE _interface_dependencies)

  foreach(_dep_pkg_name ${_private_dependencies} ${_interface_dependencies})
    generate_package_dependency_tex_doc(${_dep_pkg_name} ${FILENAME} ${_sub_level})
  endforeach()
endfunction()

#------------------------------------------------------------------------------#
function(generate_package_tex_doc pkg_name FILENAME)
  get_doc_package_name(${pkg_name} _pkg_tex)
  get_doc_label(${pkg_name} _pkg_tex_label)

  file(APPEND ${FILENAME} "\n\\begin{AkantuPackage}{${_pkg_tex}}{${_pkg_tex_label}}")

  _package_get_documentation(${pkg_name} _doc)

  if (_doc)
    file(APPEND ${FILENAME} "${_doc}")
  else()
    _package_get_filename(${pkg_name} _file_path)
    _package_get_real_name(${pkg_name} _pkg)
    get_filename_component(_file ${_file_path} NAME)
    string(REPLACE "_" "\\_" _escaped_file "${_file}")
    string(REPLACE "_" "\\_" _escaped_pkg "${_pkg}")

    set(_missing_doc
      "{\\color{red} TODO}: No Documentation in {\\color{blue} \\href{${_file_path}}{${_escaped_file}}}"
      ""
      "looking for the sequence: "
      "\\begin{cmake}"
      "\\package_declare_documentation("
      "  ${_escaped_pkg}"
      "  \"documentation text\""
      "  )"
      "\\end{cmake}")

    set(_missing_doc_str "")
    foreach(_str ${_missing_doc})
      set(_missing_doc_str "${_missing_doc_str}\n${_str}")
    endforeach()

    file(APPEND ${FILENAME} "${_missing_doc_str}")
  endif()

  _package_get_dependencies(${pkg_name} PRIVATE _private_dependencies)
  _package_get_dependencies(${pkg_name} INTERFACE _interface_dependencies)

  if(_private_dependencies OR _interface_dependencies)
    file(APPEND ${FILENAME} "\n\\begin{AkantuPackageDependencies}")
    foreach(_dep_pkg_name ${_private_dependencies} ${_interface_dependencies})
      generate_package_dependency_tex_doc(${_dep_pkg_name} ${FILENAME} 1)
    endforeach()
    file(APPEND ${FILENAME} "\n\\end{AkantuPackageDependencies}")
  endif()

  file(APPEND ${FILENAME} "\n\\end{AkantuPackage}
")

endfunction()

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
set(AKANTU_MANUAL_SOURCE_FOLDER ${CMAKE_CURRENT_BINARY_DIR}/latex)

set(DOC_DEPS_TEX_FILENAME "${AKANTU_MANUAL_SOURCE_FOLDER}/manual-packages-doc.tex")
file(WRITE ${DOC_DEPS_TEX_FILENAME} "")

find_program(RUBBER_EXECUTABLE rubber)
if (NOT RUBBER_EXECUTABLE)
  message(ERROR "Manual cannot be built without rubber latex compiler")
endif()
mark_as_advanced(RUBBER_EXECUTABLE)

package_get_all_documentation_files(_manual_files)

set(AKANTU_MANUAL_FILES_DEPEND)
set(AKANTU_MANUAL_FILES_COPY_COMMAND)

file(MAKE_DIRECTORY ${AKANTU_MANUAL_SOURCE_FOLDER})
foreach(_f  ${_manual_files})
  file(RELATIVE_PATH _rel_f ${CMAKE_CURRENT_SOURCE_DIR} "${_f}")
  list(APPEND AKANTU_MANUAL_FILES_DEPEND ${_f})
  list(APPEND AKANTU_MANUAL_FILES_COPY_COMMAND
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/${_rel_f}" "${_rel_f}")
endforeach()

set(MANUAL_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/akantu_ug-v${AKANTU_VERSION}.pdf)

add_custom_command(
  OUTPUT ${MANUAL_OUTPUT}
  DEPENDS ${AKANTU_MANUAL_FILES_DEPEND} ${DOC_DEPS_TEX_FILENAME}
  ${AKANTU_MANUAL_FILES_COPY_COMMAND}
  COMMAND ${RUBBER_EXECUTABLE} -dfq manual
  COMMAND ${CMAKE_COMMAND} -E copy manual.pdf ${MANUAL_OUTPUT}
  WORKING_DIRECTORY ${AKANTU_MANUAL_SOURCE_FOLDER}
  COMMENT "Compiling the user's manual"
  )

set_directory_properties(latex/$)

add_custom_target(manual ALL DEPENDS ${MANUAL_OUTPUT})

install(FILES ${MANUAL_OUTPUT}
  DESTINATION ${CMAKE_INSTALL_DOCDIR})

package_get_all_activated_packages(_package_list)
foreach (_pkg_name ${_package_list})
  generate_package_tex_doc(${_pkg_name} ${DOC_DEPS_TEX_FILENAME})
endforeach()

configure_file(version-definition.tex.in "${AKANTU_MANUAL_SOURCE_FOLDER}/version-definition.tex" @ONLY)
