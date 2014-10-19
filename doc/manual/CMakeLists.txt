#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Jun 09 2014
# @date last modification: Thu Jul 03 2014
#
# @brief  Build the documentation
#
# @section LICENSE
#
# Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#------------------------------------------------------------------------------#
function(get_doc_label PKG LABEL)
  string(REPLACE "_" "-" _pkg_tex_label_dash ${PKG})
  string(TOLOWER "${_pkg_tex_label_dash}" _pkg_tex_label)
  set(TOLOWER "${_pkg_tex_label_dash}" _pkg_tex_label)
  set(${LABEL} "pkg:${_pkg_tex_label}" PARENT_SCOPE)
endfunction()

#------------------------------------------------------------------------------#
function(generate_package_dependency_tex_doc PKG FILENAME LEVEL)
  string(REPLACE "_" "\\_" _pkg_tex "${PKG}")
  get_doc_label(${PKG} _pkg_tex_label)

  file(APPEND ${FILENAME} "\\AkantuPackageNameWithLabel{${_pkg_tex}}{${_pkg_tex_label}}{${LEVEL}}")

  math(EXPR _sub_level "${LEVEL}+1")
  foreach(_pkg ${${PKG}_dependencies})
    generate_package_dependency_tex_doc(${_pkg} ${FILENAME} ${_sub_level})
  endforeach()
endfunction()

#------------------------------------------------------------------------------#
function(generate_package_tex_doc PKG FILENAME)
  string(REPLACE "_" "\\_" _pkg_tex "${PKG}")
  get_doc_label(${PKG} _pkg_tex_label)

  file(APPEND ${FILENAME} "
\\begin{AkantuPackage}{${_pkg_tex}}{${_pkg_tex_label}}")

  if (${PKG}_DOCUMENTATION)
    file(APPEND ${FILENAME} "${${PKG}_DOCUMENTATION}")
  else()
    string(TOLOWER "${PKG}" _l_pkg)
    string(REPLACE "akantu_" "" _l_pkg_short ${_l_pkg})
    string(REPLACE "_" "\\_" _l_pkg_tex ${_l_pkg_short})

    file(GLOB _file_path "${PROJECT_SOURCE_DIR}/packages/*_${_l_pkg_short}.cmake")
    file(APPEND ${FILENAME} "{\\color{red} TODO}: No Documentation in {\\color{blue} \\href{${_file_path}}{${_l_pkg_tex}.cmake}}\\\\")
    file(APPEND ${FILENAME} "looking for the sequence: \\\\set(${_pkg_tex}\\_DOCUMENTATION \" documentation text \" )\\\\" )
  endif()

  set(_n_dependencies 0)
  if (${PKG}_dependencies)
    list(LENGTH ${PKG}_dependencies _n_dependencies)
  endif()

  if(_n_dependencies GREATER 0)
    file(APPEND ${FILENAME} "
\\begin{AkantuPackageDependencies}")
    foreach(_pkg ${${PKG}_dependencies})
      generate_package_dependency_tex_doc(${_pkg} ${FILENAME} 1)
    endforeach()

    file(APPEND ${FILENAME} "
\\end{AkantuPackageDependencies}")
  endif()

  file(APPEND ${FILENAME} "
\\end{AkantuPackage}
")

endfunction()

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
set(DOC_DEPS_TEX_FILENAME "${CMAKE_CURRENT_BINARY_DIR}/manual-packages-doc.tex")
file(WRITE ${DOC_DEPS_TEX_FILENAME} "")

find_program(RUBBER_EXECUTABLE rubber)
if (NOT RUBBER_EXECUTABLE)
  message(ERROR "Manual cannot be built without rubber latex compiler")
endif()
mark_as_advanced(RUBBER_EXECUTABLE)


set(AKANTU_MANUAL_FILES_DEPEND)
set(AKANTU_MANUAL_FILES_COPY_COMMAND)

set(_relative_manual_files)
foreach(_f  ${AKANTU_MANUAL_FILES})
  file(RELATIVE_PATH _rel_f ${CMAKE_CURRENT_SOURCE_DIR} "${PROJECT_SOURCE_DIR}/${_f}")
  list(APPEND _relative_manual_files ${_rel_f})
  list(APPEND AKANTU_MANUAL_FILES_DEPEND ${_rel_f})
  list(APPEND AKANTU_MANUAL_FILES_COPY_COMMAND
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/${_rel_f} ${_rel_f})
endforeach()

set(MANUAL_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/manual.pdf)

add_custom_command(
  OUTPUT ${MANUAL_OUTPUT}
  DEPENDS ${AKANTU_MANUAL_FILES_DEPEND} ${DOC_DEPS_TEX_FILENAME}
  ${AKANTU_MANUAL_FILES_COPY_COMMAND}
  COMMAND ${RUBBER_EXECUTABLE} -dfq -Wall manual
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Compiling the user's manual"
  )

add_custom_target(manual ALL DEPENDS ${MANUAL_OUTPUT})

install(FILES ${MANUAL_OUTPUT} DESTINATION share/akantu/doc)

foreach (_pkg ${AKANTU_PACKAGE_LIST})
  generate_package_tex_doc(AKANTU_${_pkg} ${DOC_DEPS_TEX_FILENAME})
endforeach()
