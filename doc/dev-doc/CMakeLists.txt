set(DOXYGEN_INPUT_DOX ${CMAKE_CURRENT_BINARY_DIR}/akantu.dox)
set(DOXYGEN_XML_DIR ${CMAKE_CURRENT_BINARY_DIR}/xml)
set(DOXYGEN_OUTPUT ${DOXYGEN_XML_DIR}/index.xml)

# configured documentation tools and intermediate build results
set(BINARY_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/_build")
# Sphinx cache with pickled ReST documents
set(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_doctrees")
# HTML output directory
set(SPHINX_HTML_DIR "${CMAKE_CURRENT_BINARY_DIR}/html")
set(SPHINX_OUTPUT "${SPHINX_HTML_DIR}/index.html")

# ---------------------------------------------------------------------------- #
# Doxygen                                                                      #
# ---------------------------------------------------------------------------- #
find_package(Doxygen REQUIRED)

set(DOXYGEN_WARNINGS NO)
set(DOXYGEN_QUIET YES)
if(CMAKE_VERBOSE_MAKEFILE)
  set(DOXYGEN_WARNINGS YES)
  set(DOXYGEN_QUIET NO)
endif(CMAKE_VERBOSE_MAKEFILE)

  package_get_all_source_files(
    AKANTU_LIBRARY_SRCS
    AKANTU_LIBRARY_PUBLIC_HDRS
    AKANTU_LIBRARY_PRIVATE_HDRS
    )

package_get_all_include_directories(
  _akantu_include_dirs
  )

package_get_all_external_informations(
  PRIVATE_INCLUDE AKANTU_PRIVATE_EXTERNAL_INCLUDE_DIR
  INTERFACE_INCLUDE AKANTU_INTERFACE_EXTERNAL_INCLUDE_DIR
  LIBRARIES AKANTU_EXTERNAL_LIBRARIES
  )

list(APPEND _akantu_include_dirs
  ${AKANTU_PRIVATE_EXTERNAL_INCLUDE_DIR}
  ${AKANTU_INTERFACE_EXTERNAL_INCLUDE_DIR}
  ${PROJECT_BINARY_DIR}/src)


file(STRINGS ${PROJECT_SOURCE_DIR}/.clang-format AKANTU_TAB_SIZE
  REGEX "^TabWidth: *([0-9]*)"
  )
string(REGEX REPLACE ".*([0-9]+)" "\\1" AKANTU_TAB_SIZE "${AKANTU_TAB_SIZE}")

if (CMAKE_VERSION VERSION_GREATER 3.9.5)
  #set(DOXYGEN_WARNINGS YES)
  #set(DOXYGEN_QUIET NO)
  set(DOXYGEN_STRIP_FROM_PATH ${PROJECT_SOURCE_DIR})
  set(DOXYGEN_STRIP_FROM_INC_PATH ${PROJECT_SOURCE_DIR})
  set(DOXYGEN_TAB_SIZE ${AKANTU_TAB_SIZE})
  set(DOXYGEN_ALIASES
    "rst=\\verbatim embed:rst"
    "endrst=\\endverbatim"
    )
  set(DOXYGEN_WARN_IF_UNDOCUMENTED NO)
  set(DOXYGEN_WARN_IF_DOC_ERROR NO)
  set(DOXYGEN_WARN_AS_ERROR NO)
  set(DOXYGEN_EXCLUDE  "${PROJECT_SOURCE_DIR}/src/common/aka_fwd.hh")
  set(DOXYGEN_RECURSIVE YES)
  set(DOXYGEN_EXCLUDE
    "aka_named_argument.hh"
    )
  set(DOXYGEN_EXAMPLE_PATH "${PROJECT_SOURCE_DIR}/examples")
  set(DOXYGEN_EXAMPLE_RECURSIVE YES)
  set(DOXYGEN_SOURCE_BROWSER NO)
  set(DOXYGEN_CLANG_ASSISTED_PARSING NO)
  #set(DOXYGEN_CLANG_OPTIONS          )
  set(DOXYGEN_CLANG_DATABASE_PATH ${CMAKE_BINARY_DIR})
  set(DOXYGEN_USE_MATHJAX YES)
  set(DOXYGEN_GENERATE_HTML NO)
  set(DOXYGEN_GENERATE_HTMLHELP NO)
  set(DOXYGEN_GENERATE_LATEX NO)
  set(DOXYGEN_GENERATE_XML YES)
  set(DOXYGEN_XML_OUTPUT xml)
  set(DOXYGEN_ENABLE_PREPROCESSING YES)
  set(DOXYGEN_MACRO_EXPANSION YES)
  set(DOXYGEN_INCLUDE_PATH ${_akantu_include_dirs})
#  set(DOXYGEN_PREDEFINED
#    ${AKANTU_DEFINITIONS}
#    "DOXYGEN"
#    "AKANTU_TO_IMPLEMENT()="
#    "DECLARE_NAMED_ARGUMENT()="
#    "OPTIONAL_NAMED_ARGUMENT(n, v)=v"
#    "REQUIRED_NAMED_ARGUMENT(n)="
#    )
  set(DOXYGEN_COLLABORATION_GRAPH NO)
  set(DOXYGEN_UML_LOOK YES)
  set(DOXYGEN_TEMPLATE_RELATIONS YES)
  set(DOXYGEN_CALL_GRAPH YES)
  set(DOXYGEN_CALLER_GRAPH YES)
  set(DOXYGEN_DOT_GRAPH_MAX_NODES 500)
  set(DOXYGEN_SHOW_FILES NO)
  set(DOXYGEN_LOOKUP_CACHE_SIZE 9)

  set(_SRCS
    ${AKANTU_LIBRARY_SRCS}
    ${AKANTU_LIBRARY_PUBLIC_HDRS}
    ${AKANTU_LIBRARY_PRIVATE_HDRS}
    )

  list(REMOVE_ITEM _SRCS
    "${PROJECT_SOURCE_DIR}/src/common/aka_named_argument.hh")

  doxygen_add_docs(doxygen-doc
    ${_SRCS}
    USE_STAMP_FILE
    COMMENT "Building XML documentation with Doxygen in ${DOXYGEN_XML_DIR}"
    )

else()
  string(REGEX REPLACE ";" " " AKANTU_DOXYGEN_DEFINTIONS "${AKANTU_DEFINITIONS};DOXYGEN")
  string(REGEX REPLACE ";" " " AKANTU_DOXYGEN_INCLUDE_DIRS "${_akantu_include_dirs}")

  make_directory(${DOXYGEN_XML_DIR})

  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/akantu.dox.in
    ${DOXYGEN_INPUT_DOX}
    )

  add_custom_command(
    OUTPUT ${DOXYGEN_OUTPUT}
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_INPUT_DOX}
    DEPENDS ${DOXYGEN_INPUT_DOX}
    COMMENT "Building XML documentation with Doxygen in ${DOXYGEN_XML_DIR}"
    )

  add_custom_target(doxygen-doc ALL
    DEPENDS ${DOXYGEN_OUTPUT}
    )

  add_custom_target(doxygen-doc-forced
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_INPUT}
    COMMENT "Building XML documentation with Doxygen (forced) in ${DOXYGEN_XML_DIR}"
    )
endif()
# ---------------------------------------------------------------------------- #
# Sphinx                                                                       #
# ---------------------------------------------------------------------------- #
find_package(Sphinx REQUIRED)

set(SPHINX_VERBOSE_FLAG "-q")
if(CMAKE_VERBOSE_MAKEFILE)
  set(SPHINX_VERBOSE_FLAG)
endif(CMAKE_VERBOSE_MAKEFILE)

# set(AKANTU_IN_READTHEDOC TRUE)
# configure_file(
#   "${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in"
#   "${CMAKE_CURRENT_SOURCE_DIR}/conf_rtd.py"
#   @ONLY)

set(AKANTU_IN_READTHEDOC FALSE)
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in"
  "${CMAKE_CURRENT_BINARY_DIR}/conf.py"
  @ONLY)

set(SPHINX_PARALLEL_FLAG)
if (SPHINX_VERSION VERSION_GREATER_EQUAL 1.7.0)
  set(SPHINX_PARALLEL_FLAG -j auto)
endif()

set(_sphinx_command ${SPHINX_BUILD_EXECUTABLE}
  ${SPHINX_PARALLEL_FLAG}
  ${SPHINX_VERBOSE_FLAG} -b html
  -c "${CMAKE_CURRENT_BINARY_DIR}"
  -d "${SPHINX_CACHE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "${SPHINX_HTML_DIR}"
  )

file(GLOB_RECURSE _SPHINX_SRCS
  "*.rst")

add_custom_command(
  OUTPUT ${SPHINX_OUTPUT}
  COMMAND ${_sphinx_command}
  DEPENDS doxygen-doc ${SPHINX_INPUT} ${_SPHINX_SRCS}
  COMMENT "Building HTML documentation with Sphinx in ${SPHINX_HTML_DIR}"
  )

add_custom_target(sphinx-doc ALL
  DEPENDS ${SPHINX_OUTPUT})

add_custom_target(sphinx-doc-forced
  COMMAND ${_sphinx_command}
  DEPENDS doxygen-doc
  COMMENT "Building HTML documentation with Sphinx (forced) in ${SPHINX_HTML_DIR}"
  )
