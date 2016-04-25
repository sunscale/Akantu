find_path(METIS_INCLUDE_DIR metis.h
  PATHS "${METIS_DIR}"
  ENV METIS_DIR
  PATH_SUFFIXES include
  )

find_library(METIS_LIBRARY NAMES metis
  PATHS "${METIS_DIR}"
  ENV METIS_DIR
  PATH_SUFFIXES lib
  )

mark_as_advanced(METIS_LIBRARY METIS_INCLUDE_DIR)

#===============================================================================
include(FindPackageHandleStandardArgs)
if(CMAKE_VERSION VERSION_GREATER 2.8.12)
  if(METIS_INCLUDE_DIR)
    file(STRINGS ${METIS_INCLUDE_DIR}/metis.h _versions
      REGEX "^#define\ +METIS_VER_(MAJOR|MINOR|SUBMINOR) .*")
    foreach(_ver ${_versions})
      string(REGEX MATCH "METIS_VER_(MAJOR|MINOR|SUBMINOR) *([0-9.]+)" _tmp "${_ver}")
      set(_metis_${CMAKE_MATCH_1} ${CMAKE_MATCH_2})
    endforeach()
    set(METIS_VERSION "${_metis_MAJOR}.${_metis_MINOR}" CACHE INTERNAL "")
  endif()

  find_package_handle_standard_args(METIS
    REQUIRED_VARS
      METIS_LIBRARY
      METIS_INCLUDE_DIR
    VERSION_VAR
      METIS_VERSION)
else()
  find_package_handle_standard_args(METIS DEFAULT_MSG
    METIS_LIBRARY METIS_INCLUDE_DIR)
endif()
