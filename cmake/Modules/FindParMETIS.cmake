find_path(PARMETIS_INCLUDE_DIR parmetis.h
  PATHS "${PARMETIS_DIR}"
  ENV PARMETIS_DIR
  PATH_SUFFIXES include
  )

find_library(PARMETIS_LIBRARY NAMES parmetis
  PATHS "${PARMETIS_DIR}"
  ENV PARMETIS_DIR
  PATH_SUFFIXES lib
  )

mark_as_advanced(PARMETIS_LIBRARY PARMETIS_INCLUDE_DIR)

#===============================================================================
include(FindPackageHandleStandardArgs)
if(CMAKE_VERSION VERSION_GREATER 2.8.12)
  if(PARMETIS_INCLUDE_DIR)
    file(STRINGS ${PARMETIS_INCLUDE_DIR}/parmetis.h _versions
      REGEX "^#define\ +PARMETIS_(MAJOR|MINOR|SUBMINOR)_VERSION .*")
    foreach(_ver ${_versions})
      string(REGEX MATCH "PARMETIS_(MAJOR|MINOR|SUBMINOR)_VERSION *([0-9.]+)" _tmp "${_ver}")
      set(_parmetis_${CMAKE_MATCH_1} ${CMAKE_MATCH_2})
    endforeach()
    set(PARMETIS_VERSION "${_parmetis_MAJOR}.${_parmetis_MINOR}" CACHE INTERNAL "")
  endif()

  find_package_handle_standard_args(ParMETIS
    REQUIRED_VARS
      PARMETIS_LIBRARY
      PARMETIS_INCLUDE_DIR
    VERSION_VAR
      PARMETIS_VERSION)
else()
  find_package_handle_standard_args(ParMETIS DEFAULT_MSG
    PARMETIS_LIBRARY PARMETIS_INCLUDE_DIR)
endif()
