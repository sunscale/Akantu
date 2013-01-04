find_library(PARADIS_LIBRARIES paradis
  PATHS ${PARADIS_DIR}
  PATH_SUFFIXES bin)

find_path(PARADIS_INCLUDE_PATH ParadisGen.h
  PATHS ${PARADIS_DIR}
  PATH_SUFFIXES include
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PARADIS DEFAULT_MSG
  PARADIS_LIBRARIES PARADIS_INCLUDE_PATH)

if (NOT PARADIS_FOUND)
  set(PARADIS_DIR "" CACHE PATH "Location of PARADIS library")
endif(NOT PARADIS_FOUND)