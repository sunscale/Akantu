find_library(SCOTCH_LIBRARY scotch
  PATHS ${SCOTCH_DIR}
  PATH_SUFFIXES src/libscotch lib
  )

find_library(SCOTCH_LIBRARY_ERR scotcherr
  PATHS ${SCOTCH_DIR}
  PATH_SUFFIXES src/libscotch lib
  )

find_path(SCOTCH_INCLUDE_PATH scotch.h
  PATHS ${SCOTCH_DIR}
  PATH_SUFFIXES include scotch src/libscotch include/scotch
  )

set(SCOTCH_LIBRARIES ${SCOTCH_LIBRARY} ${SCOTCH_LIBRARY_ERR} CACHE INTERNAL "Libraries for scotch" FORCE)

if(NOT SCOTCH_LIBRARY)
  set(SCOTCH_DIR "" CACHE PATH "Location of Scotch library.")
endif(NOT SCOTCH_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SCOTCH DEFAULT_MSG
  SCOTCH_LIBRARY SCOTCH_LIBRARY_ERR SCOTCH_INCLUDE_PATH)


