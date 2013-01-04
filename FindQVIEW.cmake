find_library(QVIEW_LIBRARIES NAME qview
  PATHS ${QVIEW_DIR} 
  PATH_SUFFIXES lib
  )

#string(REGEX REPLACE ":" ";" DEFAULT_INCLUDE_PATH $ENV{C_INCLUDE_PATH})
#MESSAGE(${DEFAULT_INCLUDE_PATH})
find_path(QVIEW_INCLUDE_PATH libqview.h
  PATHS ${QVIEW_DIR} ENV C_INCLUDE_PATH
  PATH_SUFFIXES include src
  )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(QVIEW DEFAULT_MSG
  QVIEW_LIBRARIES QVIEW_INCLUDE_PATH)


if(NOT QVIEW_FOUND)
  set(QVIEW_DIR "" CACHE PATH "Location of QVIEW library.")
endif(NOT QVIEW_FOUND)

