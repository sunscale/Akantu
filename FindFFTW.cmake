
find_library(FFTW_LIBRARIES fftw
   PATHS ${FFTW_DIR}
   PATH_SUFFIXES fftw/.libs/ lib
   )

find_path(FFTW_INCLUDE_PATH fftw.h
  PATHS ${FFTW_DIR}
  PATH_SUFFIXES include fftw
  )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW DEFAULT_MSG
  FFTW_LIBRARIES FFTW_INCLUDE_PATH)


if(NOT FFTW_FOUND)
  set(FFTW_DIR "" CACHE PATH "Location of FFTW library.")
endif(NOT FFTW_FOUND)


