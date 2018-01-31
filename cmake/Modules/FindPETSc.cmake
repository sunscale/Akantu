# - Try to find PETSc
#  PETSC_FOUND         - system has PETSc
#  PETSC_INCLUDE_DIRS  - the PETSc include directories
#  PETSC_LIBRARIES     - Link these to use PETSc
#  PETSC_VERSION       - Version string (MAJOR.MINOR.SUBMINOR)

if(PETSc_FIND_REQUIRED)
  find_package(PkgConfig REQUIRED)
else()
  find_package(PkgConfig QUIET)
  if(NOT PKG_CONFIG_FOUND)
    return()
  endif()
endif()

pkg_search_module(_petsc PETSc)

if(_petsc_FOUND AND _petsc_VERSION)
  set(PETSC_VERSION ${_petsc_VERSION})
endif()

if(_petsc_FOUND AND (NOT PETSC_LIBRARIES))
  set(_petsc_libs)
  foreach(_lib ${_petsc_LIBRARIES})
    string(TOUPPER "${_lib}" _u_lib)
    find_library(PETSC_LIBRARY_${_u_lib} ${_lib} PATHS ${_petsc_LIBRARY_DIRS})
    list(APPEND _petsc_libs ${PETSC_LIBRARY_${_u_lib}})
    mark_as_advanced(PETSC_LIBRARY_${_u_lib})
  endforeach()

  set(PETSC_LIBRARIES ${_petsc_libs} CACHE INTERNAL "")
  set(PETSC_INCLUDE_DIRS ${_petsc_INCLUDE_DIRS} CACHE INTERNAL "")
  if(NOT TARGET petsc::petsc)
    add_library(petsc::petsc INTERFACE IMPORTED)
    set_property(TARGET petsc::petsc PROPERTY INTERFACE_LINK_LIBRARIES ${PETSC_LIBRARIES})
    set_property(TARGET petsc::petsc PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${PETSC_INCLUDE_DIRS})
  endif()
endif()

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args(PETSc
  REQUIRED_VARS PETSC_LIBRARIES PETSC_INCLUDE_DIR
  VERSION_VAR PETSC_VERSION)
