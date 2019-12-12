cmake_minimum_required(VERSION 3.10)

project(AkantuIterators)

# ------------------------------------------------------------------------------
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(BUILD_SHARED_LIBS ON)
if (NOT AKANTU_ITERATORS_PYTHON_MAJOR_VERSION)
  set(AKANTU_ITERATORS_PYTHON_MAJOR_VERSION 3)
endif()
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

include(AkantuIteratorsTools)
# ------------------------------------------------------------------------------

set(AKANTU_ITERATORS_PUBLIC_HDRS
  aka_compatibilty_with_cpp_standard.hh
  aka_iterators.hh
  aka_static_if.hh
  aka_tuple_tools.hh
  iterators/aka_arange_iterator.hh
  iterators/aka_enumerate_iterator.hh
  iterators/aka_filter_iterator.hh
  iterators/aka_named_zip_iterator.hh
  iterators/aka_transform_iterator.hh
  iterators/aka_zip_iterator.hh
  )

add_library(akantu_iterators INTERFACE)

target_include_directories(akantu_iterators
  INTERFACE $<INSTALL_INTERFACE:include/akantu_iterators>
  )

# small trick for build includes in public
set_property(TARGET akantu_iterators APPEND PROPERTY
  INTERFACE_INCLUDE_DIRECTORIES $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>)

set_property(TARGET akantu_iterators PROPERTY
  INTERFACE_PUBLIC_HEADER ${AKANTU_ITERATORS_PUBLIC_HDRS})

set_target_properties(akantu_iterators
  PROPERTIES
  INTERFACE_CXX_STANDARD 14
  )
# ------------------------------------------------------------------------------
option(AKANTU_ITERATORS_TESTS "Activating tests" OFF)
mark_as_advanced(AKANTU_ITERATORS_TESTS)
if(AKANTU_ITERATORS_TESTS)
  enable_testing()
  add_external_package(GTest)
  add_subdirectory(test)
endif()

# ------------------------------------------------------------------------------
if(NOT AKANTU_ITERATORS_TARGETS_EXPORT)
  set(AKANTU_ITERATORS_TARGETS_EXPORT AkantuIteratorsTargets)
endif()

install(TARGETS akantu_iterators
  EXPORT ${AKANTU_ITERATORS_TARGETS_EXPORT}
  PUBLIC_HEADER DESTINATION include/akantu_iterators/ COMPONENT dev
  )

if("${AKANTU_ITERATORS_TARGETS_EXPORT}" STREQUAL "AkantuIteratorsTargets")
  install(EXPORT AkantuIteratorsTargets DESTINATION share/cmake/${PROJECT_NAME}
    COMPONENT dev)

  #Export for build tree
  export(EXPORT AkantuIteratorsTargets
    FILE "${CMAKE_CURRENT_BINARY_DIR}/AkantuIteratorsTargets.cmake")
  export(PACKAGE AkantuIterators)
endif()
