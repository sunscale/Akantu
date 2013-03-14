cmake_minimum_required(VERSION 2.6)
project(Akantu_SimTools)

set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries.")

set(AKANTU_SIMTOOLS_SRCS

  common/syncronized_array.cc

  ntn_contact/ntn_contact.cc
  ntn_contact/ntn_friction.cc
  ntn_contact/friction_laws/ntn_friction_coulomb.cc
  ntn_contact/friction_laws/ntn_friction_linear_slip_weakening.cc
  )

set(AKANTU_SIMTOOLS_HDRS
  common/ast_common.hh
  akantu_simtools.hh

  # headers
  common/syncronized_array.hh
  ntn_contact/ntn_contact.hh
  ntn_contact/ntn_friction.hh
  ntn_contact/friction_laws/ntn_friction_coulomb.hh
  ntn_contact/friction_laws/ntn_friction_linear_slip_weakening.hh

  # inlines
  common/syncronized_array_inline_impl.cc
  )


find_package(Akantu REQUIRED)
include_directories(${AKANTU_INCLUDE_DIR})


add_library(akantu_simtools
  ${AKANTU_SIMTOOLS_SRCS}
  )

include_directories(common/)
include_directories(ntn_contact/)
include_directories(ntn_contact/friction_laws)


target_link_libraries(akantu_simtools ${AKANTU_LIBRARIES})
set_target_properties(akantu_simtools
  PROPERTIES PUBLIC_HEADER "${AKANTU_SIMTOOLS_HDRS}")

install(TARGETS akantu_simtools
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  PUBLIC_HEADER DESTINATION include/akantu_simtools
  )