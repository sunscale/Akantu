cmake_minimum_required(VERSION 2.6)
project(Akantu_SimTools)

set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries.")

set(AKANTU_SIMTOOLS_SRCS

  common/synchronized_array.cc
  common/parameter_reader.cc
  common/manual_restart.cc

  functions/boundary_functions.cc

  ntn_contact/ntn_base_contact.cc
  ntn_contact/ntn_contact.cc
  ntn_contact/ntrf_contact.cc

  ntn_contact/ntn_base_friction.cc
  ntn_contact/friction_regularisations/ntn_fricreg_no_regularisation.cc
  ntn_contact/friction_regularisations/ntn_fricreg_rubin_ampuero.cc
  ntn_contact/friction_regularisations/ntn_fricreg_simplified_prakash_clifton.cc

  ntn_contact/ntn_initiation_function.cc
  )

set(AKANTU_SIMTOOLS_HDRS
  common/ast_common.hh
  akantu_simtools.hh

  # headers
  common/synchronized_array.hh
  common/parameter_reader.hh
  common/manual_restart.hh

  functions/boundary_functions.hh
  functions/node_filter.hh

  boundary_conditions/force_based_dirichlet.hh
  boundary_conditions/spring_bc.hh

  ntn_contact/ntn_base_contact.hh
  ntn_contact/ntn_contact.hh
  ntn_contact/ntrf_contact.hh
  ntn_contact/ntn_base_friction.hh

  ntn_contact/friction_regularisations/ntn_fricreg_no_regularisation.hh
  ntn_contact/friction_regularisations/ntn_fricreg_rubin_ampuero.hh
  ntn_contact/friction_regularisations/ntn_fricreg_simplified_prakash_clifton.hh

  ntn_contact/friction_laws/ntn_friclaw_coulomb.hh
  ntn_contact/friction_laws/ntn_friclaw_coulomb_tmpl.hh
  ntn_contact/friction_laws/ntn_friclaw_linear_slip_weakening.hh
  ntn_contact/friction_laws/ntn_friclaw_linear_slip_weakening_tmpl.hh

  ntn_contact/ntn_friction.hh
  ntn_contact/ntn_friction_tmpl.hh
  ntn_contact/ntrf_friction.hh
  ntn_contact/ntrf_friction_tmpl.hh

  ntn_contact/ntn_initiation_function.hh

  # inlines
  common/synchronized_array_inline_impl.cc
  ntn_contact/ntn_base_contact_inline_impl.cc
  )


find_package(Akantu REQUIRED)
include_directories(${AKANTU_INCLUDE_DIR})


add_library(akantu_simtools
  ${AKANTU_SIMTOOLS_SRCS}
  )

include_directories(common/)
include_directories(ntn_contact/)
include_directories(ntn_contact/friction_laws/)
include_directories(ntn_contact/friction_regularisations/)
include_directories(ntrf_contact/)
include_directories(ntrf_contact/friction_laws/)

target_link_libraries(akantu_simtools akantu)
#target_link_libraries(akantu_simtools ${AKANTU_LIBRARIES})
set_target_properties(akantu_simtools
  PROPERTIES PUBLIC_HEADER "${AKANTU_SIMTOOLS_HDRS}")

install(TARGETS akantu_simtools
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  PUBLIC_HEADER DESTINATION include/akantu_simtools
  )
