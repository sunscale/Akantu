#===============================================================================
# Energy tests
#===============================================================================
register_gtest_sources(
  SOURCES test_solid_mechanics_model_linear_elastic_potential_energy.cc
  FILES_TO_COPY test_solid_mechanics_model_linear_elastic_potential_energy_material.dat
  PACKAGE implicit
  )

register_gtest_sources(
  SOURCES test_solid_mechanics_model_kinetic_energy.cc
  FILES_TO_COPY test_solid_mechanics_model_kinetic_energy_material.dat
  PACKAGE implicit
  )

register_gtest_sources(
  SOURCES test_solid_mechanics_model_work_quasistatic.cc
  FILES_TO_COPY test_solid_mechanics_model_work_material.dat
  PACKAGE implicit
  )

register_gtest_sources(
  SOURCES test_solid_mechanics_model_work_dynamics.cc
  FILES_TO_COPY test_solid_mechanics_model_work_material.dat
  PACKAGE core
  )

register_gtest_test(test_solid_mechanics_model_energies
  FILES_TO_COPY ${PATCH_TEST_MESHES}
  DEPENDS ${PATCH_TEST_BAR_MESHES}
  #bar_segment_2 bar_segment_3
  #bar_triangle_3 bar_triangle_6
  #bar_quadrangle_4 bar_quadrangle_8
  #bar_tetrahedron_4 bar_tetrahedron_10
  #bar_hexahedron_8 bar_hexahedron_20
  #bar_pentahedron_6 bar_pentahedron_15
  PARALLEL
  )
