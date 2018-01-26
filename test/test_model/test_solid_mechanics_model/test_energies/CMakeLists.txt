#===============================================================================
# Meshes
#===============================================================================
set(_meshes)
package_get_variable(ET_ELEMENT_TYPES core _element_list)
foreach(_et ${_element_list})
  if(NOT _et STREQUAL _point_1)
    list(APPEND _meshes ../patch_tests/data/${_et}.msh)
  endif()
endforeach()

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
  FILES_TO_COPY ${_meshes}
  DEPENDS bar_segment_2 bar_segment_3
          bar_triangle_3 bar_triangle_6
	  bar_quadrangle_4 bar_quadrangle_8
	  bar_tetrahedron_4 bar_tetrahedron_10
	  bar_hexahedron_8 bar_hexahedron_20
	  bar_pentahedron_6 bar_pentahedron_15
  PARALLEL
  )
