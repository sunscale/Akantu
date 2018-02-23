add_mesh(bar_segment_2 bar_segment.geo
  DIM 1 ORDER 1 OUTPUT bar_segment_2.msh)
add_mesh(bar_segment_3 bar_segment.geo
  DIM 1 ORDER 2 OUTPUT bar_segment_3.msh)
add_mesh(bar_triangle_3 bar_triangle.geo
  DIM 2 ORDER 1 OUTPUT bar_triangle_3.msh)
add_mesh(bar_triangle_6 bar_triangle.geo
  DIM 2 ORDER 2 OUTPUT bar_triangle_6.msh)
add_mesh(bar_quadrangle_4 bar_quadrangle.geo
  DIM 2 ORDER 1 OUTPUT bar_quadrangle_4.msh)
add_mesh(bar_quadrangle_8 bar_quadrangle.geo
  DIM 2 ORDER 2 OUTPUT bar_quadrangle_8.msh)
add_mesh(bar_tetrahedron_4 bar_tetrahedron.geo
  DIM 3 ORDER 1 OUTPUT  bar_tetrahedron_4.msh)
add_mesh(bar_tetrahedron_10 bar_tetrahedron.geo
  DIM 3 ORDER 2 OUTPUT bar_tetrahedron_10.msh)
add_mesh(bar_hexahedron_8  bar_hexahedron.geo
  DIM 3 ORDER 1 OUTPUT bar_hexahedron_8.msh)
add_mesh(bar_hexahedron_20 bar_hexahedron.geo
  DIM 3 ORDER 2 OUTPUT bar_hexahedron_20.msh)
add_mesh(bar_pentahedron_6 bar_pentahedron.geo
  DIM 3 ORDER 1 OUTPUT bar_pentahedron_6.msh)
add_mesh(bar_pentahedron_15 bar_pentahedron.geo
  DIM 3 ORDER 2 OUTPUT bar_pentahedron_15.msh)

#===============================================================================
# Meshes
#===============================================================================
set(_patch_tests_meshes)
set(_patch_tests_bar_meshes)
package_get_variable(ET_ELEMENT_TYPES core _element_list)
foreach(_et ${_element_list})
  if(NOT _et STREQUAL _point_1)
    list(APPEND _patch_tests_meshes ${CMAKE_CURRENT_SOURCE_DIR}/${_et}.msh)
    list(APPEND _patch_tests_bar_meshes bar${_et})
  endif()
endforeach()

set(PATCH_TEST_MESHES ${_patch_tests_meshes}
  CACHE INTERNAL "List of mesh files")

set(PATCH_TEST_BAR_MESHES ${_patch_tests_bar_meshes}
  CACHE INTERNAL "List of bar mesh files")
