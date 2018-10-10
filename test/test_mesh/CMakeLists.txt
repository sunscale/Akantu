add_mesh(cube_periodic cube_periodic.geo DIM 3 ORDER 1)
add_mesh(square_periodic square_periodic.geo DIM 2 ORDER 1)

register_test(test_mesh_periodic
  SOURCES test_mesh_periodic.cc
  DEPENDS cube_periodic square_periodic
  PACKAGE core)
