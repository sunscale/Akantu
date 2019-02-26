add_akantu_test(test_model_solver "Test for the solvers")
add_akantu_test(test_non_local_toolbox "Test of the functionalities in the non-local toolbox")

add_mesh(dof_manager_mesh mesh.geo 3 1)
register_gtest_sources(
  SOURCES test_dof_manager.cc
  PACKAGE core)
register_gtest_test(test_dof_manager DEPENDS dof_manager_mesh)
