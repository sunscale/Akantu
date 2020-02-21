add_mesh(python_user_defined_bc_mesh fine_mesh.geo 2 1)

register_example(python_user_defined_bc
  SOURCES python_user_defined_bc.cc
  DEPENDS python_user_defined_bc_mesh
  PYTHON
  FILES_TO_COPY boundary_condition.py material.dat
  )

target_include_directories(python_user_defined_bc
  PRIVATE ${PROJECT_SOURCE_DIR}/python)

target_link_libraries(python_user_defined_bc
  PRIVATE pybind11::embed)
