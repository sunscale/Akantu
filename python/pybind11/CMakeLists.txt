add_library(pyakantu py_aka_solid_mechanics_model.cc py_aka_material.cc py_aka_boundary_conditions.cc py_aka_common.cc py_aka_mesh.cc)
target_link_libraries(pyakantu PUBLIC  akantu pybind11)
pybind11_add_module(py11_akantu py_akantu.cc)
target_link_libraries(py11_akantu PRIVATE pyakantu)

set_property(TARGET py11_akantu PROPERTY DEBUG_POSTFIX "")
