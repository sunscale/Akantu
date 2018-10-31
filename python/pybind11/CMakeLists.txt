pybind11_add_module(py11_akantu py_akantu.cc)
target_link_libraries(py11_akantu PRIVATE akantu)

set_property(TARGET py11_akantu PROPERTY DEBUG_POSTFIX "")
