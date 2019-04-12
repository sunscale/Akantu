pybind11_add_module(py11_akantu_test_common test_common.cc)
target_link_libraries(py11_akantu_test_common PRIVATE akantu)
set_property(TARGET py11_akantu_test_common PROPERTY DEBUG_POSTFIX "")

register_test(test_pybind11_common
  SCRIPT test_common.py
  PACKAGE pybind11 python_interface
  PYTHON
  )
